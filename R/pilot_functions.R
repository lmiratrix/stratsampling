#' Estimate within-cell residual variances from a coded pilot sample
#'
#' @description
#' Computes, for each treatment-by-stratum cell, the pilot sample size,
#' mean residual, and residual variance, and returns them in the format
#' \code{\link{calc_optimal_allocation}} expects for \code{var_alloc}.
#' Cells with fewer than \code{min_pilot} pilot units (or an inestimable
#' variance) fall back to the pooled within-arm pilot variance, with a
#' warning. Optionally, cell means and variances can be empirical-Bayes
#' shrunk toward pooled, within-arm values.
#'
#' @param dat Data frame
#' @param S Name of (or vector for) the 0/1 pilot indicator. Pilot units
#'   must have non-missing Yobs.
#' @param Z Name of treatment indicator column, or NULL.
#' @param strata Name of stratification column, or NULL for one stratum.
#' @param Yobs,Yhat Names of observed/predicted outcome columns.
#' @param min_pilot Minimum pilot units per cell before falling back to
#'   the pooled variance. Default 2.
#' @param shrink Logical. If \code{FALSE} (default), \code{var.alloc} and
#'   \code{mean.alloc} are just the raw per-cell values (with the
#'   \code{min_pilot} fallback applied to \code{var.alloc}, as before).
#'   If \code{TRUE}, both are empirical-Bayes shrunk, in two stages, with
#'   variances pooled fully across treatment arms and means handled as a
#'   "hybrid": variances are shrunk first, toward a single pooled
#'   variance computed across \emph{all} cells regardless of arm (a
#'   Smyth (2004)-style, df-weighted posterior mean) -- variances are a
#'   harder-to-estimate, second-order concern here, so all cells are
#'   pooled for maximum stability. Those shrunk variances are then used,
#'   in place of the raw variances, to shrink the cell means (a
#'   DerSimonian-Laird-style Fay-Herriot/James-Stein estimator): the
#'   \emph{amount} of shrinkage (the between-cell variance
#'   \code{sigma^2_u}) is likewise estimated once across all cells, but
#'   each cell is shrunk toward its \emph{own arm's} precision-weighted
#'   grand mean, not a single grand mean pooled across arms, since arms
#'   may carry a genuine mean-residual shift. See
#'   \code{scripts/pilot_shrinkage_memo.md} for the derivation and
#'   citations.
#'
#' @return Data frame with one row per cell: \code{.Z}, \code{.strata},
#'   \code{N.zk}, \code{n_pilot.zk}, \code{mean_resid_pilot} (raw),
#'   \code{var_resid_pilot} (raw), \code{var_pooled}, \code{var.alloc},
#'   and \code{mean.alloc}. When \code{shrink = FALSE}, \code{var.alloc}
#'   and \code{mean.alloc} equal the raw per-cell values (subject to the
#'   \code{min_pilot} fallback on \code{var.alloc}); when
#'   \code{shrink = TRUE} they are the shrunk values described above.
#'
#' @export
estimate_pilot_variances <- function( dat, S = "S_pilot",
                                      Z = NULL, strata = "strata",
                                      Yobs = "Yobs", Yhat = "Yhat",
                                      min_pilot = 2, shrink = TRUE ) {
    
    stopifnot( is.logical(shrink), length(shrink) == 1, !is.na(shrink) )

    dat$Yobs = get_col_or_vector( dat, Yobs )
    if ( !is.data.frame(dat) ) {
        dat = as.data.frame(dat)
    }
    if ( !is.null( S ) ) {
        dat$.S_pilot = get_col_or_vector( dat, S )
    } else {
        dat$.S_pilot = 1 # Assume data only pilot
    }
    dat$Yhat = get_col_or_vector( dat, Yhat )
    stopifnot( all( dat$.S_pilot %in% c(0, 1) ) )

    if ( is.null( strata ) ) {
        dat$.strata = "one"
    } else {
        dat$.strata = get_col_or_vector( dat, strata )
    }
    if ( is.null( Z ) ) {
        dat$.Z = 1
    } else {
        dat$.Z = get_col_or_vector( dat, Z )
    }

    if ( any( is.na( dat$Yobs[ dat$.S_pilot == 1 ] ) ) ) {
        stop( "All pilot units must have observed (coded) outcomes" )
    }
    
    dat$resid = dat$Yobs - dat$Yhat

    cell_stats = dat %>%
        group_by( .Z, .strata ) %>%
        dplyr::summarise(
            N.zk = dplyr::n(),
            n_pilot.zk = sum( .S_pilot ),
            mean_resid_pilot = mean( resid[ .S_pilot == 1 ] ),
            var_resid_pilot = var( resid[ .S_pilot == 1 ] ),
            .groups = "drop" )

    # Pooled within-arm pilot variance (df-weighted) for fallback/shrinkage
    cell_stats <- cell_stats %>%
        group_by( .Z ) %>%
        mutate( var_pooled = weighted.mean( var_resid_pilot,
                                            w = pmax( n_pilot.zk - 1, 0 ),
                                            na.rm = TRUE ) ) %>%
        ungroup()

    thin = cell_stats$n_pilot.zk < min_pilot | is.na( cell_stats$var_resid_pilot )
    if ( any( thin ) && !shrink ) {
        warning( sum(thin), " cell(s) have fewer than ", min_pilot,
                 " pilot units; using pooled within-arm variance for those cells." )
    }

    if ( !shrink ) {
        cell_stats <- cell_stats %>%
            mutate( var.alloc = ifelse( thin, var_pooled, var_resid_pilot ),
                    mean.alloc = mean_resid_pilot )
    } else {
        # Stage 1: shrink variances, pooling ALL cells across BOTH
        # arms toward one prior target/df (variances are a
        # harder-to-estimate, second-order concern here). See
        # discussion in scripts/pilot_shrinkage_memo.md.
        var_pooled_all <- weighted.mean( cell_stats$var_resid_pilot,
                                         w = pmax( cell_stats$n_pilot.zk - 1, 0 ),
                                         na.rm = TRUE )
        cell_stats$var.alloc <- .shrink_cell_variances(
            cell_stats$var_resid_pilot, cell_stats$n_pilot.zk, var_pooled_all )

        # Stage 2: shrink means, "hybrid" style -- the amount of shrinkage
        # (sigma^2_u) is estimated once, pooling both arms together for a
        # more stable estimate (as many cells as possible), but each cell
        # is shrunk toward its OWN arm's grand mean, since arms may carry
        # a genuine mean-residual shift. See scripts/pilot_shrinkage_memo.md.
        psi_all <- ifelse( cell_stats$n_pilot.zk > 0,
                           cell_stats$var.alloc / cell_stats$n_pilot.zk, Inf )
        sigma2_u <- .estimate_sigma2_u( cell_stats$mean_resid_pilot, psi_all )

        cell_stats <- cell_stats %>%
            group_by( .Z ) %>%
            mutate( mean.alloc = .shrink_toward_grand_mean(
                mean_resid_pilot,
                ifelse( n_pilot.zk > 0, var.alloc / n_pilot.zk, Inf ),
                sigma2_u ) ) %>%
            ungroup()
    }

    if ( any( !is.finite( cell_stats$var.alloc ) ) ) {
        stop( "Could not form variance parameters for all cells; pilot is too small." )
    }

    cell_stats
}


#' Empirical-Bayes shrinkage of cell residual variances (internal)
#'
#' Smyth (2004)-style shrinkage of per-cell variances toward a common
#' value \code{s0} (in \code{estimate_pilot_variances()}, this is called
#' once across all cells, pooling both treatment arms), treating
#' \code{s2} as \code{s0_k * chisq(d)/d} draws with a scaled-inverse-chi-
#' squared prior on \code{s0_k}, and estimating the prior degrees of
#' freedom by matching moments of \code{log(s2)} to their theoretical
#' variance under that model. See \code{scripts/pilot_shrinkage_memo.md}
#' for the derivation and citation.
#'
#' @param s2 Numeric vector of per-cell raw variances (\code{var_resid_pilot}).
#' @param n Numeric vector of per-cell pilot counts (\code{n_pilot.zk}).
#' @param s0 Scalar (or vector, recycled) pooled/prior variance to shrink
#'   toward.
#' @return Numeric vector of shrunk variances, same length as \code{s2}.
#' @keywords internal
#' @noRd
.shrink_cell_variances <- function( s2, n, s0 ) {

    s0 <- rep_len( s0, length(s2) )
    d <- pmax( n - 1, 0 )
    ok <- d >= 1 & is.finite( s2 ) & s2 > 0

    if ( sum(ok) < 2 ) {
        # Not enough cells to estimate a prior df; pool completely.
        return( s0 )
    }

    v_log  <- var( log( s2[ok] ) )
    avg_tg <- mean( trigamma( d[ok] / 2 ) )
    target <- v_log - avg_tg

    d0 <- if ( !is.finite(target) || target <= 0 ) {
        Inf
    } else {
        root <- tryCatch(
            uniroot( function(x) trigamma(x / 2) - target, c(1e-4, 1e6) )$root,
            error = function(e) NA_real_ )
        if ( is.na(root) ) Inf else 2 * root
    }

    if ( !is.finite(d0) ) {
        return( s0 )
    }

    d_use  <- ifelse( d <= 0 | !is.finite(s2), 0, d )
    s2_use <- ifelse( d <= 0 | !is.finite(s2), 0, s2 )

    ( d0 * s0 + d_use * s2_use ) / ( d0 + d_use )
}


#' Method-of-moments between-cell variance of cell means (internal)
#'
#' DerSimonian-Laird method-of-moments estimate of the between-cell
#' variance component \code{sigma^2_u} for a one-way random-effects model
#' on cell mean residuals, given each cell's sampling variance
#' \code{psi = var_alloc / n}. Intended to be called once across ALL
#' cells (pooling treatment arms together), so that the amount of
#' shrinkage is estimated from as many cells as possible even when the
#' shrinkage target (the grand mean) is kept arm-specific -- see
#' \code{\link{.shrink_toward_grand_mean}}. See
#' \code{scripts/pilot_shrinkage_memo.md} for the derivation and
#' citation.
#'
#' @param ebar Numeric vector of per-cell raw mean residuals.
#' @param psi Numeric vector of per-cell sampling variances (\code{var_alloc/n}).
#' @return A single non-negative numeric value.
#' @keywords internal
#' @noRd
.estimate_sigma2_u <- function( ebar, psi ) {

    ok <- is.finite(ebar) & is.finite(psi) & psi > 0
    if ( sum(ok) < 2 ) {
        return(0)
    }

    w      <- 1 / psi[ok]
    eb     <- ebar[ok]
    ebar_w <- sum(w * eb) / sum(w)
    Q      <- sum( w * (eb - ebar_w)^2 )
    K      <- sum(ok)
    denom  <- sum(w) - sum(w^2) / sum(w)

    if ( denom <= 0 ) {
        return(0)
    }
    max( 0, (Q - (K - 1)) / denom )
}


#' Empirical-Bayes shrinkage of cell mean residuals toward a
#' group-specific grand mean (internal)
#'
#' Fay-Herriot/James-Stein shrinkage of per-cell mean residuals toward
#' the precision-weighted grand mean of \emph{this} group (e.g. one
#' treatment arm) only, given a between-cell variance
#' \code{sigma2_u} that may have been estimated elsewhere (e.g. pooled
#' across groups via \code{\link{.estimate_sigma2_u}}). See
#' \code{scripts/pilot_shrinkage_memo.md} for the derivation and
#' citation.
#'
#' @param ebar Numeric vector of per-cell raw mean residuals (one group).
#' @param psi Numeric vector of per-cell sampling variances (\code{var_alloc/n}).
#' @param sigma2_u Single between-cell variance value to use for every
#'   cell in this call.
#' @return Numeric vector of shrunk means, same length as \code{ebar}.
#' @keywords internal
#' @noRd
.shrink_toward_grand_mean <- function( ebar, psi, sigma2_u ) {

    ok <- is.finite(ebar) & is.finite(psi) & psi > 0
    if ( sum(ok) < 1 ) {
        return( ebar )
    }

    w      <- ifelse( ok, 1 / psi, 0 )
    ebar_w <- sum( w * ifelse(ok, ebar, 0) ) / sum(w)

    A <- ifelse( psi == Inf, 1, psi / (psi + sigma2_u) )
    ifelse( is.na(ebar), ebar_w, (1 - A) * ebar + A * ebar_w )
}



# Playground code ----

if ( FALSE ) {

    # Playground for one_iter_pilot_select() -- generates a toy
    # population + candidate stratifications and calls the function,
    # which will drop into the browser() near the end of its body.

    library( dplyr )   # one_iter_pilot_select() uses bind_rows/mutate/relocate

    set.seed( 42 )
    dat <- generate_population( N = 1000, n_strat = 4,
                                strata_config = "unbalanced",
                                tau = 0.5, bias = "moderate",
                                residual_var = "heterogeneous",
                                R2_target = 0.6 )
    dat$X1 <- rnorm( nrow(dat) )   # a covariate beyond Yhat that is useless, for richer candidates
    table( dat$strata )
    head(dat)

    # We have the predicted outcome, garbage, and the actual stratification to pick from
    candidates <- create_stratification_candidates( dat, vars = c("Yhat", "X1", "strata"), Z = "Z" )
    names(candidates)


    # Change the strata passed over to different stratifications
    # (meaningful vs not) to see how shrinkage differently applies

    dat$S = as.numeric( sample(nrow(dat)) <= 100 )
    yo = dat$Yobs
    yo[ dat$S == 0 ] = NA
    rr <- estimate_pilot_variances( dat = NULL,
                              S=dat$S,
                              Z = dat$Z,
                              strata = candidates[["strata_levels"]],
                              Yobs = yo,
                              Yhat = dat$Yhat,
                              shrink = TRUE )

    rr
    library( tidyverse )
    ggplot( rr, aes( mean_resid_pilot, mean.alloc, col=as.factor(.Z) ) ) +
        geom_point( size = 2 ) +
        coord_fixed() + theme_minimal() +
        geom_abline( intercept = 0, slope = 1 )

    ggplot( rr, aes( var_resid_pilot, var.alloc, col=as.factor(.Z) ) ) +
        geom_point( size = 2 ) +
        coord_fixed() + theme_minimal() +
        geom_abline( intercept = 0, slope = 1 )


    # ---- Compare shrinkage across three very different candidates ----
    #
    # "strata_levels"  -- the TRUE stratification: real between-stratum
    #                     bias/variance heterogeneity, so shrinkage should
    #                     be mild (points stay near the 45-degree line).
    # "X1_quintiles"   -- a stratification on a covariate (X1) with NO
    #                     relationship to Yobs/Yhat, so any apparent
    #                     between-cell spread is pure noise; shrinkage
    #                     should be aggressive (points collapse toward
    #                     the pooled value, off the 45-degree line).
    # "Yhat_X1_2x2"    -- crosses the meaningful Yhat with the useless X1:
    #                     some real signal (from Yhat), but diluted by
    #                     crossing with noise and split into smaller cells;
    #                     shrinkage should land somewhere in between.

    three_candidates <- c( "strata_levels", "X1_quintiles", "Yhat_X1_2x2" )

    rr_three <- purrr::map_dfr( three_candidates, function( cname ) {
        estimate_pilot_variances( dat = NULL,
                                  S = dat$S, Z = dat$Z,
                                  strata = candidates[[ cname ]],
                                  Yobs = yo, Yhat = dat$Yhat,
                                  shrink = TRUE ) %>%
            mutate( candidate = factor( cname, levels = three_candidates ) )
    } )

    p_mean <- ggplot( rr_three, aes( mean_resid_pilot, mean.alloc, col = as.factor(.Z), size = N.zk ) ) +
        facet_wrap( ~ candidate, nrow = 1 ) +
        geom_point() +
        scale_size_continuous( range = c(1, 6), name = "Stratum size (N.zk)" ) +
        geom_abline( intercept = 0, slope = 1, linetype = 2 ) +
        coord_fixed() + theme_minimal() +
        labs( title = "Mean shrinkage: raw vs shrunk cell means",
              x = "raw mean_resid_pilot", y = "shrunk mean.alloc", col = "Z" )
    p_mean

    p_var <- ggplot( rr_three, aes( var_resid_pilot, var.alloc, col = as.factor(.Z), size = N.zk ) ) +
        facet_wrap( ~ candidate, nrow = 1 ) +
        geom_point() +
        scale_size_continuous( range = c(1, 6), name = "Stratum size (N.zk)" ) +
        geom_abline( intercept = 0, slope = 1, linetype = 2 ) +
        coord_fixed() + theme_minimal() +
        labs( title = "Variance shrinkage: raw vs shrunk cell variances",
              x = "raw var_resid_pilot", y = "shrunk var.alloc", col = "Z" )
    p_var

}

