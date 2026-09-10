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
#' @param dat Data frame.
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
#'   If \code{TRUE}, both are empirical-Bayes shrunk within each
#'   treatment arm: variances are shrunk first (a Smyth (2004)-style,
#'   df-weighted posterior mean toward the pooled within-arm variance),
#'   and those shrunk variances are then used, in place of the raw
#'   variances, to shrink the cell means (a DerSimonian-Laird-style
#'   Fay-Herriot/James-Stein estimator toward the precision-weighted
#'   grand mean). See \code{scripts/pilot_shrinkage_memo.md} for the
#'   derivation and citations.
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
                                      min_pilot = 2, shrink = FALSE ) {

    stopifnot( is.logical(shrink), length(shrink) == 1, !is.na(shrink) )

    dat$.S_pilot = get_col_or_vector( dat, S )
    dat$Yobs = get_col_or_vector( dat, Yobs )
    dat$Yhat = get_col_or_vector( dat, Yhat )
    stopifnot( all( dat$.S_pilot %in% c(0, 1) ) )

    if ( is.null( strata ) ) {
        dat$.strata = "one"
    } else {
        stopifnot( strata %in% colnames(dat) )
        dat$.strata = dat[[strata]]
    }
    if ( is.null( Z ) ) {
        dat$.Z = 1
    } else {
        stopifnot( Z %in% colnames(dat) )
        dat$.Z = dat[[Z]]
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
    if ( any( thin ) ) {
        warning( sum(thin), " cell(s) have fewer than ", min_pilot,
                 " pilot units; using pooled within-arm variance for those cells." )
    }

    if ( !shrink ) {
        cell_stats <- cell_stats %>%
            mutate( var.alloc = ifelse( thin, var_pooled, var_resid_pilot ),
                    mean.alloc = mean_resid_pilot )
    } else {
        # Two-stage empirical Bayes: shrink variances first, then use the
        # shrunk variances (not the raw ones) as the sampling variances
        # for shrinking the means. See scripts/pilot_shrinkage_memo.md.
        cell_stats <- cell_stats %>%
            group_by( .Z ) %>%
            mutate( var.alloc  = .shrink_cell_variances( var_resid_pilot, n_pilot.zk, var_pooled ),
                    mean.alloc = .shrink_cell_means( mean_resid_pilot, var.alloc, n_pilot.zk ) ) %>%
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
#' within-arm value \code{s0}, treating \code{s2} as
#' \code{s0_k * chisq(d)/d} draws with a scaled-inverse-chi-squared prior
#' on \code{s0_k}, and estimating the prior degrees of freedom by
#' matching moments of \code{log(s2)} to their theoretical variance
#' under that model. See \code{scripts/pilot_shrinkage_memo.md} for the
#' derivation and citation.
#'
#' @param s2 Numeric vector of per-cell raw variances (\code{var_resid_pilot}).
#' @param n Numeric vector of per-cell pilot counts (\code{n_pilot.zk}).
#' @param s0 Numeric vector (constant within group) of the pooled/prior
#'   variance to shrink toward (\code{var_pooled}).
#' @return Numeric vector of shrunk variances, same length as \code{s2}.
#' @keywords internal
.shrink_cell_variances <- function( s2, n, s0 ) {

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


#' Empirical-Bayes shrinkage of cell mean residuals (internal)
#'
#' DerSimonian-Laird-style Fay-Herriot/James-Stein shrinkage of per-cell
#' mean residuals toward the precision-weighted grand mean, using
#' \code{var_alloc / n} as each cell's sampling variance. Pass the
#' shrunk variances (not the raw ones) for the two-stage estimator. See
#' \code{scripts/pilot_shrinkage_memo.md} for the derivation and
#' citation.
#'
#' @param ebar Numeric vector of per-cell raw mean residuals.
#' @param var_alloc Numeric vector of per-cell variances to use for the
#'   sampling variance of \code{ebar} (\code{var_alloc / n}).
#' @param n Numeric vector of per-cell pilot counts.
#' @return Numeric vector of shrunk means, same length as \code{ebar}.
#' @keywords internal
.shrink_cell_means <- function( ebar, var_alloc, n ) {

    psi <- ifelse( n > 0, var_alloc / n, Inf )
    ok  <- is.finite(ebar) & is.finite(psi) & psi > 0

    if ( sum(ok) < 2 ) {
        return( ebar )
    }

    w      <- ifelse( ok, 1 / psi, 0 )
    ebar_w <- sum( w * ifelse(ok, ebar, 0) ) / sum(w)
    Q      <- sum( w * ( ifelse(ok, ebar, ebar_w) - ebar_w )^2 )
    K      <- sum(ok)
    denom  <- sum(w) - sum(w^2) / sum(w)

    sigma2_u <- if ( denom > 0 ) max( 0, (Q - (K - 1)) / denom ) else 0

    A <- ifelse( psi == Inf, 1, psi / (psi + sigma2_u) )
    ifelse( is.na(ebar), ebar_w, (1 - A) * ebar + A * ebar_w )
}


