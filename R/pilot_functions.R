#' Estimate within-cell residual variances from a coded pilot sample
#'
#' @description
#' Computes, for each treatment-by-stratum cell, the pilot sample size,
#' mean residual, and residual variance, and returns them in the format
#' \code{\link{calc_optimal_allocation}} expects for \code{var_alloc}.
#' Cells with fewer than \code{min_pilot} pilot units (or an inestimable
#' variance) fall back to the pooled within-arm pilot variance, with a
#' warning. Optionally, all cell variances can be shrunk toward the
#' pooled value.
#'
#' @param dat Data frame.
#' @param S Name of (or vector for) the 0/1 pilot indicator. Pilot units
#'   must have non-missing Yobs.
#' @param Z Name of treatment indicator column, or NULL.
#' @param strata Name of stratification column, or NULL for one stratum.
#' @param Yobs,Yhat Names of observed/predicted outcome columns.
#' @param min_pilot Minimum pilot units per cell before falling back to
#'   the pooled variance. Default 2.
#' @param shrink Shrinkage weight between 0-1 pulling each cell variance
#'   toward the pooled within-arm variance. Default 0 (no shrinkage).
#'
#' @return Data frame with one row per cell: \code{.Z}, \code{.strata},
#'   \code{N.zk}, \code{n_pilot.zk}, \code{mean_resid_pilot},
#'   \code{var_resid_pilot}, \code{var_pooled}, and \code{var.alloc}.
#'
#' @export
estimate_pilot_variances <- function( dat, S = "S_pilot",
                                      Z = NULL, strata = "strata",
                                      Yobs = "Yobs", Yhat = "Yhat",
                                      min_pilot = 2, shrink = 0 ) {
    
    stopifnot( shrink >= 0, shrink <= 1 )
    
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
    
    cell_stats <- cell_stats %>%
        mutate( var.alloc = (1 - shrink) * var_resid_pilot + shrink * var_pooled,
                var.alloc = ifelse( thin, var_pooled, var.alloc ) )
    
    if ( any( !is.finite( cell_stats$var.alloc ) ) ) {
        stop( "Could not form variance parameters for all cells; pilot is too small." )
    }
    
    cell_stats
}