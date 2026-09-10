# Largest-remainder apportionment: integer counts summing to `total`,
# proportional to `target`, capped elementwise at `cap`.
apportion_counts <- function(target, total, cap) {
    n <- pmin( floor(target), cap )
    frac <- target - floor(target)
    rem <- total - sum(n)
    while ( rem > 0 ) {
        room <- cap - n
        cand <- ifelse( room > 0, frac, -Inf )
        i <- which.max( cand )
        if ( !is.finite( cand[i] ) ) break   # every cell full
        n[i] <- n[i] + 1
        frac[i] <- frac[i] - 1               # deprioritize once topped up
        rem <- rem - 1
    }
    n
}


#' Generate Stratified or Simple Random Sample Indicators
#'
#' Given a population data frame, generates a binary vector indicating
#' which units are selected for human coding. Supports simple random
#' sampling (SRS), stratified proportional allocation, and custom
#' per-stratum sample sizes.
#'
#' @param dat Data frame containing the population.
#' @param n_code Sample size specification. If a single number, units
#'   are allocated proportionally across strata (and treatment arms if
#'   \code{Z} is specified). If a vector matching the number of
#'   strata-by-treatment cells, specifies exact counts per cell.
#' @param Z Character string naming the treatment indicator column, or
#'   \code{NULL} to ignore treatment assignment. When specified,
#'   sampling is done within substrata defined by crossing strata with
#'   treatment. Default is \code{NULL}.
#' @param strata Character string naming the stratification column, or
#'   \code{NULL} for SRS (all units in one stratum). Default is
#'   \code{NULL}.
#' @param p_tx Numeric between 0 and 1. Proportion of the total
#'   sample allocated to the treatment arm when \code{Z} is specified
#'   and \code{n_code} is a single number. If \code{NULL}, uses the
#'   observed treatment proportion. Default is 0.5.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An integer vector of length \code{nrow(dat)} with values 0
#'   (not sampled) or 1 (sampled).
#'
#' @seealso \code{\link{generate_observed_sample}} for a wrapper that
#'   also masks unsampled outcomes,
#'   \code{\link{calc_sampling_weights}} for computing the implied
#'   sampling probabilities.
#'
#' @export
generate_sample_index <- function( dat,
                                   n_code,
                                   Z = NULL,
                                   strata = NULL,
                                   p_tx = 0.5,
                                   seed=NULL ) {
  
  # First calculate number to sample within arms
  
  by_tx = !is.null(Z)
  
  if ( !is.null( strata ) ) {
    stopifnot( strata %in% colnames(dat) )
    dat$.strata = get_col_or_vector(dat,strata)
  } else {
    dat$.strata = "one"
  }
  if ( by_tx ) {
    stopifnot( Z %in% colnames(dat) )
    dat$.Z = dat[[Z]]
  } else {
    dat$.Z = 1
  }
  n_strat = dplyr::n_distinct(dat$.strata)
  dat$.micro_strata = interaction( dat$.strata, dat$.Z, drop=TRUE )
  
  counts <- dat %>%
    dplyr::group_by( .micro_strata, .strata, .Z ) %>%
    dplyr::summarise( N_k = n(), .groups="drop" )  %>%
    dplyr::arrange( .Z, .strata )
  
  if ( by_tx ) {
    if ( is.null(p_tx) ) {
      p_tx = mean( dat$.Z )
    }
  } else {
    p_tx = 1
  }
  stopifnot( p_tx >= 0 & p_tx <= 1 )
  
  if ( length(n_code) == 1 ) {
  counts <- counts %>%
      dplyr::group_by( .Z ) %>%
      dplyr::mutate( p = ifelse( .Z, p_tx, 1-p_tx ),
                     target = n_code * p * N_k / sum(N_k),
                     n_k = apportion_counts( target,
                                             total = round( n_code * first(p) ),
                                             cap = N_k ) ) %>%
      dplyr::ungroup()
  } else if ( is.list(n_code) || length(n_code) == nrow(counts) ) {
    counts$n_k = n_code
  } else if ( 2*length(n_code) == nrow(counts) ) {
    stopifnot( by_tx )
    counts <- counts %>%
      group_by( .Z ) %>%
      dplyr::mutate( p = ifelse( .Z, p_tx, 1-p_tx ),
                     n_k = round( p * n_code ) ) %>%
      ungroup()
  } else {
    stop( "n_code not correct size for call, given Z and number of strata" )
  }
  
  N = nrow(dat)
  
  if ( any( counts$n_k > counts$N_k ) ) {
    stop(glue::glue( "Cannot sample more than available in a strata: see strata {counts$.strata[ which(counts$n_k > counts$N_k ) ]}" ) )
  }
  
  K = nrow(counts)
  
  
  if (!is.null(seed)){set.seed(seed)}
  
  dat$.ID = 1:N
  samp_ids = unlist( lapply( 1:K, function(k) {
    stratum = counts$.micro_strata[[k]]
    n_k = counts$n_k[[k]]
    ids_k = dat$.ID[ dat$.micro_strata == stratum ]
    ids_k[ sample.int( length(ids_k), size=n_k ) ]
  } ) )
  
  dat$S = 0
  dat$S[ samp_ids ] = 1
  
  dat$S
}



#' Generate Observed Dataset for Simulation
#'
#' Wrapper around \code{\link{generate_sample_index}} that draws a
#' (stratified) sample and sets unsampled outcomes to \code{NA}. This
#' simulates the real-world scenario where only a subsample of units
#' receives human coding.
#'
#' This function is intended for simulation studies only; it requires
#' fully observed outcomes for the entire population.
#'
#' @inheritParams generate_sample_index
#' @param by_tx Logical. If \code{TRUE}, sample separately within
#'   treatment arms. Default is \code{FALSE}.
#'
#' @return The input data frame with two additional/modified columns:
#'   \describe{
#'     \item{S}{Sampling indicator (1 = coded, 0 = not coded)}
#'     \item{Yobs}{Observed outcome; set to \code{NA} for unsampled units}
#'   }
#'
#' @seealso \code{\link{generate_sample_index}},
#'   \code{\link{generate_population}}
#'
#' @export
generate_observed_sample <- function( dat, n_code,
                                      Z = NULL,
                                      strata = NULL,
                                      by_tx = FALSE,
                                      seed=NULL ) {
  
  S = generate_sample_index( dat=dat, n_code = n_code,
                             Z = Z, strata = strata, seed = seed )
  
  dat = dat %>%
    mutate( S = S,
            Yobs = ifelse(S==1, Yobs, NA) )
  
  return(dat)
}





#' Optimal (Neyman) allocation with oracle residual variances or user-supplied variance estimates
#'
#' @description
#' Computes an optimal-allocation stratified sample for the model-assisted
#' estimator. By default this is an oracle: it requires \code{Yobs} for all
#' units and calculates Neyman weights using the actual within-cell residual
#' variances. Alternatively, users can specify \code{var_alloc} (e.g., estimates from a
#' coded pilot) to allocate without full outcome data.
#'
#' If \code{S_pilot} is supplied, those units are treated as already coded:
#' \code{n_code} is interpreted as the TOTAL coding budget (pilot included),
#' the optimal total allocation is computed for that budget, realized pilot
#' counts are subtracted (never un-sampled), and the second wave is drawn
#' only from un-piloted units. The returned \code{S} marks the combined
#' coded sample. Within each cell the combined sample is a valid
#' within-stratum SRS, so \code{\link{est_ma_strat}} can be applied.
#'
#' The \code{distort} argument purposefully corrupts the variance
#' parameters used for allocation (for simulations of misestimated
#' variances); the actual values used are stored in \code{var.alloc}.
#'
#' @param dat Data frame.
#' @param n_code Total number of units for coding, INCLUDING any pilot units
#'   flagged by \code{S_pilot}.
#' @param Z Name of treatment indicator column, or NULL to ignore
#'   treatment. If passed, allocation is computed within arm.
#' @param strata Name of stratification column. If NULL, all units form a
#'   single stratum.
#' @param Yobs,Yhat Names of observed/predicted outcome columns (or
#'   vectors). Only required when \code{var_alloc} is NULL (oracle mode).
#' @param var_alloc Optional data frame of variance parameters with
#'   columns \code{.Z}, \code{.strata}, \code{var.alloc}, one row per
#'   treatment-by-stratum cell (use .Z = 1 when Z is NULL). Output of
#'   \code{\link{estimate_pilot_variances}} works directly.
#' @param S_pilot Optional name of (or vector for) a 0/1 indicator of
#'   already-coded pilot units.
#' @param distort "none" (default), "permute" (shuffle variance parameters
#'   across cells within arm), "random" (generate random variance parameters) or "reverse" (rank-reverse within arm).
#' @param min_n Minimum coded units per cell (capped at cell size).
#'   Default 2 so within-cell variances are estimable. Set 0 for no floor.
#' @param return_counts If TRUE, return the allocation table instead of
#'   the sampled data.
#' @param seed Optional random seed.
#'
#' @return If \code{return_counts} is TRUE, returns a data frame with one row per cell
#'   (sizes, pilot counts, variance parameters, weights, total allocation
#'   \code{n.zk}, and second-wave counts \code{n_new.zk}). Otherwise, returns a data frame
#'   \code{dat} augmented with \code{S} (combined coded sample) and
#'   \code{pi} (realized sampling probabilities).
#'
#' @export
calc_optimal_allocation <- function( dat, n_code,
                                     Z = NULL,
                                     strata = NULL,
                                     Yobs = "Yobs",
                                     Yhat = "Yhat",
                                     var_alloc = NULL,
                                     S_pilot = NULL,
                                     distort = c("none", "permute", "random", "reverse"),
                                     min_n = 2,
                                     return_counts = FALSE,
                                     seed = NULL ) {
    
    distort = match.arg( distort )
    
    
    if ( !is.null(seed) ) { set.seed(seed) }
    
    # make strata depending on strata variable and treatment
    if ( is.null( strata ) ) {
        dat$.strata = "one"
    } else {
        stopifnot( strata %in% colnames(dat) )
        dat$.strata = dat[[strata]]
    }
    
    by_tx = !is.null( Z )
    if ( by_tx ) {
        stopifnot( Z %in% colnames(dat) )
        dat$.Z = dat[[Z]]
    } else {
        dat$.Z = 1
    }
    N = nrow(dat)
    
    # Pilot indicator (all zeros when no pilot passed)
    if ( !is.null( S_pilot ) ) {
        dat$.S_pilot = get_col_or_vector( dat, S_pilot )
        stopifnot( all( dat$.S_pilot %in% c(0, 1) ) )
    } else {
        dat$.S_pilot = 0
    }
    
    # Cell sizes and realized pilot counts
    strat_counts = dat %>%
        group_by( .Z, .strata ) %>%
        dplyr::summarise( N.zk = dplyr::n(),
                          n_pilot.zk = sum( .S_pilot ),
                          .groups = "drop" )
    
    # Variance parameters for the Neyman weights.
    # var.resid holds the oracle truth (NA when unavailable); var.alloc
    # holds the values actually used for allocation.
    if ( !is.null( var_alloc ) ) {
        stopifnot( all( c(".Z", ".strata", "var.alloc") %in% colnames(var_alloc) ) )
        strat_counts <- left_join( strat_counts,
                                   var_alloc %>% dplyr::select( .Z, .strata, var.alloc ),
                                   by = c(".Z", ".strata") )
        if ( any( is.na( strat_counts$var.alloc ) ) ) {
            stop( "var_alloc is missing or NA for some treatment-by-stratum cells" )
        }
        strat_counts$var.resid = NA_real_
    } else {
        dat$Yobs = get_col_or_vector( dat, Yobs )
        dat$Yhat = get_col_or_vector( dat, Yhat )
        dat$resid = dat$Yobs - dat$Yhat
        oracle_vars = dat %>%
            group_by( .Z, .strata ) %>%
            dplyr::summarise( var.resid = var( resid ), .groups = "drop" )
        strat_counts <- left_join( strat_counts, oracle_vars,
                                   by = c(".Z", ".strata") )
        strat_counts$var.alloc = strat_counts$var.resid
    }
    
    # Optional distortion of the allocation variances (within arm)
    if ( distort != "none" ) {
        strat_counts <- strat_counts %>%
            group_by( .Z ) %>%
            dplyr::mutate( var.alloc = switch( distort,
                                               "permute" = sample( var.alloc ),
                                               "reverse" = sort( var.alloc, decreasing = TRUE )[ rank( var.alloc, ties.method = "first" ) ],
                                               "random"  = {
                                                   # Uninformative draws with realistic magnitude:
                                                   # iid lognormal matched to the within-arm mean and spread of the log variance parameters.
                                                   lv = log( pmax( var.alloc, 1e-12 ) )
                                                   s  = sd( lv )
                                                   if ( is.na( s ) ) s = 0
                                                   exp( rnorm( dplyr::n(), mean = mean(lv), sd = s ) )
                                               }
            ) ) %>%
            ungroup()
    }
    
    # Neyman weights and initial (total-budget) allocation
    strat_counts <- strat_counts %>%
        mutate( weight = N.zk * sqrt( var.alloc ) )
    
    if ( by_tx ) {
        strat_counts <- strat_counts %>%
            group_by( .Z ) %>%
            dplyr::mutate( prop = weight / sum(weight),
                           n.zk = round( prop * (n_code / 2) ) ) %>%
            ungroup()
    } else {
        strat_counts <- strat_counts %>%
            mutate( prop = weight / sum(weight),
                    n.zk = round( prop * n_code ) )
    }
    
    # Per-arm cleanup: cap at cell size, apply floors (min_n and realized
    # pilot counts, since pilot units cannot be un-sampled), reconcile so
    # total coded (pilot + second wave) matches the arm budget.
    for ( z in unique( strat_counts$.Z ) ) {
        sc_z <- strat_counts %>% filter( .Z == z )
        n_code_z = round( n_code / length( unique( strat_counts$.Z ) ) )
        
        # Iteratively adjust allocations when n.zk > N.zk
        while ( any( sc_z$n.zk > sc_z$N.zk ) ) {
            over_allocated <- sc_z$n.zk >= sc_z$N.zk
            if ( !any( over_allocated ) ) break
            
            sc_z$n.zk[over_allocated] <- sc_z$N.zk[over_allocated]
            n_remaining <- n_code_z - sum( sc_z$n.zk[over_allocated] )
            
            sc_remaining <- sc_z %>% filter( !over_allocated )
            total_weight <- sum( sc_remaining$weight )
            
            sc_z$n.zk[!over_allocated] <- round( n_remaining * (sc_remaining$weight / total_weight) )
        }
        
        # Floors
        floor_k <- pmax( pmin( min_n, sc_z$N.zk ), sc_z$n_pilot.zk )
        sc_z$n.zk <- pmax( sc_z$n.zk, floor_k )
        
        # Reconcile totals: trim cells with the most slack above their
        # floor; top up (by weight) cells that still have room.
        excess <- sum( sc_z$n.zk ) - n_code_z
        while ( excess > 0 ) {
            slack <- sc_z$n.zk - floor_k
            if ( max( slack ) <= 0 ) break
            i <- which.max( slack )
            sc_z$n.zk[i] <- sc_z$n.zk[i] - 1
            excess <- excess - 1
        }
        if ( excess > 0 ) {
            warning( "Arm ", z, ": budget is below the per-cell floors (pilot + min_n); ",
                     "allocations left at their floors." )
        }
        while ( excess < 0 ) {
            room <- sc_z$N.zk - sc_z$n.zk
            if ( max( room ) <= 0 ) break
            i <- which.max( ifelse( room > 0, sc_z$weight, -Inf ) )
            sc_z$n.zk[i] <- sc_z$n.zk[i] + 1
            excess <- excess + 1
        }
        
        strat_counts$n.zk[ strat_counts$.Z == z ] <- sc_z$n.zk
    }
    
    strat_counts <- strat_counts %>%
        mutate( n_new.zk = n.zk - n_pilot.zk ) %>%
        dplyr::arrange( .Z, .strata )
    
    if ( return_counts ) {
        return( strat_counts )
    }
    
    # Draw the second wave. If S_pilot was supplied (pooled design),
    # .S_pilot flags already-coded units, so we sample only the
    # remainder and the returned S = pilot + top-up. If no pilot was
    # supplied, .S_pilot is all zeros and this is an ordinary fresh
    # stratified draw from the full population (the discard designs
    # rely on this by never passing S_pilot).    
    dat$.ID = 1:N
    S_new = rep( 0, N )
    for ( k in 1:nrow( strat_counts ) ) {
        if ( strat_counts$n_new.zk[[k]] > 0 ) {
            ids_k = dat$.ID[ dat$.Z == strat_counts$.Z[[k]] &
                                 dat$.strata == strat_counts$.strata[[k]] &
                                 dat$.S_pilot == 0 ]
            S_new[ ids_k[ sample.int( length(ids_k), size = strat_counts$n_new.zk[[k]] ) ] ] = 1
        }
    }
    dat$S = as.integer( dat$.S_pilot == 1 | S_new == 1 )
    
    # Realized sampling probabilities
    wts = strat_counts %>%
        dplyr::select( .Z, .strata, N.zk, n_pilot.zk, n.zk ) %>%
        mutate( pi = n.zk / N.zk )
    
    dat.strat = left_join( dat, wts, by = c(".Z", ".strata") ) %>%
        dplyr::arrange( .Z, .strata )
    
    dat.strat
}


#' Oracle Optimal (Neyman) Allocation Calculator (Wrapper)
#'
#' @description
#' Backwards-compatible wrapper for \code{\link{calc_optimal_allocation}}
#' in oracle mode (variance parameters computed from the actual residuals,
#' requiring Yobs for all units, and no pilot).
#'
#' Computes the optimal sample allocation across strata using Neyman
#' allocation, which allocates more units to strata with higher
#' residual variance. This is an oracle method because it requires
#' knowing the true outcomes (Yobs) for all units in order to compute
#' the actual residuals.
#'
#' @inheritParams calc_optimal_allocation
#' @export
calc_optimal_allocation_oracle <- function( dat, n_code,
                                            Z = NULL,
                                            strata = NULL,
                                            Yobs = "Yobs",
                                            Yhat = "Yhat",
                                            distort = c("none", "permute", "random", "reverse"),
                                            min_n = 2,
                                            return_counts = FALSE,
                                            seed = NULL ) {
    
    distort = match.arg( distort )
    
    calc_optimal_allocation( dat = dat, n_code = n_code,
                             Z = Z, strata = strata,
                             Yobs = Yobs, Yhat = Yhat,
                             var_alloc = NULL, S_pilot = NULL,
                             distort = distort, min_n = min_n,
                             return_counts = return_counts, seed = seed )
}


#' Calculate Sampling Weights from a Stratified Sample
#'
#' Given a data frame with a sampling indicator and (optional)
#' stratification, computes the inclusion probability for each unit as
#' the proportion sampled within its stratum (and treatment arm, if
#' specified).
#'
#' @param dat Data frame containing the population.
#' @param Z Character string naming the treatment indicator column, or
#'   \code{NULL} to ignore treatment. Default is \code{NULL}.
#' @param S Character string naming the sampling indicator column
#'   (values 0/1), or a numeric vector of indicators. Default is
#'   \code{"S"}.
#' @param strata Character string naming the stratification column, or
#'   \code{NULL} for unstratified weights. Default is \code{NULL}.
#'
#' @return Numeric vector of sampling probabilities (one per row of
#'   \code{dat}), where each value equals n_k / N_k for the unit's
#'   stratum.
#'
#' @seealso \code{\link{generate_sample_index}},
#'   \code{\link{est_ma_strat}}
#'
#' @export
calc_sampling_weights <- function( dat,
                                   Z = NULL,
                                   S = "S",
                                   strata = NULL ) {
  
  dat$S = get_col_or_vector( dat, S )
  stopifnot( all( dat$S %in% c(0,1) ) )
  
  # make strata depending on strata variable and treatment
  if ( is.null( strata ) ) {
    dat$.strata = "one"
  } else {
    stopifnot( strata %in% colnames(dat) )
    dat$.strata = dat[[strata]]
  }
  
  by_tx = !is.null( Z )
  if ( by_tx ) {
    stopifnot( Z %in% colnames(dat) )
    dat$.Z = dat[[Z]]
  } else {
    dat$.Z = 1
  }
  N = nrow(dat)
  
  wts_prop <- dat %>%
    group_by(.Z, .strata) %>%
    dplyr::summarise(N.zk = n(),
                     n.zk = sum(S==1),
                     pi = n.zk / N.zk,
                     .groups = "drop")
  
  
  wts = left_join(dat, wts_prop, by = c(".Z", ".strata"))
  
  wts$pi
}




