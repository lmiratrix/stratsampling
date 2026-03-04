

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
                     n_k = round( n_code * p * N_k / sum(N_k) ) ) %>%
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





#' Oracle Optimal (Neyman) Allocation Calculator
#'
#' Computes the optimal sample allocation across strata using Neyman
#' allocation, which allocates more units to strata with higher
#' residual variance. This is an oracle method because it requires
#' knowing the true outcomes (Yobs) for all units in order to compute
#' the actual residuals.
#'
#' @param dat Data frame containing the full population with observed
#'   outcomes and surrogate predictions.
#' @param n_code Integer. Total number of units to sample for coding.
#' @param Z Character string naming the treatment indicator column, or
#'   \code{NULL} to ignore treatment. Default is \code{NULL}.
#' @param strata Character string naming the stratification column, or
#'   \code{NULL} (single stratum). Default is \code{NULL}.
#' @param Yobs Character string or vector. The observed outcome
#'   variable. Default is \code{"Yobs"}.
#' @param Yhat Character string or vector. The surrogate prediction
#'   variable. Default is \code{"Yhat"}.
#' @param return_counts Logical. If \code{TRUE}, returns only the
#'   stratum-level allocation table instead of the full dataset with
#'   sampling indicators. Default is \code{FALSE}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return If \code{return_counts = FALSE} (default), returns the
#'   input data frame augmented with sampling indicator \code{S} and
#'   sampling probability \code{pi}. If \code{return_counts = TRUE},
#'   returns a data frame with the optimal allocation counts per
#'   stratum.
#'
#' @seealso \code{\link{generate_sample_index}} for proportional
#'   allocation, \code{\link{est_ma_strat}} for the stratified
#'   model-assisted estimator.
#'
#' @export
calc_optimal_allocation_oracle <- function( dat, n_code,
                                            Z = NULL,
                                            strata = NULL,
                                            Yobs = "Yobs",
                                            Yhat = "Yhat",
                                            return_counts = FALSE,
                                            seed = NULL ) {
  
  dat$Yobs = get_col_or_vector( dat, Yobs )
  dat$Yhat = get_col_or_vector( dat, Yhat )
  
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
    #dat$.strata = interaction( dat$.strata, dat[[Z]], drop=TRUE )
  } else {
    dat$.Z = 1
  }
  N = nrow(dat)
  
  #dat$Z = sample(rep(0:1, N/2))
  
  #dat = dat %>% mutate(Yobs = (1-Z)*Y0 + Z*Y1,
  #                     Yhat = (1-Z)*Yhat0 + Z*Yhat1) %>%
  #  select(ID, Yobs, Yhat, Z, strata)
  
  
  n.strat = length(unique(dat$.strata))
  
  dat$resid = dat$Yobs - dat$Yhat
  
  
  # Variance of residuals in each treatment x strata
  strat_counts = dat %>%
    group_by(.Z, .strata) %>%
    dplyr::summarise( N.zk = dplyr::n(),
                      var.resid = var(resid),
                      .groups = "drop") %>%
    mutate( weight = N.zk * sqrt(var.resid))
  
  if ( by_tx ) {
    strat_counts <- strat_counts %>%
    group_by(.Z) %>%
    dplyr::mutate( prop = weight/sum(weight),
                   n.zk = round(prop * (n_code/2))) %>%
    ungroup()
  } else {
    strat_counts <- strat_counts %>%
      mutate( prop = weight/sum(weight),
              n.zk = round(prop * n_code))
  }
  
  # Iteratively adjust allocations when n.zk > N.zk
  for(z in unique(strat_counts$.Z)) {
    sc_z <- strat_counts %>% filter(.Z == z)
    n_code_z=n_code / length(unique(strat_counts$.Z))
    
    while (any(sc_z$n.zk > sc_z$N.zk)){
      over_allocated <- sc_z$n.zk >= sc_z$N.zk
      if(!any(over_allocated)) break
      
      sc_z$n.zk[over_allocated] <- sc_z$N.zk[over_allocated]
      n_remaining <- n_code_z - sum(sc_z$n.zk[over_allocated])
      
      sc_remaining <- sc_z %>% filter(!over_allocated)
      total_weight <- sum(sc_remaining$weight)
      
      sc_z$n.zk[!over_allocated] <- round(n_remaining * (sc_remaining$weight / total_weight)) # need to fix rounding issue
      # where sometimes after rounding we get sum(n.zk)>n_code (by 1-2 units)
    }
    
    strat_counts$n.zk[strat_counts$.Z == z] <- sc_z$n.zk
  }
  
  strat_counts <- dplyr::arrange( strat_counts,
                                  .Z, .strata )
  
  if ( return_counts ) {
    return( strat_counts )
  }
  
  dat$S = generate_sample_index( dat,
                                 n_code = strat_counts$n.zk,
                                 strata=".strata",
                                 Z=".Z" )
  
  # Calculate sampling weights for each strata in each treatment arm & merge with full data
  wts = strat_counts %>%
    select(.Z, .strata, N.zk, n.zk) %>%
    mutate(pi = n.zk / N.zk)  # Compute sampling probability
  
  dat.strat = left_join( dat, wts,
                         by = c(".Z", ".strata")) %>%
    dplyr::arrange( .Z, .strata )
  
  dat.strat
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




