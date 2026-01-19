

#' Given a full population, generate a vector of indicators for who is
#' sampled.
#'
#' @param p_tx Proportion of sample to allocate to treatment arm, if
#'   treatment arm is defined and n_code doesn't specify exact amounts
#'   for each strata.  If NULL, then will allocation proportional to
#'   proportion actually treated in each strata.
#' @param Z If NULL, sample ignoring treatment.  If not null,
#'   sample within substrata defined by interacting strata with
#'   treatment.
#' @param strata Name of stratification to use.  If null will do
#'   SRS (i.e., treat all units as a single strata).
#' @param n_code Number to code.  If single number, will roughly split
#'   across all strata under stratified proportional to the strata
#'   sizes.  If a vector of length strata, will code that many in each
#'   strata (ignoring treatment assignment unless tx col passed). If
#'   tx col passed, then it is best to make sure the passed strata all
#'   are all treatment or all control, although it will interact
#'   strata with treatment indicator, and then order the strata as all
#'   control strata followed by all treatment strata.  So if there is
#'   strata 1,2,3, you would get six strata of 1-0,2-0,3-0, and then
#'   1-1,2-1,3-1. You would then need to pass a vector of length 6 for
#'   n_code (or a single number for total sample size). For SRS, will
#'   sample ignoring treatment assignment unless treatment column
#'   passed, and then either make two strata (a co and tx strata) and
#'   spread n_code across it, or take a vector of length 2 for n_code
#'   for control, then treatment.
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



#' Given a fully observed population, sample and generate observed
#' outcomes, given possibly stratification, etc.  This method would be
#' used for simulation.
#' @param strata Name of stratification to use.  If null will do
#'   SRS.
#' @param n_code Number to code.  If single number, will roughly split
#'   across all strata under stratified proportional to the strata
#'   sizes.  If a vector of length strata, will code that many in each
#'   strata (ignoring treatment assignment unless by_tx=TRUE).  For
#'   srs, will sample ignoring treatment assignment unless by_tx =
#'   TRUE
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





#' This method assumes we have Yobs known for everyone, and thus is an
#' oracle.  It will look at the actual residuals (Yobs - Yhat) to
#' determine optimal allocation of samples across strata (and
#' treatment arms if Z passed) to minimize variance of the
#' model-assisted estimator.
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


#' Given a sampling indicator and stratification, calculate the
#' associated sampling weights.
#'
#' @param S Name of column in dat that is the sampling indicator
#'   (0/1).  Or vector of values.
#'
#' @return Vector of sampling probabilities (weights) for each unit in
#'   dat.
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




if ( FALSE ) {
  
  
  dta = generate_population(N=1000, n_strat=4, oracle=FALSE )
  head( dta )
  hits = dta$strata == 1 & dta$Z == 1
  dta$Yobs[ hits ] = dta$Yobs[ hits ] + rnorm( sum(hits), sd=10 )
  
  calc_optimal_allocation_oracle( dta, n_code = 100, strata="strata", by_tx = TRUE, Z = "Z",
                                  return_counts = TRUE )
  
  dd <- calc_optimal_allocation_oracle( dta, n_code = 100, strata="strata", by_tx = TRUE, Z = "Z",
                                        return_counts = FALSE )
  head( dd )
  table( dd$.strata, dd$S )
  
  sdta = generate_observed_sample( dta,
                                   Z="Z", strata="strata",
                                   n_code=200 )
}
