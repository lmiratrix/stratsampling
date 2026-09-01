

library(stratsampling)

#' One iteration calculating the simple subset estimator under SRS sampling
one_iter_subset = function(dat, n_code, finite_pop=FALSE){
    
    if (!"Z"%in%names(dat)){
        dat$Z=1
    }
    n_groups = length(unique(dat$Z))
    # SRS select subsample to code in each treatment group
    samp.srs = dat %>% group_by(Z) %>% sample_n(n_code/n_groups)
    dat.srs = dat %>% mutate(S=ifelse(ID%in%samp.srs$ID,1,0), pi=1)


    # Estimation with equal sampling probabilities
    out.subset=est_subset(dat.srs, S=dat.srs$S, finite_pop=finite_pop)

    out.subset
}



#' One iteration calculating the model-assisted estimator under SRS sampling
one_iter_srs = function(dat, n_code, finite_pop=FALSE){
  require( stratsampling ) # for parallel code
  N = nrow(dat)
  #dat$Z = sample(rep(0:1, N/2))
  if (!"Z"%in%names(dat)){
      dat$Z=1
  }
  # SRS select subsample to code in each treatment group
  n_groups = length(unique(dat$Z))
  
  samp.srs = dat %>% group_by(Z) %>% sample_n(n_code/n_groups)
  dat.srs = dat %>% mutate(S=ifelse(ID%in%samp.srs$ID,1,0), pi=1)

  # Estimation with equal sampling probabilities
  out.srs=est_ma_srs(dat.srs, Yhat=dat.srs$Yhat, S=dat.srs$S, finite_pop=finite_pop)

  out.srs
}



# One iteration calculating the model-assisted estimator under stratified sampling
# with number sampled in each strata proportional to strata size

# Estimation with sampling probabilities proportional to strata size
# within treatment arm

one_iter_strat = function(dat, n_code, Z="Z", finite_pop=FALSE){
  require( stratsampling ) # for parallel code

  dat$S = generate_sample_index(dat, n_code,
                                Z=Z, strata="strata" )

  out.strat = est_ma_strat( dat,
                            Yhat = "Yhat", S = "S", Z = Z,
                            strata = "strata" , finite_pop=finite_pop)

  out.strat
}



# One iteration calculating the model-assisted estimator under stratified sampling using RANDOMLY PERMUTED (incorrect) stratification variable
# with number sampled in each strata proportional to strata size

one_iter_random_strat = function(dat, n_code, Z="Z", finite_pop=FALSE){
    require( stratsampling ) # for parallel code
    
    if (!is.null(Z)) {
        dat = dat %>%
            group_by(.data[[Z]]) %>%
            mutate(strata_random = sample(strata)) %>%
            ungroup()
    } else {
        dat$strata_random = sample(dat$strata)
    }
    
    dat$S = generate_sample_index(dat, n_code,
                                  Z=Z, strata="strata_random" )
    
    out.strat = est_ma_strat( dat,
                              Yhat = "Yhat", S = "S", Z = Z,
                              strata = "strata_random" , finite_pop=finite_pop)
    
    out.strat
}



# One iteration calculating the model-assisted estimator under
# stratified sampling with optimal number sampled in each strata
# (based on variance of residuals--an oracle knowledge quantity)
one_iter_optimal_strat = function(dat, n_code, Z="Z", finite_pop=FALSE){

  require( stratsampling ) # for parallel code
  dat.strat <- calc_optimal_allocation_oracle( dat, n_code,
                                               strata = "strata",
                                               Z = Z, distort="none")

  # Estimation with sampling probabilities proportional to strata size
  out.strat = est_ma_strat( dat.strat,
                            Yhat = "Yhat", S = "S", Z = Z,
                            strata = "strata" , finite_pop=finite_pop)

  out.strat
}


# One iteration of the optimal-allocation estimator computed on the CORRECT
# strata but with wrong variance parameters. Sampling is still a valid
# stratified design, so the estimator stays unbiased with valid SEs; the
# cost of misallocation is pure efficiency loss.
one_iter_optimal_wrong = function(dat, n_code, Z="Z", 
                                  distort = c("permute", "reverse","random"), 
                                  finite_pop=FALSE){
    
    require( stratsampling ) # for parallel code
    distort = match.arg(distort)
    
    dat.strat <- calc_optimal_allocation_oracle( dat, n_code,
                                                 strata = "strata",
                                                 Z = Z, distort=distort)
    
    # Estimation with sampling probabilities proportional to strata size
    out.strat = est_ma_strat( dat.strat,
                              Yhat = "Yhat", S = "S", Z = Z,
                              strata = "strata" , finite_pop=finite_pop)
    
    out.strat
}


# One iteration of the pilot-based optimal allocation with the stratification FIXED a priori. 
# A proportional stratified pilot is drawn on the known strata, and estimated within-cell residual variances 
# from the pilot are used for the Neyman allocation of the remaining budget.
#
# pool_pilot = TRUE (default): pilot units are retained in the estimation sample 
# This is classical double sampling with inference conditional on the realized allocation (i.e., allocation depends on pilot variances, but not pilot means)
# Exactly unbiased conditional on the allocation if pilot within-stratum residual variances are independent of residual means
# This is not as clean as the option below, but means the pilot is essentially free since the estimate uses all n_code coded observations and no coded units are discarded
#
# pool_pilot = FALSE: the pilot informs the variance parameters only. 
# In this case, the estimation sample is a fresh optimal-allocation draw of remaining budget (n_code - n_pilot) 
# from the full population, analyzed on its own (exactly unbiased given the pilot sample). 
# This is the cleaner option in terms of conditioning/exact unbiasedness, but costs the size of the pilot sample in terms of efficiency.
#
# pilot_strata: stratification used to DRAW the pilot. Defaults to
#   "strata" (proportional stratified, appropriate when strata are known); set to NULL for an SRS pilot.
one_iter_optimal_pilot = function(dat, n_code, pilot_frac = 0.05,
                                  Z = "Z", finite_pop = FALSE,
                                  pool_pilot = TRUE,
                                  pilot_strata = "strata") {
    require( stratsampling )
    
    N = nrow(dat)
    stopifnot( pilot_frac > 0, pilot_frac < 1 )
    n_pilot = round( pilot_frac * n_code )
    stopifnot( n_pilot >= 1 )
    
    # Wave 1: pilot draw
    dat$S_pilot = generate_sample_index( dat, n_pilot,
                                         Z = Z, strata = pilot_strata )
    
    # Estimate within-cell residual variances (cells of the final strata)
    vt = estimate_pilot_variances( dat, S = "S_pilot",
                                   Z = Z, strata = "strata" )
    
    if ( pool_pilot ) {
        # Wave 2 tops the pilot up to the total budget; estimate on the
        # combined sample.
        dat.strat = calc_optimal_allocation( dat, n_code,
                                             Z = Z, strata = "strata",
                                             var_alloc = vt,
                                             S_pilot = "S_pilot" )
    } else {
        # Fresh stratified sample of the remaining budget from the full
        # population; the pilot is discarded from estimation.
        n_pilot_realized = sum( dat$S_pilot )
        dat.strat = calc_optimal_allocation( dat, n_code - n_pilot_realized,
                                             Z = Z, strata = "strata",
                                             var_alloc = vt )
    }
    
    est_ma_strat( dat.strat,
                  Yhat = "Yhat", S = "S", Z = Z,
                  strata = "strata", finite_pop = finite_pop )
}

# Function to compute all three types of model-assisted estimators in one pass
one_iter_all <- function(dat, n_code, Z="Z", finite_pop=FALSE, 
                         include_rand=TRUE, pilot_frac=NULL) {
  require( stratsampling ) # for parallel code


    #### Subset ESTIMATOR ####

    out.subset=one_iter_subset(dat, n_code, finite_pop=finite_pop)
    out.subset$method <- "Subset"


  #### SRS ESTIMATOR ####

  out.srs=one_iter_srs(dat, n_code, finite_pop=finite_pop)
  out.srs$method <- "SRS"

  #### PROPORTIONAL STRATIFIED ESTIMATOR ####
  out.strat_prop = one_iter_strat(dat, n_code, Z=Z, finite_pop=finite_pop)
  out.strat_prop$method = "Stratified"

  #### OPTIMAL STRATIFIED ESTIMATOR ####
  out.strat_opt = one_iter_optimal_strat(dat, n_code, Z=Z, finite_pop=finite_pop)
  out.strat_opt$method = "Optimal"

  
  #### PROPORTIONAL STRATIFIED ESTIMATOR (WRONG STRATA) ####
  out.strat_random = one_iter_random_strat(dat, n_code, Z=Z, finite_pop=finite_pop)
  out.strat_random$method = "Stratified-Random"
  
  #### OPTIMAL STRATIFIED ESTIMATOR (WRONG VARIANCE PARAMETERS - PERMUTED) ####
  out.strat_opt_perm = one_iter_optimal_wrong(dat, n_code, Z=Z, distort="permute", finite_pop=finite_pop)
  out.strat_opt_perm$method = "Optimal-Permuted"
  
  #### OPTIMAL STRATIFIED ESTIMATOR (WRONG VARIANCE PARAMETERS - REVERSED) ####
  out.strat_opt_rev = one_iter_optimal_wrong(dat, n_code, Z=Z, distort="reverse", finite_pop=finite_pop)
  out.strat_opt_rev$method = "Optimal-Reversed"
  
  
  
  out_combined = bind_rows(out.subset, out.srs, out.strat_prop, out.strat_opt,
                           out.strat_random, out.strat_opt_perm, out.strat_opt_rev) 
  
  if (include_rand){
  #### OPTIMAL STRATIFIED ESTIMATOR (WRONG VARIANCE PARAMETERS - RANDOM) ####
    out.strat_opt_rand = one_iter_optimal_wrong(dat, n_code, Z=Z, distort="random", finite_pop=finite_pop)
    out.strat_opt_rand$method = "Optimal-Random"
  
    out_combined = bind_rows(out_combined, out.strat_opt_rand)
  
  }

  
  if ( !is.null( pilot_frac ) ) {
      out.strat_pilot_pool = one_iter_optimal_pilot(dat, n_code,
                                               pilot_frac = pilot_frac,
                                               Z = Z, finite_pop = finite_pop,
                                               pool_pilot = T)
      out.strat_pilot_pool$method = "Optimal-Pilot-Pooled"

      out.strat_pilot_ind = one_iter_optimal_pilot(dat, n_code,
                                                    pilot_frac = pilot_frac,
                                                    Z = Z, finite_pop = finite_pop,
                                                    pool_pilot = F)
      out.strat_pilot_ind$method = "Optimal-Pilot-Ind"
      out_combined = bind_rows(out_combined, out.strat_pilot_pool,out.strat_pilot_ind)
      out_combined$pilot_frac = pilot_frac
  }
  
  out_combined = out_combined %>%
      relocate( method )

  # Return the combined results
  return(out_combined)
}






## Another possibility: pilot-SELECTED stratification. The stratification is not
# fixed a priori: an SRS pilot is drawn, the pilot-mode diagnostic
# ranks the candidate stratifications, and the winner is used for both
# the proportional and the optimal estimators, with the optimal-pilot estimator determining
# allocation using the pilot variance estimates. The pilot is then discarded from estimation step; 
# wave 2 sample of size n_code-n_pilot is drawn from the full population (see note below for logic).

# Discarding pilot here is needed because selection maximizes the spread of pilot cell residuals, 
# the same quantities the estimator's bias correction averages; drawing wave 2 from the full population
# (rather than the un-piloted remainder) is therefore necessary to make the estimator exactly conditionally unbiased. 

# Any wave-2 unit that was already coded in the pilot reuses its existing code, so no document is coded twice in practice; 
# but in the simulations we still charge the full budget, so simulation results will be slightly conservative under this strategy.
#
# candidates: named list of stratification vectors (build once with
#   create_stratification_candidates(); depends only on covariates).
# min_cell: optional screen requiring every stratum of a candidate to
#   have at least this many units in the full data.
#
# Returns two rows (proportional + optimal), with the chosen
# stratification recorded in strat_choice.
one_iter_pilot_select = function(dat, n_code, candidates,
                                 pilot_select_frac = 0.05,
                                 Z = "Z", finite_pop = FALSE,
                                 min_cell = 0) {
    require( stratsampling )
    
    if (!is.null(Z)){
        n_pilot = round( pilot_select_frac * nrow(dat) /2)*2
    }
    else{
        n_pilot = round( pilot_select_frac * nrow(dat))
        
    }

    n_main  = n_code - n_pilot
    stopifnot( n_pilot >= 1, n_main >= 1 )
    
    # Wave 1: SRS pilot (within arm); no stratification exists yet
    dat$S_pilot = generate_sample_index( dat, n_pilot, Z = Z, strata = NULL )
    
    # Choose the stratification via the pilot-mode diagnostic.
    # Candidates the pilot cannot evaluate (NA pilot spread) sort last
    # and are never chosen.
    cmp = compare_stratifications_precoding( dat, candidates, Z = Z,
                                             Yobs = "Yobs", S = "S_pilot" )
    if ( min_cell > 0 ) {
        cmp_f = cmp %>% filter( min_stratum_size >= min_cell )
        if ( nrow(cmp_f) > 0 ) { cmp = cmp_f }
    }
    
    best = cmp$stratification[[1]]
    if ( is.na( cmp$mean_spread_pilot[[1]] ) ) {
        # No candidate evaluable by this pilot: fall back to the
        # pre-coding proxy ranking
        best = cmp$stratification[[ which.max( cmp$mean_spread ) ]]
    }
    dat$pilot_strata = candidates[[ best ]]
    
    # Pilot variance estimates on the chosen strata (optimal arm only)
    vt = estimate_pilot_variances( dat, S = "S_pilot",
                                   Z = Z, strata = "pilot_strata" )
    
    # Wave 2: fresh draws of the remaining budget from the FULL
    # population; pilot discarded from estimation.
    
    # (a) Proportional allocation (pilot used only for selection)
    dat$S = generate_sample_index( dat, n_main, Z = Z, strata = "pilot_strata" )
    out.prop = est_ma_strat( dat, Yhat = "Yhat", S = "S", Z = Z,
                             strata = "pilot_strata", finite_pop = finite_pop )
    out.prop$method = "Stratified-PilotSelect"
    
    # (b) Optimal allocation with pilot-estimated variances
    dat.opt = calc_optimal_allocation( dat, n_main, Z = Z, strata = "pilot_strata",
                                       var_alloc = vt )
    out.opt = est_ma_strat( dat.opt, Yhat = "Yhat", S = "S", Z = Z,
                            strata = "pilot_strata", finite_pop = finite_pop )
    out.opt$method = "Optimal-PilotSelect"
    
    bind_rows( out.prop, out.opt ) %>%
        mutate( strat_choice = best,
                n_pilot = n_pilot ) %>%
        relocate( method )
}


