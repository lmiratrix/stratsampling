

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

# One iteration calculating the model-assisted estimator under stratified sampling
# where the strata are QUANTILE CUTS OF YHAT (a feasible design: uses only
# pre-coding information) instead of the true strata, with the number sampled
# in each stratum proportional to stratum size.
#
# Cuts are made within treatment arm (make_quantile_strata_by_arm(), the same
# construction create_stratification_candidates() uses with cut_by_arm = TRUE),
# so each arm is split into n_groups equal-size strata. With Z = NULL the cuts
# are pooled. n_groups defaults to the number of true strata, so the comparison
# with "Stratified" holds the number of strata fixed.
# If the cuts cannot be formed (e.g., a discrete Yhat with fewer than n_groups
# distinct values), returns an NA row with a warning instead of stopping.
one_iter_yhat_strat = function(dat, n_code, Z="Z", n_groups=NULL, finite_pop=FALSE){
    require( stratsampling ) # for parallel code
    
    if (is.null(n_groups)) {
        n_groups = dplyr::n_distinct(dat$strata)
    }
    z_vec = if (!is.null(Z)) dat[[Z]] else NULL
    
    strata_yhat = stratsampling:::make_quantile_strata_by_arm(dat$Yhat, z_vec,
                                                              n_groups = n_groups,
                                                              var_name = "Yhat")
    if (is.null(strata_yhat)) {
        warning("Stratified-Yhat: could not form ", n_groups,
                " Yhat quantile strata; returning NA.")
        return(data.frame(est = NA_real_, SE = NA_real_, SEhat = NA_real_,
                          SEhat.full = NA_real_, SEhat.samp = NA_real_))
    }
    dat$strata_yhat = strata_yhat
    
    dat$S = generate_sample_index(dat, n_code,
                                  Z=Z, strata="strata_yhat" )
    
    out.strat = est_ma_strat( dat,
                              Yhat = "Yhat", S = "S", Z = Z,
                              strata = "strata_yhat" , finite_pop=finite_pop)
    
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
        n_est = n_code
        dat.strat = calc_optimal_allocation( dat, n_code,
                                             Z = Z, strata = "strata",
                                             var_alloc = vt,
                                             S_pilot = "S_pilot" )
    } else {
        # Fresh stratified sample of the remaining budget from the full
        # population; the pilot is discarded from estimation.
        n_pilot_realized = sum( dat$S_pilot )
        n_est = n_code-n_pilot_realized
        
        dat.strat = calc_optimal_allocation( dat, n_code - n_pilot_realized,
                                             Z = Z, strata = "strata",
                                             var_alloc = vt )
    }
    
    out = est_ma_strat( dat.strat,
                  Yhat = "Yhat", S = "S", Z = Z,
                  strata = "strata", finite_pop = finite_pop ) %>%
        mutate(n_pilot = n_pilot,
               n_est = n_est)
    out
}





# Function to compute all three types of model-assisted estimators, plus bad-design versions and pilot-based variants in one pass
one_iter_all <- function(dat, n_code, Z="Z", finite_pop=FALSE, 
                         include_rand=TRUE, pilot_frac=NULL, pilot_select_frac=NULL, candidates=NULL,
                         n_groups_yhat=NULL,...) {
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
    
    #### PROPORTIONAL STRATIFIED ESTIMATOR (YHAT-QUANTILE STRATA) ####
    out.strat_yhat = one_iter_yhat_strat(dat, n_code, Z=Z, n_groups=n_groups_yhat, finite_pop=finite_pop)
    out.strat_yhat$method = "Stratified-Yhat"
    
    #### PROPORTIONAL STRATIFIED ESTIMATOR (WRONG STRATA) ####
    out.strat_random = one_iter_random_strat(dat, n_code, Z=Z, finite_pop=finite_pop)
    out.strat_random$method = "Stratified-Random"
    
    #### OPTIMAL STRATIFIED ESTIMATOR (WRONG VARIANCE PARAMETERS - PERMUTED) ####
    out.strat_opt_perm = one_iter_optimal_wrong(dat, n_code, Z=Z, distort="permute", finite_pop=finite_pop)
    out.strat_opt_perm$method = "Optimal-Permuted"
    
    #### OPTIMAL STRATIFIED ESTIMATOR (WRONG VARIANCE PARAMETERS - REVERSED) ####
    out.strat_opt_rev = one_iter_optimal_wrong(dat, n_code, Z=Z, distort="reverse", finite_pop=finite_pop)
    out.strat_opt_rev$method = "Optimal-Reversed"
    
    
    
    out_combined = bind_rows(out.subset, out.srs, out.strat_prop, out.strat_opt,out.strat_yhat,
                             out.strat_random, out.strat_opt_perm, out.strat_opt_rev) 
    
    if (include_rand){
        #### OPTIMAL STRATIFIED ESTIMATOR (WRONG VARIANCE PARAMETERS - RANDOM) ####
        out.strat_opt_rand = one_iter_optimal_wrong(dat, n_code, Z=Z, distort="random", finite_pop=finite_pop)
        out.strat_opt_rand$method = "Optimal-Random"
        
        out_combined = bind_rows(out_combined, out.strat_opt_rand)
        
    }
    
    if (is.null(candidates)){candidates=NULL}
    
    if ( !is.null( pilot_frac ) ) {
        out.strat_pilot_pool = one_iter_optimal_pilot(dat, n_code,
                                                      pilot_frac = pilot_frac,
                                                      Z = Z, finite_pop = finite_pop,
                                                      pool_pilot = T)
        out.strat_pilot_pool$method = "Pilot-Allocated-Pooled"
        
        out.strat_pilot_ind = one_iter_optimal_pilot(dat, n_code,
                                                     pilot_frac = pilot_frac,
                                                     Z = Z, finite_pop = finite_pop,
                                                     pool_pilot = F)
        out.strat_pilot_ind$method = "Pilot-Allocated-Discard"
        

        if (is.null(pilot_select_frac)){pilot_select_frac=pilot_frac}
        
        out.pilot_select = suppressWarnings(one_iter_pilot_select(dat, n_code, candidates=candidates, pilot_select_frac=pilot_select_frac,
                                                                  Z=Z, finite_pop=finite_pop,...))
        
        out_combined = bind_rows(out_combined, out.strat_pilot_pool,out.strat_pilot_ind, out.pilot_select)
        out_combined$pilot_frac = pilot_frac
    }
    
    out_combined = out_combined %>%
        relocate( method )
    
    # Return the combined results
    return(out_combined)
}



#' Precision-weighted combination of two independent treatment-effect estimates
#'
#' Combines two estimator outputs (e.g. a pilot-wave estimate and a
#' second-wave estimate) of the SAME estimand into one pooled estimate,
#' using inverse-variance ("precision") weighting: the estimate with the
#' smaller \code{SEhat^2} gets the larger weight. \code{est}, \code{SEhat},
#' \code{SEhat.full}, \code{SEhat.samp}, and \code{SE} are all combined as
#' the variance of a weighted sum of independent estimators
#' (\code{sqrt(w1^2*var1 + w2^2*var2)}); \code{SEhat.full}/\code{SEhat.samp}
#' reuse the exact weights derived from \code{SEhat} (not their own
#' variances) so that the additive decomposition
#' \code{SEhat^2 = SEhat.full^2 + SEhat.samp^2} is preserved in the pooled
#' result. \code{SE} (the oracle SE) is combined with its own,
#' separately-estimated weight. \code{S2.0}/\code{S2.1} are combined as a
#' plain precision-weighted average (they are already variance-scale, and
#' reuse the \code{SEhat}-derived weight as a practical approximation
#' rather than their own precision).
#'
#' @param est1,est2 One-row estimator outputs (as returned by
#'   e.g. \code{\link{est_ma_srs}}/\code{\link{est_ma_strat}}), each with
#'   \code{est}, \code{SEhat}, \code{SEhat.full}, \code{SEhat.samp},
#'   \code{S2.0}, \code{S2.1}, \code{SE}, and \code{method}.
#' @return A one-row list/data frame in the same format as \code{est1},
#'   with all fields replaced by their pooled values and
#'   \code{method} set to \code{"<est1$method>+<est2$method>"}.
#'
#' @details
#' This assumes \code{est1} and \code{est2} are independent. If the two
#' waves share coded units (e.g. a pilot pooled into the estimation
#' sample), that independence is only approximate, and the pooled SEs
#' here will be mildly anti-conservative (too small); see the discussion
#' above \code{\link{one_iter_pilot_select}}.
combine_estimators <- function( est1, est2 ) {

  var1 = est1$SEhat^2
  var2 = est2$SEhat^2
  
  prec = var2 / (var1 + var2)
  
  est1$est = prec * est1$est + (1-prec)*est2$est
  est1$SEhat = sqrt( prec^2 * var1 + (1-prec)^2 * var2 )
  est1$SEhat.full = sqrt( prec^2 * est1$SEhat.full^2 + (1-prec)^2 * est2$SEhat.full^2 )
  est1$SEhat.samp = sqrt( prec^2 * est1$SEhat.samp^2 + (1-prec)^2 * est2$SEhat.samp^2 )
  est1$S2.0 = prec*est1$S2.0 + (1-prec)*est2$S2.0
  est1$S2.1 = prec*est1$S2.1 + (1-prec)*est2$S2.1

  var1T = est1$SE^2
  var2T = est2$SE^2
  prec = var2T / (var1T + var2T)
  est1$SE = sqrt( prec^2 * var1T + (1-prec)^2 * var2T )
  
  est1$method = paste0( est1$method, "+", est2$method )
  return( est1 )
}





## Pilot-Selected stratification
#
# Another possibility: pilot-SELECTED stratification. The stratification is not
# fixed a priori: an SRS pilot is drawn, the pilot-mode diagnostic
# ranks the candidate stratifications, and the winner is used for both
# the proportional and the optimal estimators, with the optimal-pilot estimator determining
# allocation using the pilot variance estimates.
#
# For second step there are two options:
# (1) The pilot is then discarded from estimation step; 
# wave 2 sample of size n_code-n_pilot is drawn from the full population (see note below for logic).
# (2) The pilot estimate and the main estimate are combined with a precision weighting.

# Arguably, Discarding pilot here is needed because selection maximizes the spread of pilot cell residuals, 
# the same quantities the estimator's bias correction averages; drawing wave 2 from the full population
# (rather than the un-piloted remainder) is therefore necessary to make the estimator exactly conditionally unbiased. 
#
# Counterargument: The pilot is an unbiased estimate of the effect.  The second sample is an independent draw (possibly not optimally allocated and possibly on the wrong strata) that is also unbiased.  One wrinkle: If we disallow reusing units, then we have correlation.  For now we resample and possibly have overlap.  Second wrinkle: This means the two impact estimates are not actually independent. Hopefully this correlation term is small and thus ignorable.  Simulation will show.  Future work here.

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

#' @param pilot_select_frac Fraction of _FULL SAMPLE_ for pilot.
one_iter_pilot_select = function(dat, n_code, candidates,
                                 pilot_select_frac = 0.05,
                                 Z = "Z", Yhat = "Yhat", Yobs = "Yobs",
                                 finite_pop = FALSE,
                                 min_cell = 0,
                                 min_pilot_cell=0,
                                 combine_estimators=FALSE,
                                 add_selection_covs = FALSE) {
    require( stratsampling )
    
    if (!is.null(Z)){
        n_pilot = round( pilot_select_frac * n_code / 2 )*2
    }
    else{
        n_pilot = round( pilot_select_frac * n_code)
        
    }
  
    n_main  = n_code - n_pilot
    stopifnot( n_pilot >= 1, n_main >= 1 )
    
    # Wave 1: SRS pilot (within arm); no stratification exists yet
    dat$S_pilot = generate_sample_index( dat, n_pilot, Z = Z, strata = NULL )
    
    
    # Choose the stratification via the pilot-mode diagnostic.
    # Candidates the pilot cannot evaluate (NA pilot spread) sort last
    # and are never chosen.
    cust=NULL
    
    if (add_selection_covs){
        cust = list("strata_true" = dat$strata,
                    "strata_noisy" = dat$strata_noisy,
                    "strata_random" = dat$strata_random)
    }
    if (is.null(candidates)){
        candidates = create_stratification_candidates(dat,vars="Yhat",Z=Z,
                                                      custom=cust,cross_custom = F,
                                                      verbose=F)
    }
    cmp = compare_stratifications_precoding( dat, candidates, Z = Z,
                                             Yobs = Yobs, S = "S_pilot",shrink=TRUE )
    # What proxy mode (Yhat only, no outcomes) would have picked from the same candidates
    proxy_best = cmp$stratification[[ which.max( cmp$mean_spread ) ]]
   
    # Possibly subset to strata that have at least a given size of units?
    # TODO/QUESTION: Not following this logic
    if ( min_cell > 0 ) {
        cmp_f = cmp %>% filter( min_stratum_size >= min_cell )
        if ( nrow(cmp_f) > 0 ) { 
          cmp = cmp_f 
        }
    }
    
    # (b) after the existing `if ( min_cell > 0 ) {...}` block
    if ( min_pilot_cell > 0 && "min_pilot_cell" %in% names(cmp) ) {
        cmp_f = cmp %>% filter( !is.na(min_pilot_cell), min_pilot_cell >= !!min_pilot_cell )
        if ( nrow(cmp_f) > 0 ) cmp = cmp_f
    }
    
    # Select stratification
    best = cmp$stratification[[1]]
    if ( is.na( cmp$mean_spread_pilot[[1]] ) ) {
        # No candidate evaluable by this pilot: fall back to the
        # pre-coding proxy ranking
        best = cmp$stratification[[ which.max( cmp$mean_spread ) ]]
    }
    dat$pilot_strata = candidates[[ best ]]
    
    # Get Wave 1 (Pilot) estimates
    out.wave1 = est_ma_srs( dat,
                            Yhat = Yhat, Yobs = Yobs,
                            S = "S_pilot", Z = Z,
                            finite_pop=finite_pop )
    out.wave1$method = "Pilot"
    
    # Wave 2: fresh draws of the remaining budget from the FULL
    # population; pilot discarded from estimation.
    
    # (a) Proportional allocation (pilot used only for selection)
    dat$S = generate_sample_index( dat, n_main, Z = Z, strata = "pilot_strata" )
    out.prop = est_ma_strat( dat, Yhat = Yhat, S = "S", Z = Z,
                             strata = "pilot_strata", finite_pop = finite_pop )
    out.prop$method = "Pilot-Selected-Prop"
    
    # (b) Optimal allocation with pilot-estimated variances
    vt = estimate_pilot_variances( dat, S = "S_pilot",
                                   Z = Z, strata = "pilot_strata", 
                                   shrink = TRUE )
    dat.opt = calc_optimal_allocation( dat, n_main, Z = Z, strata = "pilot_strata",
                                       var_alloc = vt )
    out.opt = est_ma_strat( dat.opt, Yhat = Yhat, S = "S", Z = Z,
                            strata = "pilot_strata", finite_pop = finite_pop )
    out.opt$method = "Pilot-Selected-Optimal"
    
    
    # Pooling Wave (1) and Wave (2)
    if (combine_estimators){
    out.prop.pool = combine_estimators( out.wave1, out.prop )
    out.opt.pool = combine_estimators( out.wave1, out.opt )
    
    bind_rows( out.prop, out.opt, out.prop.pool, out.opt.pool ) %>%
        mutate( strat_choice = best,       
                proxy_choice = proxy_best,
                n_pilot = n_pilot ) %>%
        relocate( method )
    }
    else{
        bind_rows( out.prop, out.opt ) %>%
            mutate( strat_choice = best,
                    proxy_choice = proxy_best,
                    n_pilot = n_pilot,
                    n_est = n_main) %>%
            relocate( method )
    }
}


if ( FALSE ) {

    # Playground for one_iter_pilot_select() -- generates a toy
    # population + candidate stratifications and calls the function,
    # which will drop into the browser() near the end of its body.

    library( dplyr )   # one_iter_pilot_select() uses bind_rows/mutate/relocate

    set.seed( 42 )
    dat <- generate_population( N = 1000, n_strat = 4,
                                strata_config = "balanced-exact",
                                tau = 0.5, bias = "moderate",
                                residual_var = "heterogeneous",
                                R2_target = 0.6 )
    dat$X1 <- rnorm( nrow(dat) )   # a covariate beyond Yhat, for richer candidates

    candidates <- create_stratification_candidates( dat, vars = c("Yhat", "X1"), Z = "Z" )
    names(candidates)
    table( candidates[[7]] )
    result <- one_iter_pilot_select( dat, n_code = 200, candidates = candidates,
                                     pilot_select_frac = 0.05,
                                     Z = "Z", Yhat = "Yhat", Yobs = "Yobs" )
}




