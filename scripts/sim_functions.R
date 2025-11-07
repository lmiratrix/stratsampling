

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




# One iteration calculating the model-assisted estimator under
# stratified sampling with optimal number sampled in each strata
# (based on variance of residuals--an oracle knowledge quantity)
one_iter_optimal_strat = function(dat, n_code, Z="Z", finite_pop=FALSE){

  require( stratsampling ) # for parallel code
  dat.strat <- calc_optimal_allocation_oracle( dat, n_code,
                                               strata = "strata",
                                               Z = Z )

  # Estimation with sampling probabilities proportional to strata size
  out.strat = est_ma_strat( dat.strat,
                            Yhat = "Yhat", S = "S", Z = Z,
                            strata = "strata" , finite_pop=finite_pop)

  out.strat
}




# Function to compute all three types of model-assisted estimators in one pass
one_iter_all <- function(dat, n_code, Z="Z", finite_pop=FALSE) {
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
  out.strat_opt$method = "Stratified-Opt"

  # Combine results from all three methods
  out_combined = bind_rows(out.subset, out.srs, out.strat_prop, out.strat_opt) %>%
    relocate( method )

  # Return the combined results
  return(out_combined)
}
