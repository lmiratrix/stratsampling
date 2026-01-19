
# Run single scenario to try and debug the optimal allocation with a
# single strata.

library( tidyverse )
library( stratsampling )
source( "scripts/sim_functions.R" )



params <- list(
    N = 1000,
    coding_budget = 0.2,
    n_strat = c(4),
    strata_config="balanced-exact",
    tau=0,
    bias=c("large"),
    residual_var = "heterogeneous",
    sd_y=c(3),
    R2_target = 0.85
)

params



one_iter = function(params) {
    
    dat = generate_population( N=params$N, 
                               n_strat=params$n_strat,
                               strata_config=params$strata_config,
                               tau=params$tau,
                               bias=params$bias,
                               residual_var=params$residual_var,
                               sd_y=params$sd_y, 
                               R2_target=params$R2_target,
                               no_tx = TRUE )
    
    n_code = round(params$N * params$coding_budget)
    
    res <- one_iter_all(dat, n_code,Z=NULL) %>%
        mutate(N=params$N,
               n_code = n_code,
               n_strat=params$n_strat,
               strata_config=params$strata_config,
               tau=params$tau,
               bias=params$bias,
               residual_var=params$residual_var,
               pop.sd = params$sd_y,
               R2_target=params$R2_target)
    
    
    res
    
    
}

if ( FALSE ) {
    # Test one_iter()
    one_iter( params )
}


n_iter=100

set.seed(1234)

res <- map_dfr( 1:n_iter, \(.) one_iter(params),
               .progress=TRUE )



head( res )

res %>%
    group_by( method ) %>%
    summarise( SE_true = sd(est),
               ESE_hat = sqrt( mean(SEhat^2) ) )



# Looking at stratification code more closely ----

# Make sample dataset
dat = generate_population( N=params$N, 
                           n_strat=params$n_strat,
                           strata_config=params$strata_config,
                           tau=params$tau,
                           bias=params$bias,
                           residual_var=params$residual_var,
                           sd_y=params$sd_y, 
                           R2_target=params$R2_target,
                           no_tx = TRUE )
head( dat )
nrow( dat )
table( dat$strata, dat$Z )

n_code = round(params$N * params$coding_budget)
n_code

res <- one_iter_all(dat, n_code,Z=NULL) 
res


# Check the call itself
dat.strat <- calc_optimal_allocation_oracle( dat, n_code,
                                             strata = "strata",
                                             Z = NULL )
head( dat.strat )
table( dat.strat$pi, dat.strat$strata )

ds <- dat.strat %>%
    group_by( strata, pi ) %>%
    summarise( sdYobs = sd( Yobs ),
               sdYhat = sd( Yhat ),
               sdresid = sd( resid ),
               Esamp = sum( pi ),
               n = sum(S),
               N = n(), .groups = "drop" )

# What did we sample?
ds

# Sampled correct number of units:
sum( ds$Esamp )
n_code




# Look at subset code more closely ----

n_groups = length(unique(dat$Z))
n_groups

# SRS select subsample to code in each treatment group
samp.srs = dat %>% group_by(Z) %>% sample_n(n_code/n_groups)
dat.srs = dat %>% 
    mutate(S=ifelse(ID%in%samp.srs$ID,1,0), pi=1) %>%
    as_tibble()
samp.srs
dat.srs

# Estimation with equal sampling probabilities
out.subset=est_subset(dat.srs, S=dat.srs$S, finite_pop=TRUE)

out.subset

out.subset=est_subset(dat.srs, S=dat.srs$S, finite_pop=FALSE)

out.subset

