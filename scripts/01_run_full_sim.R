
# Run multifactor simulation looking at benefits of stratification

param_grid <- expand.grid(
  N = 1000,
  coding_budget = seq(0.1,0.9,by=0.1),
  n_strat = c(4),
  strata_config="balanced-exact",
  tau=0,
  bias=c("none","small", "moderate","large","extreme-contrast"),
  residual_var = c("homogeneous","heterogeneous","extreme-contrast"),
  sd_y=c(3),
  R2_target = c(0.4,0.85),
  stringsAsFactors = FALSE
)

cli::cli_inform( "Total simulation scenarios: {nrow(param_grid)}\n" )

#source("scripts/dat_gen_functions.R" )
#source( "scripts/sampling_functions.R" )
#source("scripts/est_functions.R")
library(dplyr)
library(stratsampling)
source("scripts/sim_functions.R")


one_iter = function(params, pilot_frac=NULL) {
  
  dat = generate_population(N=params$N, n_strat=params$n_strat,
                            strata_config=params$strata_config,
                            tau=params$tau,
                            bias=params$bias,
                            residual_var=params$residual_var, sd_y=params$sd_y, R2_target=params$R2_target)

  n_code = round(params$N * params$coding_budget)

  res <- one_iter_all(dat, n_code, include_rand=TRUE, pilot_frac=pilot_frac) %>%
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
  params = param_grid[1,]
  one_iter( params )
  
  one_iter(params, pilot_frac=0.4)
}


library(doParallel)
library(foreach)
library(doRNG)


cores <- min(parallel::detectCores() - 1, nrow(param_grid))
cl <- makeCluster(cores)
registerDoParallel(cl)

n_iter=1000

set.seed(1234)

system.time(
    results <- foreach(i=1:nrow(param_grid),  .packages = c("dplyr","stratsampling")) 
    %dorng% {
        plyr::rdply(n_iter, one_iter(param_grid[i,])) 
        }
)

results = dplyr::bind_rows(results)

saveRDS(results,file="results/sim_results_balanced_exact_strata_all.rds")


## Run again with balanced-approx configuration
set.seed(1234)
param_grid$strata_config="balanced-approx"
system.time(
    results2 <- foreach(i=1:nrow(param_grid), .packages = c("dplyr","stratsampling")) %dorng% {
        plyr::rdply(n_iter, one_iter( param_grid[i,]))  # Run 1000 iterations
    }
)
results2 = dplyr::bind_rows(results2)

saveRDS(results2,file="results/sim_results_balanced_approx_strata_all.rds")


## Run again with unbalanced configuration
set.seed(1234)
param_grid$strata_config="unbalanced"
system.time(
    results3 <- foreach(i=1:nrow(param_grid), .packages = c("dplyr","stratsampling")) %dorng% {
        plyr::rdply(n_iter, one_iter(param_grid[i,]))  # Run 1000 iterations
    }
)
results3 =  dplyr::bind_rows(results3)

saveRDS(results3,file="results/sim_results_unbalanced_strata_all.rds")




# Simulate oracle estimator (i.e., full coding everything) for this data generating process for comparison

params_oracle = expand.grid(
    N=1000,
    n_strat=4,
    strata_config=c("balanced-exact","balanced-approx","unbalanced"),
    tau=0,
    bias=c("none","small", "moderate","large","extreme-contrast"),
    residual_var = c("homogeneous","heterogeneous","extreme-contrast"),
    sd_y=3,
    R2_target=unique(param_grid$R2_target),
    stringsAsFactors = FALSE
)

sim_oracle = function(params) {
    dat = generate_population(N=params$N, n_strat=params$n_strat,
                              strata_config=params$strata_config,
                              tau=params$tau,
                              bias=params$bias,
                              residual_var=params$residual_var, sd_y=params$sd_y, R2_target=params$R2_target)



    est.part <- dat %>% group_by(Z) %>%
        summarise(Ybar=mean(Yobs),
                  var.y=var(Yobs),
                  N=n())

    est=diff(est.part$Ybar)
    var = sum(est.part$var.y/est.part$N)

    res <- data.frame(method="Oracle",
                      est=est,
                      SE = sqrt(var)) %>%
        mutate(N=params$N,
           n_strat=params$n_strat,
           strata_config=params$strata_config,
           tau=params$tau,
           bias=params$bias,
           residual_var=params$residual_var,
           pop.sd = params$sd_y,
           R2_target=params$R2_target)

    res

}

set.seed(1234)
results.oracle <- foreach(i=1:nrow(params_oracle), .combine = 'rbind', .packages = c("dplyr","stratsampling")) %dorng% {
    params <- params_oracle[i,]
    res <- plyr::rdply(n_iter, sim_oracle(params))  # Run 1000 iterations
    return(res)
}
saveRDS(results.oracle,file="results/sim_results_oracle.rds")

stopCluster(cl)

