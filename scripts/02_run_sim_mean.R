
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


one_iter = function(params) {
    dat = generate_population(N=params$N, n_strat=params$n_strat,
                              strata_config=params$strata_config,
                              tau=params$tau,
                              bias=params$bias,
                              residual_var=params$residual_var, sd_y=params$sd_y, R2_target=params$R2_target)
    
    dat$Z=NULL
    
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
    params = param_grid[1,]
    one_iter( params )
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
    results <- foreach(i=1:nrow(param_grid), .combine = 'rbind', .packages = c("dplyr","stratsampling")) %dorng% {
        params <- param_grid[i,]
        res <- plyr::rdply(n_iter, one_iter(params))  # Run 1000 iterations
        return(res)
    }
)

saveRDS(results,file="results/sim_results_mean_only.rds")






results$method[results$method=="Stratified-Opt"]="Optimal"



analysis_results <- results %>%
    filter(!is.na(est), !is.na(SEhat)) %>%  
    mutate(coding_frac=n_code/N) %>% 
    group_by(N, n_strat, strata_config, tau, bias, residual_var, R2_target, pop.sd, n_code, coding_frac, method) %>%
    mutate(
        CI_lower = est - 2 * SEhat,
        CI_upper = est + 2 * SEhat,
        covers   = if_else(CI_lower < 0 & CI_upper > 0, 1L, 0L)
    ) %>%
    summarise(
        mean_est   = mean(est),
        bias_est   = mean(est - tau),
        coverage   = mean(covers),
        emp_var    = var(est),
        emp_se     = sd(est),
        mean_se    = mean(SE),
        mean_sehat = mean(SEhat),
        mse        = mean((est - tau)^2),
        .groups = "drop"
    )

analysis_results <- analysis_results %>%
    mutate(
        method        = fct_relevel(method, "Subset", "SRS", 
                                    "Stratified", "Optimal"),
        residual_var  = fct_relevel(factor(residual_var), "homogeneous",
                                    "heterogeneous", "extreme-contrast"),
        bias          = fct_relevel(factor(bias), "none", "small", "moderate",
                                    "large", "extreme-contrast")    )


source("scripts/plotting_functions.R")

make_frontier_plot(analysis_results, r2_val = 0.85)
make_frontier_plot(analysis_results, r2_val=0.4)
