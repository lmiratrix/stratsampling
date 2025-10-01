
## Run again with approximate and unbalanced strata


param_grid <- expand.grid(
  N = 1000,
  coding_budget = seq(0.1,0.9,by=0.1),
  n_strat = c(4,5,8),
  strata_config=c("balanced-approx", "unbalanced"),
  tau=0,
  bias=c("none","moderate","large"),
  residual_var = c("homogeneous","heterogeneous","extreme-contrast"),
  stringsAsFactors = FALSE
)


generate_population <- function(N, n_strat, strata_config, tau, bias, residual_var,seed=NULL) {
  
  if (!is.null(seed)){set.seed(seed)}
  
  # Generate strata
  strat = switch(strata_config,
                 "balanced-exact" = sample(rep(1:n_strat, N/n_strat)),
                 "balanced-approx" = sample(1:n_strat, size=N, replace=TRUE, prob=rep(1/n_strat, n_strat)),
                 "unbalanced" = sample(1:n_strat, size = N, prob = runif(n_strat, min=0.1,max=0.8), replace = TRUE))
  
  
  # Specify bias levels per scenario
  b = switch(bias,
             "none" = rep(0, n_strat),
             "moderate" = seq(-1, 1, length=n_strat),
             "large" = seq(-3, 3, length=n_strat),
             "extreme-contrast"=c(-3, rep(0, n_strat-2), 3))
  
  # Specify residual variance per scenario
  var = switch(residual_var,
               "homogeneous" = rep(1, n_strat),
               "heterogeneous" = seq(0.5, 2, length=n_strat),
               "extreme-contrast"=c(0.5, rep(1, n_strat-2), 3))
  
  # Outcome generation
  Y0 = round(rtruncnorm(N, a=0, b=6, mean=3, sd=1.25)) # human-coded outcome
  Y1 = Y0 + tau
  
  # ML predictions with bias
  yhat0 = Y0 + b[strat] + rnorm(N, mean=0, sd=var[strat])
  yhat1 = Y1 + b[strat] + rnorm(N, mean=0, sd=var[strat])
  
  Z = sample(rep(0:1, N/2))
  
  Yobs = ifelse(Z==1,Y1, Y0)
  yhat=ifelse(Z==1, yhat1,yhat0)
  
  #dat = data.frame(ID=1:N, Y0=Y0, Y1=Y1, yhat0=yhat0, yhat1=yhat1, strata=strat)
  dat = data.frame(ID=1:N, Yobs=Yobs, yhat=yhat, Z=Z, strata=strat)
  
  return(dat)
}

source("scripts/est_functions.R")
source("scripts/sim_functions.R")


one_iter = function(params){
  dat = generate_population(N=params$N, n_strat=params$n_strat,
                            strata_config=params$strata_config,
                            tau=params$tau,
                            bias=params$bias,
                            residual_var=params$residual_var)
  
  n_code = round(params$N * params$coding_budget)
  
  res <- one_iter_all(dat, n_code) %>%
    mutate(N=params$N, 
           n_code = n_code,
           n_strat=params$n_strat,
           strata_config=params$strata_config,
           tau=params$tau,
           bias=params$bias,
           residual_var=params$residual_var)
  
  res
  
}



library(doParallel)
library(foreach)

cores <- min(parallel::detectCores() - 1, nrow(param_grid))
cl <- makeCluster(cores)
registerDoParallel(cl)

results = foreach(i=1:nrow(param_grid), .combine = bind_rows, .packages = c("plyr", "dplyr","truncnorm")) %dopar% {
  params <- param_grid[i,]
  res <- plyr::rdply(500, one_iter(params))  # Run 1000 iterations
  return(res)
}

stopCluster(cl)


saveRDS(results,file="results/sim_results_expanded.rds")


results$method[results$method=="Stratified-Opt"]="Optimal"


# Calculate performance metrics
analysis_results = results %>%
  group_by(N, n_strat, strata_config, tau, bias, residual_var, n_code, method) %>%
  summarize(
    # Bias metrics
    mean_est = mean(est),
    bias_est = mean(est - tau),
    
    # Variance metrics
    emp_var = var(est),
    emp_se = sd(est),
    mean_se = mean(SE),
    mean_sehat = mean(SEhat),
    
    # MSE metrics
    mse = mean((est - tau)^2)
  ) %>% ungroup()

analysis_results$bias.var.config = sprintf("Bias: %s, Var: %s", analysis_results$bias, analysis_results$residual_var)

# Empirical SE across configurations
ggplot(analysis_results[analysis_results$n_strat==4,], aes(x = n_code, y = emp_se, color = method)) +
  geom_smooth(se=F) +
  geom_point() +
  facet_grid(bias.var.config ~ strata_config) +
  labs(x = "Proportion of Documents Coded", y = "Empirical SE of impact estimate",
       title = "Efficiency Gains by Bias/Variance Configuration") +
  theme_minimal()



ggplot(analysis_results, aes(x = n_code, y = emp_se, color = method)) +
  geom_smooth(se=F) +
  geom_point() +
  facet_wrap(~bias.var.config,ncol=3) +
  labs(x = "Proportion of Documents Coded", y = "Empirical SE of impact estimate",
       title = "Efficiency Gains by Bias/Variance Configuration") +
  theme_minimal()

# Calculate relative efficiency 
re_analysis <- analysis_results %>%
  select(N, n_strat, strata_config, tau, bias, residual_var, n_code, method, emp_se, mean_se) %>%
  pivot_wider(names_from = method, values_from = c(emp_se, mean_se)) %>%
  mutate(
    re_strat_vs_srs = ((emp_se_SRS^2 - emp_se_Stratified^2)/emp_se_SRS^2) * 100,
    re_optim_vs_srs =((emp_se_SRS^2 - emp_se_Optimal^2)/emp_se_SRS^2) * 100,
    re_optim_vs_strat = ((emp_se_Stratified^2 - emp_se_Optimal^2)/emp_se_Stratified^2) * 100
  )



# 1. Relative Efficiency across different bias/variance configurations
ggplot(re_analysis, aes(x = n_code, y = re_optim_vs_srs)) +
  geom_line() +
  geom_point() +
  facet_grid(bias ~ residual_var,labeller = "label_both") +
  labs(x = "Proportion of Documents Coded", y = "Relative Efficiency (%) vs. SRS",
       title = "Efficiency Gains by Bias/Variance Configuration") +
  theme_minimal()+
  geom_hline(yintercept = 0,lty=2) +
  geom_line(data=re_analysis, aes(x=n_code, y=re_strat_vs_srs),col="blue")+
  geom_point(data=re_analysis, aes(x=n_code, y=re_strat_vs_srs),col="blue")

