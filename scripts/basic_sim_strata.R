

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = F)
setwd("../")


library(tidyverse)
library(truncnorm)

N = 1000 # Total sample size
n_code = 200 # fixed coding budget
n.strat = 4 # number of strata


# Make dataset ----
set.seed(123)
strat = sample(rep(1:n.strat, N/n.strat)) # strata indicator
# strat = sample(1:n.strat, size=N, replace=T) # create strata (without forcing equal sample sizes)
Y = round(rtruncnorm(N, a=0, b=6, mean=3, sd=1.25)) # human-coded outcome

b = c(0,0,2,2) # bias in predictions introduced for each strata
v = c(1, 3, 1, 3) # error in predictions for each strata

Yhat = Y + b[strat] + rnorm(N, mean=0, sd=v[strat])
Z = sample(rep(0:1, N/2)) # fake treatment indicator s.t. true tx=0

dat = data.frame(ID=1:N, Yobs=Y, Yhat=Yhat, Z=Z, strata=strat)



# SRS with equal sample sizes in each treatment arm ----


library( stratsampling )
source(here::here( "scripts/sim_functions.R") )



if ( FALSE ) {
  one_iter_srs(dat,n_code=200)
}

res.srs = plyr::rdply( 1000, one_iter_srs(dat,n_code=100), .progress="text")
res.srs %>% dplyr::summarise(mean.est = mean(est), empirical.sd.est = sd(est),
                             SE.est = mean(SE),
                             mean.SEhat.est = mean(SEhat))

# Stratified random sampling within each treatment arm ----

one_iter_strat(dat,n_code=200)
res.strat = plyr::rdply(1000, one_iter_strat(dat,n_code=100), .progress="text")
res.strat %>% summarise(mean.est = mean(est), empirical.sd.est = sd(est),
                        SE.est = mean(SE),
                        mean.SEhat.est = mean(SEhat))


# Larger sim in parallel  ----

cli::cli_h1("Launching parallel simulation" )

library(doParallel)
library(foreach)

cores = min(10, detectCores() - 1)  # Use one less than total cores to avoid overloading
cl = makeCluster(cores)
registerDoParallel(cl)

# Define the range of n_code values
n_code.values <- c(seq(100, 900, by=100),950)

# Run simulations in parallel for each value of n_code
res.out.srs = foreach(nc = n_code.values, .combine = bind_rows, .packages = c("plyr", "dplyr")) %dopar% {
  res <- plyr::rdply(1000, one_iter_srs(dat, n_code=nc))  # Run 1000 iterations
  res$n_code <- nc  # Add n_code identifier
  return(res)
}
res.out.srs$method="SRS"

res.out.strat = foreach(nc = n_code.values, .combine = bind_rows, .packages = c("plyr", "dplyr")) %dopar% {
  res <- plyr::rdply(1000, one_iter_strat(dat, n_code=nc))  # Run 1000 iterations
  res$n_code <- nc  # Add n_code identifier
  return(res)
}
res.out.strat$method="Stratified"

# Stop the parallel cluster
stopCluster(cl)


# Analyze the simulation results ----

res.out=rbind(res.out.srs,res.out.strat)
res_summary = res.out %>% group_by(n_code, method) %>%
  summarise(mean.est = mean(est), empirical.sd.est = sd(est),
            SE.est = mean(SE),
            mean.SEhat.est = mean(SEhat))

RE=res_summary %>% select(n_code, method, empirical.sd.est, SE.est) %>%
  rename(emp.SE=empirical.sd.est, mean.SE=SE.est)  %>%
  tidyr::pivot_wider(names_from=method, values_from=c(emp.SE, mean.SE)) %>%
  mutate(RE_emp_SE = 100*(emp.SE_SRS - emp.SE_Stratified)/emp.SE_SRS,
         RE_mean_SE = 100*(mean.SE_SRS - mean.SE_Stratified)/mean.SE_Stratified)

plot(RE$n_code, RE$RE_mean_SE,type='l',
     ylim=c(0,20),xlim=c(100,1000),
     ylab="% reduction in SE (Strat vs. SRS)",
     xlab="Number of documents sampled/coded")
lines(RE$n_code, RE$RE_emp_SE,type='l',col="blue")
