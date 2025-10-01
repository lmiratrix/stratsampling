setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = F)
setwd("../")


library(tidyverse)
library(truncnorm)

N = 1000 # Total sample size
n.code = 200 # fixed coding budget
n.strat = 4 # number of strata


#set.seed(123)
strat = sample(rep(1:n.strat, N/n.strat)) # strata indicator
# strat = sample(1:n.strat, size=N, replace=T) # create strata (without forcing equal sample sizes)
Y = round(rtruncnorm(N, a=0, b=6, mean=3, sd=1.25)) # human-coded outcome

b = c(0,0,0,0) # bias in predictions introduced for each strata - try none, small, large x constant var, hetero var
v = c(2,2,2,2) # error in predictions for each strata

yhat = Y + b[strat] + rnorm(N, mean=0, sd=v[strat])
Z = sample(rep(0:1, N/2)) # fake treatment indicator s.t. true tx=0

dat = data.frame(ID=1:N, Yobs=Y, yhat=yhat, Z=Z, strata=strat)

## SRS with equal sample sizes in each treatment arm

source("scripts/est_functions.R")
source("scripts/sim_functions.R")


if (FALSE){
one_iter_optimal_strat(dat,n.code=200)
res.opt.strat = plyr::rdply(1000, one_iter_optimal_strat(dat,n.code=100), .progress="text")
res.opt.strat %>% summarise(mean.est = mean(est), empirical.sd.est = sd(est),
                        SE.est = mean(SE),
                        mean.SEhat.est = mean(SEhat))
}

## Do things in parallel to speed it up

library(doParallel)
library(foreach)

cores = min(15, detectCores() - 1)  # Use one less than total cores to avoid overloading
cl = makeCluster(cores)
registerDoParallel(cl)

# Define the range of n.code values
n.code.values <- c(seq(100, 900, by=100))

# Run simulations in parallel for each value of n.code
res.out.srs = foreach(nc = n.code.values, .combine = bind_rows, .packages = c("plyr", "dplyr")) %dopar% {
  res <- plyr::rdply(1000, one_iter_srs(dat, n.code=nc))  # Run 1000 iterations
  res$n.code <- nc  # Add n.code identifier
  return(res)
}
res.out.srs$method="SRS"

res.out.strat = foreach(nc = n.code.values, .combine = bind_rows, .packages = c("plyr", "dplyr")) %dopar% {
  res <- plyr::rdply(1000, one_iter_strat(dat, n.code=nc))  # Run 1000 iterations
  res$n.code <- nc  # Add n.code identifier
  return(res)
}
res.out.strat$method="Stratified"


res.out.opt = foreach(nc = n.code.values, .combine = bind_rows, .packages = c("plyr", "dplyr")) %dopar% {
  res <- plyr::rdply(1000, one_iter_optimal_strat(dat, n.code=nc))  # Run 1000 iterations
  res$n.code <- nc  # Add n.code identifier
  return(res)
}
res.out.opt$method="StratOpt"


# Stop the parallel cluster
stopCluster(cl)


res.out=rbind(res.out.srs,res.out.strat,res.out.opt)
res_summary = res.out %>% group_by(n.code, method) %>%
  summarise(mean.est = mean(est), empirical.sd.est = sd(est),
            SE.est = mean(SE), 
            mean.SEhat.est = mean(SEhat))


rs = res_summary
rs$method[rs$method=="StratOpt"]="Stratified (Optimal n)"
ggplot(rs,aes(x=n.code,y=empirical.sd.est,col=method))+
  geom_smooth(se=F)+geom_point()+
  scale_x_continuous(limits=c(0,950),breaks=seq(0,1000,by=100))+
  scale_y_continuous(limits=c(0,1.2),breaks=seq(0,1,by=0.1))+
  theme_minimal()+
  labs(x="Total number of documents coded", y="SE of estimated treatment effect",
       title=sprintf("Bias_k=(%s), Var_k=(%s)", paste0(b,collapse=","),paste0(v,collapse=",")))+
  theme(legend.position=c(0.8, 0.6))


RE=res_summary %>% select(n.code, method, empirical.sd.est, SE.est) %>% 
  rename(emp.SE=empirical.sd.est, mean.SE=SE.est)  %>% 
  tidyr::pivot_wider(names_from=method, values_from=c(emp.SE, mean.SE)) %>%
  mutate(RE_emp_SE = 100*(emp.SE_SRS - emp.SE_Stratified)/emp.SE_SRS, 
         RE2_emp_SE = 100*(emp.SE_SRS - emp.SE_StratOpt)/emp.SE_SRS,
         RE_mean_SE = 100*(mean.SE_SRS - mean.SE_Stratified)/mean.SE_SRS,
         RE_mean_SE2 = 100*(mean.SE_SRS - mean.SE_StratOpt)/mean.SE_SRS)

plot(RE$n.code, RE$RE_mean_SE,type='l',ylim=c(0,80),xlim=c(100,1000),
     ylab="Percent reduction in SE (Stratified vs. SRS)",xlab="Number of documents sampled/coded")
lines(RE$n.code, RE$RE_emp_SE,type='l',col="blue")
lines(RE$n.code, RE$RE2_emp_SE,type='l',col="red")
lines(RE$n.code, RE$RE_mean_SE2,type='l',col="red",lty=2)


