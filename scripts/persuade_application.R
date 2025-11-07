load("generated data/persuade-data-with-LLM-predictions-v2.RData")

table(dat_persuade$LLM_score)

library(tidyverse)
dat = dat_persuade %>% rename(score=holistic_essay_score)

#strat = cut(dat$LLM_score,breaks=c(0,2,3,4,5, 6),include.lowest = T) 
dat$strat = as.factor(dat$LLM_score)
table(dat$strat)

dat %>% group_by(strat) %>% mutate(resid=score-LLM_score) %>%
    summarise(n=n(), bias=mean(resid), mse=mean(resid^2), sd.resid=sd(resid))


# True mean
top = dat %>% summarise(mean=mean(score), var=var(score), N=n())
tau = top$mean
SE.tau.hat = sqrt(sum(top$var/top$N))


dat = dat %>% rename(Yobs=score, Yhat=LLM_score,strata=strat)
library( stratsampling )

source("scripts/sim_functions.R")
dat$ID=1:nrow(dat)

library(doParallel)
library(foreach)
library(doRNG)

params=expand.grid(n_code=ceiling(nrow(dat)*seq(0.1,0.9,by=0.1)))

cores <- min(parallel::detectCores() - 2, nrow(params))
cl <- makeCluster(cores)
registerDoParallel(cl)

n_iter=1000


set.seed(1234)

system.time(
    results <- foreach(i=1:nrow(params), .combine = 'rbind', .packages = c("dplyr","stratsampling")) %dorng% {
        nc = params$n_code[i]
        res <- plyr::rdply(n_iter, one_iter_all(dat,n_code=nc, Z=NULL, finite_pop=TRUE))  # Run 1000 iterations
        res$n_code=nc
        return(res)
    }
)

stopCluster(cl)


load("results/persuade-application-results.RData")

res.out = results %>% mutate(coding_frac = n_code/nrow(dat), LL=est-2*SEhat, UL=est+2*SEhat, covers=ifelse(tau>=LL & tau<=UL,1,0)) %>%
    group_by(method,coding_frac, n_code)%>%
    summarise(mean.est=mean(est), sd.est=sd(est), mean.SE=mean(SE), mean.SEhat=mean(SEhat), mean.varhat=mean(SEhat^2),
              coverage=mean(covers),VR_to_full=1+mean.varhat/SE.tau.hat^2)

SE.persuade1=ggplot(res.out,aes(x=coding_frac,y=sd.est,col=method))+geom_point()+geom_line()+theme_minimal()+    
    scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),labels=scales::percent) + labs(y="Empirical SE", x="Coding Fraction",col="Estimator")


ggplot(res.out,aes(x=n_code,y=mean.SEhat,col=method))+geom_point()+geom_line()+theme_minimal()

VR.persuade1=ggplot(res.out,aes(x=coding_frac,y=VR_to_full,col=method))+
    geom_point()+geom_line()+
    theme_minimal()+geom_hline(yintercept = 1,lty=2)+
    ggplot2::labs(x="Coding Fraction", y="Variance ratio to full coding",x="Coding Fraction", col="Estimator")+
    scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),labels=scales::percent) + 
    ggplot2::ggtitle("Persuade application - Stratifying on LLM predictions")

## What if we stratify on the predictions + residuals themselves (as a "best case" scenario -- not possible in practice)
dat = dat %>% mutate(resid=Yobs-Yhat)
sub = dat %>% select(Yhat, resid) %>% scale()


k=kmeans(sub,4)

dat$strata=as.factor(k$cluster)

dat %>% group_by(strata) %>%
    summarise(n=n(), bias=mean(resid), mse=mean(resid^2), sd.resid=sd(resid))


cl <- makeCluster(cores)
registerDoParallel(cl)

set.seed(1234)


system.time(
    results2 <- foreach(i=1:nrow(params), .combine = 'rbind', .packages = c("dplyr","stratsampling")) %dorng% {
        nc = params$n_code[i]
        res <- plyr::rdply(n_iter, one_iter_all(dat,n_code=nc, Z=NULL, finite_pop=TRUE))  # Run 1000 iterations
        res$n_code=nc
        return(res)
    }
)


stopCluster(cl)

res.out2 = results2 %>% mutate(coding_frac=n_code/nrow(dat), LL=est-2*SEhat, UL=est+2*SEhat, covers=ifelse(tau>=LL & tau<=UL,1,0)) %>%
    group_by(method,coding_frac, n_code)%>%
    summarise(mean.est=mean(est), sd.est=sd(est), mean.SE=mean(SE), mean.SEhat=mean(SEhat), mean.varhat=mean(SEhat^2),
              coverage=mean(covers),VR_to_full=1+mean.varhat/SE.tau.hat^2)

SE.persuade2=ggplot(res.out2,aes(x=coding_frac,y=sd.est,col=method))+geom_point()+geom_line()+theme_minimal()+    
    scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),labels=scales::percent) + labs(y="Empirical SE", x="Coding Fraction",col="Estimator")

ggplot(res.out2,aes(x=n_code,y=mean.SEhat,col=method))+geom_point()+geom_line()+theme_minimal()

VR.persuade2=ggplot(res.out2,aes(x=coding_frac,y=VR_to_full,col=method))+
    geom_point()+geom_line()+
    theme_minimal()+geom_hline(yintercept = 1,lty=2)+
    ggplot2::labs(x="Coding Fraction", y="Variance ratio to full coding",x="Coding Fraction", col="Estimator")+
    scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),labels=scales::percent) + 
    ggplot2::ggtitle("Persuade application - Stratifying on LLM predictions + residuals")
