load("~/Dropbox/text and RCT sampling math notes/stratsampling/generated data/MOREdata-with-LLM-predictions-v2.RData")


table(dat$LLM_score)

library(tidyverse)

dat0=dat

dat = dat0 %>% filter(grade==1, subject=="science")

strat = cut(dat$LLM_score,breaks=c(0,1,2,4,8),include.lowest = T) 
#strat=cut(dat$LLM_score,breaks=c(0,2,3,4,8),include.lowest=T)
dat$strat = as.factor(strat)
table(dat$strat)

dat %>% group_by(more,strat) %>% mutate(resid=score-LLM_score) %>%
    summarise(bias=mean(resid), mse=mean(resid^2), sd.resid=sd(resid))


# True treatment effect (averaging over all grades and subjects)
top = dat %>% group_by(more) %>% summarise(mean=mean(score), var=var(score), N=n())
tau = diff(top$mean)
SE.tau.hat = sqrt(sum(top$var/top$N))


dat = dat %>% rename(Yobs=score, Yhat=LLM_score, Z=more,strata=strat)
library( stratsampling )

source("scripts/sim_functions.R")
dat$ID=1:nrow(dat)

library(doParallel)
library(foreach)
library(doRNG)


cores <- min(parallel::detectCores() - 2, 10)
cl <- makeCluster(cores)
registerDoParallel(cl)

n_iter=1000

params=expand.grid(n_code=seq(100,900,by=100))

set.seed(1234)

system.time(
    results <- foreach(i=1:nrow(params), .combine = 'rbind', .packages = c("dplyr","stratsampling")) %dorng% {
        nc = params$n_code[i]
        res <- plyr::rdply(n_iter, one_iter_all(dat,n_code=nc, finite_pop=TRUE))  # Run 1000 iterations
        res$n_code=nc
        return(res)
    }
)

stopCluster(cl)


load("results/reads-application-results.RData")

res.out = results %>% mutate(coding_frac=n_code/nrow(dat), LL=est-2*SEhat, UL=est+2*SEhat, covers=ifelse(tau>=LL & tau<=UL,1,0)) %>%
    group_by(method,coding_frac, n_code)%>%
    summarise(mean.est=mean(est), sd.est=sd(est), mean.SE=mean(SE), mean.SEhat=mean(SEhat), mean.varhat=mean(SEhat^2),
              coverage=mean(covers),VR_to_full=1+mean.varhat/SE.tau.hat^2)

SE.reads1=ggplot(res.out,aes(x=coding_frac,y=sd.est,col=method))+geom_point()+geom_line()+theme_minimal()+
    scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),labels=scales::percent)+
    labs(x="Coding Fraction", y="Empirical SE", col="Estimator")+
    ggplot2::ggtitle("Reads application - Stratifying on LLM predictions")


ggplot(res.out,aes(x=n_code,y=mean.SEhat,col=method))+geom_point()+geom_line()+theme_minimal()

VR.reads1=ggplot(res.out,aes(x=coding_frac,y=VR_to_full,col=method))+
    geom_point()+geom_line()+
    theme_minimal()+geom_hline(yintercept = 1,lty=2)+
    ggplot2::labs(y="Variance ratio to full coding",x="Coding Fraction", col="Estimator")+
    scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),labels=scales::percent)+
    ggplot2::ggtitle("Reads application - Stratifying on LLM predictions")


## What if we stratify on the predictions + residuals themselves (as a "best case" scenario -- not possible in practice)
dat = dat %>% mutate(resid=Yobs-Yhat)
sub0 = dat %>% filter(Z==0) %>% select(Yhat, resid) %>% scale()
sub1 = dat %>% filter(Z==1) %>% select(Yhat, resid) %>% scale()


k0=kmeans(sub0,4)
k1=kmeans(sub1,4)

strat=rep(NA,nrow(dat))
strat[dat$Z==0]=k0$cluster
strat[dat$Z==1]=k1$cluster
dat$strata=as.factor(strat)

dat %>% group_by(Z, strata) %>%
    summarise(bias=mean(resid), mse=mean(resid^2), sd.resid=sd(resid))


cl <- makeCluster(cores)
registerDoParallel(cl)

set.seed(1234)


system.time(
    results2 <- foreach(i=1:nrow(params), .combine = 'rbind', .packages = c("dplyr","stratsampling")) %dorng% {
        nc = params$n_code[i]
        res <- plyr::rdply(n_iter, one_iter_all(dat,n_code=nc, finite_pop=TRUE))  # Run 1000 iterations
        res$n_code=nc
        return(res)
    }
)


stopCluster(cl)

res.out2 = results2 %>% mutate(coding_frac=n_code/nrow(dat), LL=est-2*SEhat, UL=est+2*SEhat, covers=ifelse(tau>=LL & tau<=UL,1,0)) %>%
    group_by(method,coding_frac, n_code)%>%
    summarise(mean.est=mean(est), sd.est=sd(est), mean.SE=mean(SE), mean.SEhat=mean(SEhat), mean.varhat=mean(SEhat^2),
              coverage=mean(covers),VR_to_full=1+mean.varhat/SE.tau.hat^2)

SE.reads2=ggplot(res.out2,aes(x=coding_frac,y=sd.est,col=method))+geom_point()+geom_line()+theme_minimal()+
    scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),labels=scales::percent)+
    labs(x="Coding Fraction", y="Empirical SE", col="Estimator")+
    ggplot2::ggtitle("Reads application - Stratifying on LLM predictions + residuals")


ggplot(res.out2,aes(x=n_code,y=mean.SEhat,col=method))+geom_point()+geom_line()+theme_minimal()

VR.reads2=ggplot(res.out2,aes(x=coding_frac,y=VR_to_full,col=method))+
    geom_point()+geom_line()+
    theme_minimal()+geom_hline(yintercept = 1,lty=2)+
    ggplot2::labs(y="Variance ratio to full coding",x="Coding Fraction", col="Estimator")+
    scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),labels=scales::percent)+
    ggplot2::ggtitle("Reads application - Stratifying on LLM predictions + residuals")
