setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = F)
setwd("../")

library(tidyverse)
library(caret)
library(doParallel)
library(randomForest)

# Load full corpus and all text features
load("generated data/persuade_text_covs.RData")

# Choose "pilot" sample of 5,000 (for initial ML training) to generated predicted scores and predicted uncertainty
set.seed(1234)
samp = sample(1:nrow(dat), size=5000)

dat.pilot = dat[samp,]
dat.main = dat[-samp,]


# Train a rf model to predict human-coded scores in pilot sample
registerDoParallel(20)
rf=randomForest(dat.pilot[,-c(1:2)] , dat.pilot$Yobs , mtry=5,
                ntree = 1000 , importance = TRUE ,
                proximity = TRUE)
varImpPlot(rf , sort = TRUE , 10 )


v = as.data.frame(varImpPlot(rf , sort = TRUE , 10 )) %>%
  rename(MSE="%IncMSE",Purity=IncNodePurity) %>% 
  arrange(desc(MSE))


top.vars = rownames(v)[1:10]


pred = predict(rf, dat.main,predict.all=T)


pred.out = data.frame(id=dat.main$id, Yobs=dat.main$Yobs, 
                      yhat.pilot=pred$aggregate, var.yhat.pilot=apply(pred$individual,1,var))
cor(pred.out[,-c(1)])

plot(pred.out$yhat.pilot, pred.out$var.yhat.pilot) # most uncertainty for medium-range scores

sim.dat = pred.out
save(sim.dat, file="generated data/persuade_main_pilot_predictions.RData")

