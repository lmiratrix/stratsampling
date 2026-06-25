


## Function that inputs an LLM prompt and essay text and iteratively evaluates to generate numerical quality scores
### This version uses the new ellmer package, which is easier and more convenient but also slightly more expensive to run


Sys.setenv(
    OPENAI_API_KEY = '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'
)

library(ellmer)
library(tidyverse)
dat=read.csv("raw data/persuade_2.0_human_scores_demo_id_github.csv")

dat_indep = dat %>% filter(task=="Independent")
dat_dep = dat %>% filter(task=="Text dependent")

rubric_indep=paste(readLines("raw data/independent.txt", warn = FALSE), collapse = "\n")
rubric_dep=paste(readLines("raw data/source-dependent.txt", warn = FALSE), collapse = "\n")


type_score <- type_object(
    "You are an expert essay grader for students in grades 6-12. Given an essay text, your job is to evaluate and score the overall quality of the essay.
       Use a scale from 1 to 6 , where higher score indicates greater quality. Present your response as a numerical value. 
    
    In addition, provide a confidence rating for each essay indicating how confident you are in your quality score on a scale from 0-100, where: 
   - 86-100: Very confident, clear-cut case
   - 70-85: Fairly confident, minor ambiguity
   - 50-69: Moderate confidence, could reasonably be adjacent score
   - Below 50: Low confidence, difficult to assess
   
   Provide only the overall score and confidence rating for each essay, without any explanation or additional text.",
    quality_score = type_number(
        "Holistic essay quality score, ranging from 1 to 6"
    ),
    confidence_rating=type_number(
        "Confidence rating, ranging from 0 to 100"
    )
)

## First score the independent essays with the grading rubric as the system prompt
text=dat_indep$full_text
grade=dat_indep$grade_level
prom=dat_indep$assignment

p <- interpolate("Please score the following essay, which was written by a student in grade {{grade}} in response to the prompt \"{{prom}}\": {{text}}")


chat_indep <- chat_openai(model="gpt-5-mini", system_prompt=rubric_indep)
out1 <-parallel_chat_structured(chat_indep, p[1:5000], type_score, max_active=2000,rpm=2000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat_indep, p[5001:10000], type_score, max_active=200,rpm=1000, include_cost = T) 
out3 <-parallel_chat_structured(chat_indep, p[10001:nrow(dat_indep)], type_score, max_active=200,rpm=1000, include_cost = T) 

out_indep = rbind(out1,out2,out3)


sum(out_indep$cost) # how much did running this cost?

dat_indep$LLM_score=out_indep$quality_score
dat_indep$LLM_confidence=out_indep$confidence_rating

dat_indep %>% mutate(resid=holistic_essay_score-LLM_score) %>% 
    rename(Yobs=holistic_essay_score, Yhat=LLM_score, ehat=LLM_confidence) %>%
    summarise(cor_Y_Yhat=cor(Yobs, Yhat), 
              cor_Y_ehat=cor(Yobs, ehat),
              cor_Yhat_ehat=cor(Yhat, ehat), 
              cor_resid_Yhat=cor(resid,Yhat),
              cor_resid_ehat=cor(resid,ehat),cor_Yhat_wc=cor(Yhat,word_count),cor_resid_wc=cor(resid,word_count),cor_ehat_wc=cor(ehat,word_count))

rm(out1,out2,out3)

## Do the same for the source-dependent essays
text=dat_dep$full_text
grade=dat_dep$grade_level
prom=dat_dep$assignment

p <- interpolate("Please score the following essay, which was written by a student in grade {{grade}} in response to the prompt \"{{prom}}\": {{text}}")


chat_dep <- chat_openai(model="gpt-4o-mini", system_prompt=rubric_dep)
out1 <-parallel_chat_structured(chat_dep, p[1:6000], type_score, max_active=200,rpm=1000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat_dep, p[6001:nrow(dat_dep)], type_score, max_active=200,rpm=1000, include_cost = T) 

out_dep = rbind(out1,out2)


sum(out_dep$cost) # how much did running this cost?

dat_dep$LLM_score=out_dep$quality_score
dat_dep$LLM_confidence=out_dep$confidence_rating

dat_dep %>% mutate(resid=holistic_essay_score-LLM_score) %>% 
    rename(Yobs=holistic_essay_score, Yhat=LLM_score, ehat=LLM_confidence) %>%
    summarise(cor_Y_Yhat=cor(Yobs, Yhat), 
              cor_Y_ehat=cor(Yobs, ehat),
              cor_Yhat_ehat=cor(Yhat, ehat), 
              cor_resid_Yhat=cor(resid,Yhat),
              cor_resid_ehat=cor(resid,ehat),cor_Yhat_wc=cor(Yhat,word_count),
              cor_resid_wc=cor(resid,word_count),cor_ehat_wc=cor(ehat,word_count))

# Put everything back together
dat_all = rbind(dat_indep, dat_dep)

dat_all %>% mutate(resid=holistic_essay_score-LLM_score) %>% 
    rename(Yobs=holistic_essay_score, Yhat=LLM_score, ehat=LLM_confidence) %>%
    summarise(cor_Y_Yhat=cor(Yobs, Yhat), 
              cor_Y_ehat=cor(Yobs, ehat),
              cor_Yhat_ehat=cor(Yhat, ehat), 
              cor_resid_Yhat=cor(resid,Yhat),
              cor_resid_ehat=cor(resid,ehat),cor_Yhat_wc=cor(Yhat,word_count),cor_resid_wc=cor(resid,word_count),cor_ehat_wc=cor(ehat,word_count))


dat_persuade=dat_all

save(dat_persuade, file="generated data/persuade-data-with-LLM-predictions-v3.RData")


library(ggplot2)

dat_all = dat_all %>% rename(Yobs=holistic_essay_score, Yhat=LLM_score, ehat=LLM_confidence) %>%
    mutate(resid=Yobs-Yhat, log_wc=log(word_count))

scatter_with_r2 <- function(data, xvar, yvar) {
    x <- data[[xvar]]
    y <- data[[yvar]]
    
    fit <- lm(y ~ x)
    r2  <- summary(fit)$r.squared
    
    ggplot(data, aes(x = x, y = y,col=as.factor(Yobs))) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        annotate("text",
                 x = Inf, y = Inf,
                 label = paste0("R² = ", round(r2, 3)),
                 hjust = 1.1, vjust = 1.5, size = 5) +
        labs(x=xvar, y=yvar,col="Human score")+
        theme_minimal()
}

 scatter_with_r2(dat_all, "Yhat", "Yobs")
 scatter_with_r2(dat_all, "ehat", "Yobs")
 scatter_with_r2(dat_all, "Yhat", "ehat")
 scatter_with_r2(dat_all, "word_count", "ehat")
 scatter_with_r2(dat_all, "log_wc", "ehat")
 
 
 
