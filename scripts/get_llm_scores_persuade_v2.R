


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
    "Given an essay text, your job is to evaluate and score the overall quality of the essay.
  Present your response as a numerical value. Provide only the overall score for each essay, without any explanation or additional text.",
    quality_score = type_number(
        "Holistic essay quality score, ranging from 1 to 6"
    )
)

## First score the independent essays with the grading rubric as the system prompt
text=dat_indep$full_text
grade=dat_indep$grade_level
prom=dat_indep$assignment

p <- interpolate("Please score the following essay, which was written by a student in grade {{grade}} in response to the prompt \"{{prom}}\": {{text}}")


chat_indep <- chat_openai(model="gpt-4o-mini", system_prompt=rubric_indep)
out1 <-parallel_chat_structured(chat_indep, p[1:5000], type_score, max_active=200,rpm=1000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat_indep, p[5001:10000], type_score, max_active=200,rpm=1000, include_cost = T) 
out3 <-parallel_chat_structured(chat_indep, p[10001:nrow(dat_indep)], type_score, max_active=200,rpm=1000, include_cost = T) 

out_indep = rbind(out1,out2,out3)


sum(out_indep$cost) # how much did running this cost?

dat_indep$LLM_score=out_indep$quality_score

cor(dat_indep$holistic_essay_score, out_indep$quality_score)
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

cor(dat_dep$holistic_essay_score, out_dep$quality_score)


dat_all = rbind(dat_indep, dat_dep)
cor(dat_all$holistic_essay_score, dat_all$LLM_score)



dat_persuade=dat_all

save(dat_persuade, file="generated data/persuade-data-with-LLM-predictions-v2.RData")

