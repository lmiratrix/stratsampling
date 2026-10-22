

library(ellmer)
library(tidyverse)
dat=read.csv("raw data/persuade_2.0_human_scores_demo_id_github.csv")


dat = dat %>% arrange(task) %>% mutate(ID=row_number())

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

## Separate independent prompt and text dependent prompt essays and instructions
text=dat_indep$full_text
grade=dat_indep$grade_level
prom=dat_indep$assignment
gtxt = ifelse(is.na(grade), "", paste0("in grade ", grade, " "))

p_indep = interpolate("Please score the following essay, which was written by a student {{gtxt}}in response to the prompt: \"{{prom}}\"\n\n{{text}}")

rm(text,grade,prom,gtxt)

## Dependent interpolated prompts
text=dat_dep$full_text
grade=dat_dep$grade_level
prom=dat_dep$assignment
gtxt = ifelse(is.na(grade), "", paste0("in grade ", grade, " "))

p_dep = interpolate("Please score the following essay, which was written by a student {{gtxt}}in response to the prompt: \"{{prom}}\"\n\n{{text}}")
rm(text,grade,prom,gtxt)



### Now process essays by task

### --------------------------------------------- ###
### ------ SCORE WITH OPENAI (LLM MODEL 1) ------ ### 
### --------------------------------------------- ###
set.seed(1234)

chat_indep <- chat_openai(model="gpt-4o-mini", system_prompt=rubric_indep)
out1 <-parallel_chat_structured(chat_indep, p_indep[1:5000], type_score, max_active=2000,rpm=1000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat_indep, p_indep[5001:10000], type_score, max_active=2000,rpm=1000, include_cost = T) 
out3 <-parallel_chat_structured(chat_indep, p_indep[10001:nrow(dat_indep)], type_score, max_active=2000,rpm=1000, include_cost = T) 

out_indep = rbind(out1,out2,out3)
rm(out1,out2,out3)

chat_dep <- chat_openai(model="gpt-4o-mini", system_prompt=rubric_dep)
out1 <-parallel_chat_structured(chat_dep, p_dep[1:6000], type_score, max_active=2000,rpm=1000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat_dep, p_dep[6001:nrow(dat_dep)], type_score, max_active=2000,rpm=1000, include_cost = T) 

out_dep = rbind(out1,out2)

out = rbind(out_indep, out_dep)
rm(out1,out2)

dat.out = dat %>% select(ID, essay_id_comp:grade_level) %>% 
    mutate(LLM_score=out$quality_score,
           LLM_confidence=out$confidence_rating,
           provider=chat_indep$get_provider()@name,
           model=chat_indep$get_model() )


### --------------------------------------------- ###
### ------ SCORE WITH MISTRAL (LLM MODEL 2) ----- ### 
### --------------------------------------------- ###
set.seed(1234)

chat2_indep <- chat_mistral(model="ministral-3b-latest", system_prompt=rubric_indep)
out1 <-parallel_chat_structured(chat2_indep, p_indep[1:5000], type_score, max_active=500, rpm=1000, include_cost=T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat2_indep, p_indep[5001:10000], type_score, max_active=500,rpm=1000, include_cost = T) 
out3 <-parallel_chat_structured(chat2_indep, p_indep[10001:nrow(dat_indep)], type_score, max_active=500,rpm=1000, include_cost = T) 

out_indep = rbind(out1,out2,out3)
rm(out1,out2,out3)

chat2_dep <- chat_mistral(model="ministral-3b-latest", system_prompt=rubric_dep)
out1 <-parallel_chat_structured(chat2_dep, p_dep[1:6000], type_score, max_active=500,rpm=1000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat2_dep, p_dep[6001:nrow(dat_dep)], type_score, max_active=500,rpm=1000, include_cost = T) 

out_dep = rbind(out1,out2)

out = rbind(out_indep, out_dep)
rm(out1,out2)

dat.out2 = dat %>% select(ID, essay_id_comp:grade_level) %>% 
    mutate(LLM_score=out$quality_score,
           LLM_confidence=out$confidence_rating,
           provider=chat2_dep$get_provider()@name,
           model=chat2_dep$get_model() )


### --------------------------------------------- ###
### ------ SCORE WITH GEMINI (LLM MODEL 3) ------ ### 
### --------------------------------------------- ###
set.seed(1234)

chat3_indep <- chat_google_gemini(model="gemini-3.1-flash-lite", system_prompt=rubric_indep)
out1 <-parallel_chat_structured(chat3_indep, p_indep[1:5000], type_score, max_active=2000,rpm=1000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat3_indep, p_indep[5001:10000], type_score, max_active=2000,rpm=1000, include_cost = T) 
out3 <-parallel_chat_structured(chat3_indep, p_indep[10001:nrow(dat_indep)], type_score, max_active=2000,rpm=1000, include_cost = T) 

out_indep = rbind(out1,out2,out3)
rm(out1,out2,out3)

chat3_dep <- chat_google_gemini(model="gemini-3.1-flash-lite", system_prompt=rubric_dep)
out1 <-parallel_chat_structured(chat3_dep, p_dep[1:6000], type_score, max_active=1000,rpm=1000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat3_dep, p_dep[6001:nrow(dat_dep)], type_score, max_active=1000,rpm=1000, include_cost = T) 

out_dep = rbind(out1,out2)

out = rbind(out_indep, out_dep)
rm(out1,out2)

dat.out3 = dat %>% select(ID, essay_id_comp:grade_level) %>% 
    mutate(LLM_score=out$quality_score,
           LLM_confidence=out$confidence_rating,
           provider=chat3_dep$get_provider()@name,
           model=chat3_dep$get_model() )




## -- Combine scores/output from all models -- ## 

dat.all = rbind(dat.out,dat.out2,dat.out3)
dat.all$model_short = sapply(1:nrow(dat.all), function(x) strsplit(dat.all$model[x], "-")[[1]][1])

dat.all = dat.all %>% rename(Yhat=LLM_score, Yhat.conf=LLM_confidence)


all.pred = dat.all %>% select(-provider, -model) %>% pivot_wider(names_from=model_short, values_from=c(Yhat, Yhat.conf))


names(all.pred)


all.pred$n_distinct = apply(all.pred %>% select(starts_with("Yhat_")),1,n_distinct)
all.pred$Yhat.mean=rowMeans(all.pred %>% select(starts_with("Yhat_")))


persuade.model.info=dat.all%>%group_by(provider,model,model_short)%>%summarise(n=n()) %>% mutate(date_run=Sys.Date())

save(all.pred, persuade.model.info, file="generated data/persuade-data-with-all-LLM-predictions.RData")

 
