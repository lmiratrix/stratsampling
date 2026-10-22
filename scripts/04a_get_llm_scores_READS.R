
library(ellmer)
library(tidyverse)
dat=read.csv("raw data/MORE_study_main.csv")
dat = dat %>% arrange(grade, subject)

prompts = data.frame(grade=c(1,1,2,2), subject=rep(c("science","social"),2),
                     prompt=c("Should people be allowed to cut down trees in the rainforest?",
                              "Both Amelia Earhart and Sally Ride deserve to be celebrated. But, if you had to pick just one of these women explorers to celebrate, which one would you choose - Amelia Earhart or Sally Ride? Why?",
                              "Some scientists think the dinosaurs died out after an asteroid struck Earth. But, not everyone agrees. Do you think that an asteroid killed the dinosaurs? Why or why not?",
                              "Both young Leonardo da Vinci and Henry Ford deserve to be celebrated. But, if you had to pick just one of these young inventors to celebrate, which one would you choose - Leonardo da Vinci or Henry Ford? Why?"
                     ))

reads_rubric = " You are an expert essay grader for students in grades 1-2. Given an essay text, your job is to evaluate and score the overall quality of the essay. Use a scale from 0 to 8, where a higher score indicates greater quality. 

Your evaluation should consider the following criteria:
  1. Claim (0-2 points)
  2. Evidence (0-4 points)
  3. Ending (0-1 points)
  4. Clarity/Voice (0-1 points)

  When evaluating an essay, first code for the 4 dimensions above. Then calculate an overall quality score as a function of the dimension scores (claims, reason, ending, clarity/voice). 
  Note that the term `XXX` in an essay indicates a word that was not legible or could not be transcribed from the handwritten essays.
  Present your response as a numerical value.   Provide only the overall score and confidence rating for each essay, without any explanation or additional text."

dat1 = select(dat, s_id, grade, subject, text)
dat1 = merge(dat, prompts, by=c("grade","subject"))
input = dat1 %>% dplyr::select(s_id, prompt, text) %>% dplyr::rename(essay_prompt=prompt, essay_text=text)

prom=input$essay_prompt
text=input$essay_text

p <- interpolate("Please score the following essay, which was written in response to the prompt \"{{prom}}\": {{text}}")




type_score <- type_object(
  "Given an essay text, your job is to evaluate and score the overall quality of the essay.
   Use a scale from 0 to 8, where higher score indicates greater quality.
   Your evaluation should consider the following criteria:
  1. Claim (0-2 points)
  2. Evidence (0-4 points)
  3. Ending (0-1 points)
  4. Clarity/Voice (0-1 points)

  When evaluating an essay, first code for the 4 dimensions above. Then derive an overall quality score as a function of the dimension scores. Present your response as a numerical value. 
  
  In addition, provide a confidence rating for each essay indicating how confident you are in your quality score on a scale from 0-100, where: 
   - 86-100: Very confident, clear-cut case
   - 70-85: Fairly confident, minor ambiguity
   - 50-69: Moderate confidence, could reasonably be adjacent score
   - Below 50: Low confidence, difficult to assess
   
   Provide only the overall score and confidence rating for each essay, without any explanation or additional text. ",
  quality_score = type_number(
    "Holistic essay quality score, ranging from 0 to 8"
  ),
  confidence_rating=type_number(
      "Confidence rating, ranging from 0 to 100"
  )
)

set.seed(1234)
chat <- chat_openai(model="gpt-4o-mini", system_prompt=reads_rubric)
out <-parallel_chat_structured(chat, p, type_score, max_active=2000,rpm=2000, include_cost = T) # takes about 2 min to run
                                                                              # can go faster by increasing max_active, but this increases cost
dat.out = select(dat, s_id, grade, subject, more, score, text) %>% mutate(LLM_score=out$quality_score, 
                                                        LLM_confidence=out$confidence_rating,
                                                        provider=chat$get_provider()@name,
                                                        model=chat$get_model() )


set.seed(1234)

chat2 = chat_mistral(model="ministral-3b-latest", system_prompt=reads_rubric)
out2 <-parallel_chat_structured(chat2, p, type_score, max_active=2000,rpm=2000,include_cost=T) 


dat.out2 = select(dat, s_id, grade, more, subject, score, text) %>% mutate(LLM_score=out2$quality_score, 
                                                       LLM_confidence=out2$confidence_rating,
                                                       provider=chat2$get_provider()@name,
                                                       model=chat2$get_model() )


set.seed(1234)
chat3 <- chat_google_gemini(model="gemini-3.1-flash-lite", system_prompt=reads_rubric)
out3 <-parallel_chat_structured(chat3, p, type_score, max_active=2000, rpm=2000, include_cost=T) 


dat.out3 = select(dat, s_id, grade, subject, more, score, text) %>% mutate(LLM_score=out3$quality_score, 
                                                        LLM_confidence=out3$confidence_rating,
                                                        provider=chat3$get_provider()@name,
                                                        model=chat3$get_model() )


dat.all = rbind(dat.out,dat.out2,dat.out3)
dat.all$model_short = sapply(1:nrow(dat.all), function(x) strsplit(dat.all$model[x], "-")[[1]][1])

dat.all = dat.all %>% rename(Yhat=LLM_score, Yhat.conf=LLM_confidence)


all.pred = dat.all %>% select(-provider, -model) %>% pivot_wider(names_from=model_short, values_from=c(Yhat, Yhat.conf))



all.pred$n_distinct = apply(all.pred %>% select(starts_with("Yhat_")),1,n_distinct)
all.pred$Yhat.mean=rowMeans(all.pred %>% select(starts_with("Yhat_")))


reads.model.info=dat.all%>%group_by(provider,model,model_short)%>%summarise(n=n()) %>% mutate(date_run=Sys.Date())

save(all.pred, reads.model.info, file="generated data/MOREdata-with-all-LLM-predictions.RData")


