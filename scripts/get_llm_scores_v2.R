


## Function that inputs an LLM prompt and essay text and iteratively evaluates to generate numerical quality scores
### This version uses the new ellmer package, which is easier and more convenient but also slightly more expensive to run


Sys.setenv(
  OPENAI_API_KEY = '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'
)

library(ellmer)
library(tidyverse)
dat=read.csv("~/Dropbox/Text grant - Methods paper/power-rep/raw data/MORE_study_main.csv")
dat = dat %>% arrange(grade, subject)

prompts = data.frame(grade=c(1,1,2,2), subject=rep(c("science","social"),2),
                     prompt=c("Should people be allowed to cut down trees in the rainforest?",
                              "Both Amelia Earhart and Sally Ride deserve to be celebrated. But, if you had to pick just one of these women explorers to celebrate, which one would you choose - Amelia Earhart or Sally Ride? Why?",
                              "Some scientists think the dinosaurs died out after an asteroid struck Earth. But, not everyone agrees. Do you think that an asteroid killed the dinosaurs? Why or why not?",
                              "Both young Leonardo da Vinci and Henry Ford deserve to be celebrated. But, if you had to pick just one of these young inventors to celebrate, which one would you choose - Leonardo da Vinci or Henry Ford? Why?"
                     ))


dat1 = select(dat, s_id, grade, subject, text)
dat1 = merge(dat, prompts, by=c("grade","subject"))
input = dat1 %>% dplyr::select(s_id, prompt, text) %>% dplyr::rename(essay_prompt=prompt, essay_text=text)

prom=input$essay_prompt
text=input$essay_text

p <- interpolate("Please score the following essay, which was written in response to the prompt \"{{prom}}\": {{text}}")


type_score <- type_object(
  "You are an expert essay grader for students in grades 1-2.
        Given an essay text, your job is to evaluate and score the overall quality of the essay.
   Use a scale from 0 to 8, where higher score indicates greater quality.
   Your evaluation should consider the following criteria:
  1. Claim (0-2 points)
  2. Evidence (0-4 points)
  3. Ending (0-1 points)
  4. Clarity/Voice (0-1 points)

  When evaluating an essay, first code for the 4 dimensions above. Then derive an overall quality score as a function of the dimension scores (claims, reason, ending,
  clarity/voice). Note that the term `XXX` in an essay indicates a word that was not legible or could not be transcribed from the handwritten essays.
  Present your response as a numerical value. Provide only the overall score for each essay, without any explanation or additional text.

  ",
  quality_score = type_number(
    "Holistic essay quality score, ranging from 0 to 8"
  )
)

chat <- chat_openai(model="gpt-4o-mini")
out <-parallel_chat_structured(chat, p, type_score, max_active=200,rpm=2000, include_cost = T) # takes about 2 min to run
                                                                              # can go faster by increasing max_active, but this increases cost

sum(out$cost) # how much did running this cost?

all.equal(dat$s_id, input$s_id)
dat$LLM_score=out$quality_score

cor(dat$score, out$quality_score)


save(dat, file="generated data/MOREdata-with-LLM-predictions-v2.RData")



library(reticulate)

# Install Python tiktoken if not present
reticulate::py_install("tiktoken", pip = TRUE)

tiktoken <- import("tiktoken")

# Helper function: count tokens and compute cost
estimate_tokens_cost <- function(text, model = "gpt-4o-mini", price_per_1M = 0.15) {
    enc <- tiktoken$encoding_for_model(model)
    all_text = paste(text, collapse = " ")

    tokens <- enc$encode(all_text)
    n_tokens <- length(tokens)
    cost <- (n_tokens / 1000000) * price_per_1M
    list(tokens = n_tokens, cost = sprintf("$%.6f", cost))
}

# Example usage in your script:
prompt <- "This is a test message to send to the API."
estimate <- estimate_tokens_cost(prompt, model = "gpt-4o-mini")

cat("Prompt tokens:", estimate$tokens, "\n")
cat("Estimated cost (USD):", estimate$cost, "\n")



