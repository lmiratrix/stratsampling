
library(ellmer)
library(tidyverse)

# Load the article-level Barbera finite population built by
# scripts/00_prep_barbera_data.R. Each row is one NYT article with
# cleaned text, metadata, and the human-coded mean positivity (Yobs).
load("generated data/barbera-article-data.RData")
dat = dat_barbera

# -----------------------------------------------------------------------------
# Build coder-task prompt ----
# -----------------------------------------------------------------------------


rubric_barbera <- "
In this job, you will be presented with excerpts from newspaper articles. We will ask you to read each excerpt and then answer a few questions about whether the excerpt provides information about how the US economy is doing.

Process
1. Read the article
2. Determine if the article gives an indication of how the United States economy is performing.

For example, an article with this headline probably gives you an indication of how the U.S. economy is performing: \"Unemployment just went up by 2%.\"

Note that you are only being asked whether the article tells you something about how the United States economy is doing -- so articles that only contain information about the economy in other countries are not relevant.

3. If the article gives you an indication of how the economy is performing, tell us whether the indication is negative or positive, using a scale ranging from 1 (negative) to 9 (positive). Use the full range of the scale. Pick any value from 1 to 9 that best captures
where the article falls between the endpoints.

If the text contains typos or is otherwise corrupted in some way, do your best to make out the message of the article.

Positive tone reflects coverage emphasizing growth, strength, recovery, improvement, optimism, or gains. Negative tone reflects coverage emphasizing decline, weakness, recession, deterioration, pessimism, or losses.
Judge the tone of the article as a typical human reader would perceive it, not the tone intended by the author.

Always provide your best-guess tone score, even when relevance is 'No' or 'Not Sure', so the tone field is never missing.



You are not trying to be correct in an abstract sense. You are trying to match
the judgment a typical CrowdFlower worker would make after reading the segment
once. Use the two worked examples below as your anchor for how workers used
the scale.

When responding, provide only the the relevance label, single-digit positivity score, and confidence rating, without any explanation or additional text.

--- Worked example 1 (target: Q1 = Yes, Q2 = 5) -----------------------------

Headline: Administration Using Study To Push Elderly Drug Plan: Pressing
Republicans to Reach Compromise

Date: 2000-04-26

Segment: WASHINGTON, April 25 - Prices for the 50 drugs most frequently used
by older Americans increased last year at nearly twice the rate of inflation,
a study to be issued on Wednesday at the White House says. The study found a
similar trend over the last six years. President Clinton planned to cite the
data in arguing that drug costs constitute a growing burden for the elderly
and that Congress should pass legislation to reduce the burden. The Democratic
leaders of Congress, Senator Tom Daschle of South Dakota and Representative
Richard A Gephardt of Missouri, intend to join Mr Clinton at the White House.

Reasoning: The article touches on a cost-of-living concern (drug prices rising
faster than inflation), which makes it relevant to the US economy. The tone is
mixed: it reports a negative fact (rising prices) but also a policy response
aimed at relieving the burden. A typical worker reads this as relevant but
not directionally positive or negative; the modal coding is 5.

--- Worked example 2 (target: Q1 = Yes, Q2 = 2) -----------------------------

Headline: Nation's Governors See a Dismal Economic Outlook and a Slow Recovery

Date: 2009-02-22

Segment: WASHINGTON - The nation's governors said Saturday that passage of a
$787 billion bill to stimulate the economy might help them avert draconian
budget cuts, but that they did not expect to see signs of an economic recovery
until late this year or early 2010. The officials, arriving here for the
winter meeting of the National Governors Association, said that state revenues
were coming in far below their projections and that the new federal measure,
while helpful, would not be a panacea. Gov. Jon Huntsman Jr. of Utah, where
the economy is better than in most states, said the revenue figures were
'still dismal'.

Reasoning: The segment is explicitly about the state of the US national
economy (relevant). The language is strongly negative: 'dismal outlook',
'draconian budget cuts', 'far below projections', 'would not be a panacea'.
A typical worker reads this as clearly negative but not at the extreme of
the scale; the modal coding is 2."


type_score <- type_object(
    "Given an article, make two judgments:

     (1) Relevance: is the article about the US national economy? Answer
         'Yes', 'No', or 'Not Sure'.
     (2) Tone: Evaluate the overall tone of the article on a continuous (single-digit) 9-point scale 
        where 1 is very negative and 9 is very positive.

     Also provide a confidence rating for your positivity score on a scale
     from 0-100, where:
       - 86-100: Very confident, clear-cut case
       - 70-85: Fairly confident, minor ambiguity
       - 50-69: Moderate confidence, could reasonably be adjacent score
       - Below 50: Low confidence, difficult to assess

     Provide only the the relevance label, single-digit positivity score, and confidence rating, without any explanation or additional text.",
    relevance = type_enum(
        "Is the article about the US national economy?",
        values = c("Yes", "No", "Not Sure")
    ),
    positivity = type_number(
        "Article-level tone toward the US national economy, on a continuous (single-digit) 9-point scale from 1 (negative) to 9 (positive)."
    ),
    confidence_rating=type_number(
        "Confidence rating, ranging from 0 to 100"
    )
)



# -----------------------------------------------------------------------------
# Prepare per-article user prompts ----
# -----------------------------------------------------------------------------

SENT_ABBREV <- c("Mr", "Mrs", "Ms", "Mx", "Dr", "Jr", "Sr", "St",
                 "Prof", "Rev", "Hon", "Gov", "Sen", "Rep",
                 "Capt", "Lt", "Sgt", "Col", "Gen", "Adm",
                 "Inc", "Corp", "Ltd", "Co", "No", "vs", LETTERS)

first_segment <- function(x, n_sents = 5, max_words = 175) {
    if (is.na(x) || !nzchar(x)) return(x)
    
    # Mask periods inside known abbreviations so they are not flagged as
    # sentence terminators. The replacement character is a section sign, which
    # preserves character positions so offsets map back to the original text.
    pattern <- paste0("\\b(", paste(SENT_ABBREV, collapse = "|"), ")\\.")
    masked  <- str_replace_all(x, pattern, "\\1\u00A7")
    
    # Also mask periods inside single-letter initial sequences such as U.S.,
    # U.K., U.S.A., and J.R.R. so "U.S. economy" is not split mid-token.
    masked <- str_replace_all(masked, "(?<=\\b[A-Z])\\.(?=[A-Z]\\b)", "\u00A7")
    
    # Locate end-of-sentence boundaries (. ! ?) followed by whitespace or EOS.
    ends <- str_locate_all(masked, "[.!?]+(?=\\s|$)")[[1]]
    
    if (nrow(ends) >= n_sents) {
        cut <- ends[n_sents, "end"]
        seg <- str_sub(x, 1, cut)
    } else {
        seg <- x
    }
    
    # Word-count safeguard for unusually long leads.
    words <- str_split(seg, "\\s+")[[1]]
    if (length(words) > max_words) {
        seg <- paste(words[seq_len(max_words)], collapse = " ")
    }
    
    str_trim(seg)
}


dat <- dat %>%
    mutate(text_segment = vapply(text, first_segment, character(1)))

head    = dat$headline
dt      = dat$date
segment = dat$text_segment


p <- interpolate("Please read and code the following article segment.
First, judge whether the article is about the US national economy (Yes / No / Not Sure).
Then rate its overall tone toward the US national economy on the 9-point scale described in the instructions.
Base both judgments only on the article segment provided below, which
reproduces the opening of the article exactly as the original CrowdFlower
workers saw it. Match the style of the two worked examples in the system
prompt.

Headline: {{head}}
Date: {{dt}}

Article segment:
{{segment}}")


# -----------------------------------------------------------------------------
# Score all articles in batches ----
# -----------------------------------------------------------------------------

chat <- chat_openai(model="gpt-4o-mini", system_prompt=rubric_barbera)

out1 <-parallel_chat_structured(chat, p[1:2000], type_score, max_active=1000,rpm=2000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat, p[2001:nrow(dat)], type_score, max_active=1000,rpm=2000, include_cost = T)

out = rbind(out1, out2)

dat.out = select(dat, ID, article_id, date, year, decade, headline, 
                 text, word_count, n_coders, Yobs, text_segment) %>% 
    mutate(LLM_relevance=out$relevance,
           LLM_positivity=out$positivity, 
           LLM_confidence=out$confidence_rating,
           provider=chat$get_provider()@name,
           model=chat$get_model() )

rm(out1,out2)


set.seed(1234)

chat2 <- chat_mistral(model="ministral-3b-latest",system_prompt=rubric_barbera)
out1 <-parallel_chat_structured(chat2, p[1:2000], type_score, max_active=1000,rpm=500) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat2, p[2001:nrow(dat)], type_score, max_active=1000,rpm=500)
out2 = rbind(out1, out2)


dat.out2 = select(dat, ID, article_id, date, year, decade, headline, 
                  text, word_count, n_coders, Yobs, text_segment) %>% 
    mutate(LLM_relevance=out2$relevance,
           LLM_positivity=out2$positivity, 
           LLM_confidence=out2$confidence_rating,
           provider=chat2$get_provider()@name,
           model=chat2$get_model() )

set.seed(1234)

chat3 <- chat_google_gemini(model="gemini-3.1-flash-lite", system_prompt=rubric_barbera)
out3 <-parallel_chat_structured(chat3, p, type_score, max_active=1000, rpm=2000) 


dat.out3 = select(dat, ID, article_id, date, year, decade, headline, 
                  text, word_count, n_coders, Yobs, text_segment) %>% 
    mutate(LLM_relevance=out3$relevance,
           LLM_positivity=out3$positivity, 
           LLM_confidence=out3$confidence_rating,
           provider=chat3$get_provider()@name,
           model=chat3$get_model() )


## -- Combine scores/output from all models -- ## 

dat.all = rbind(dat.out,dat.out2,dat.out3)
dat.all$model_short = sapply(1:nrow(dat.all), function(x) strsplit(dat.all$model[x], "-")[[1]][1])

dat.all = dat.all %>% rename(rel=LLM_relevance, Yhat=LLM_positivity, Yhat.conf=LLM_confidence)


all.pred = dat.all %>% select(-provider, -model) %>% pivot_wider(names_from=model_short, values_from=c(rel, Yhat, Yhat.conf))


names(all.pred)


all.pred$n_distinct = apply(all.pred %>% select(starts_with("Yhat_")),1,n_distinct)
all.pred$Yhat.mean=rowMeans(all.pred %>% select(starts_with("Yhat_")))


barbera.model.info=dat.all%>%group_by(provider,model,model_short)%>%summarise(n=n()) %>% mutate(date_run=Sys.Date())

save(all.pred, barbera.model.info, file="generated data/barbera-data-with-all-LLM-predictions.RData")

