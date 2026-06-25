

Sys.setenv(
    OPENAI_API_KEY = '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'
)

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

chat_barbera <- chat_openai(model="gpt-4.1", system_prompt=rubric_barbera)

out1 <-parallel_chat_structured(chat_barbera, p[1:2000], type_score, max_active=1000,rpm=2000, include_cost = T) # takes 3-5 min to run
out2 <-parallel_chat_structured(chat_barbera, p[2001:nrow(dat)], type_score, max_active=1000,rpm=2000, include_cost = T)

out_barbera = rbind(out1, out2)




sum(out_barbera$cost) # how much did running this cost?
dat$LLM_relevance=out_barbera$relevance
dat$LLM_score=out_barbera$positivity
dat$LLM_confidence=out_barbera$confidence_rating

dat$model="gpt-5.4-mini"
save(dat_barbera, file="generated data/barbera-data-with-LLM-predictions-5.4-mini-v2.RData")



dat %>% mutate(resid=Yobs-LLM_score) %>%
    rename(Yhat=LLM_score, ehat=LLM_confidence) %>%
    summarise(cor_Y_Yhat=cor(Yobs, Yhat),
              cor_Y_ehat=cor(Yobs, ehat),
              cor_Yhat_ehat=cor(Yhat, ehat),
              cor_resid_Yhat=cor(resid,Yhat),
              cor_resid_ehat=cor(resid,ehat),
              cor_Yhat_wc=cor(Yhat,word_count),
              cor_resid_wc=cor(resid,word_count),
              cor_ehat_wc=cor(ehat,word_count))

dat %>% group_by(LLM_relevance) %>% mutate(resid=Yobs-LLM_score) %>%
    rename(Yhat=LLM_score, ehat=LLM_confidence) %>%
    summarise(n=n(), cor_Y_Yhat=cor(Yobs, Yhat),
              cor_Y_ehat=cor(Yobs, ehat),
              cor_Yhat_ehat=cor(Yhat, ehat),
              cor_resid_Yhat=cor(resid,Yhat),
              cor_resid_ehat=cor(resid,ehat),
              cor_Yhat_wc=cor(Yhat,word_count),
              cor_resid_wc=cor(resid,word_count),
              cor_ehat_wc=cor(ehat,word_count))

dat %>% filter(n_coders>=2) %>% group_by(LLM_relevance) %>% mutate(resid=Yobs-LLM_score) %>%
    rename(Yhat=LLM_score, ehat=LLM_confidence) %>%
    summarise(n=n(), cor_Y_Yhat=cor(Yobs, Yhat),
              cor_Y_ehat=cor(Yobs, ehat),
              cor_Yhat_ehat=cor(Yhat, ehat),
              cor_resid_Yhat=cor(resid,Yhat),
              cor_resid_ehat=cor(resid,ehat),
              cor_Yhat_wc=cor(Yhat,word_count),
              cor_resid_wc=cor(resid,word_count),
              cor_ehat_wc=cor(ehat,word_count))


rm(out1,out2)


if (FALSE){
    
    ## Try multiple draws
    
    # Score each article K = 3 times, then average
    K <- 5
    draws <- vector("list", K)
    for (k in 1:5) {
        set.seed(1000+k)
        chat_k <- chat_openai(model = "gpt-4o-mini",
                              system_prompt = rubric_barbera,
                              params=params(temperature=1.8),service_tier = "priority")   # distinct seed per draw
        draws[[k]] <- parallel_chat_structured(chat_k, p, type_score,
                                               max_active = 2000, rpm = 2000,
                                               include_cost = TRUE)
        print(k)
    }
    
    dat$LLM_rel <- rowMeans(sapply(draws, \(d) ifelse(d$relevance=="Yes",1,0)))
    
    dat$LLM_score <- rowMeans(sapply(draws, \(d) d$positivity))
    dat$LLM_confidence <- rowMeans(sapply(draws, \(d) d$confidence_rating))
    
    
    dat %>% mutate(resid=Yobs-LLM_score) %>%
        rename(Yhat=LLM_score, ehat=LLM_confidence) %>%
        summarise(cor_Y_Yhat=cor(Yobs, Yhat),
                  cor_Y_ehat=cor(Yobs, ehat),
                  cor_Yhat_ehat=cor(Yhat, ehat),
                  cor_resid_Yhat=cor(resid,Yhat),
                  cor_resid_ehat=cor(resid,ehat),
                  cor_Yhat_wc=cor(Yhat,word_count),
                  cor_resid_wc=cor(resid,word_count),
                  cor_ehat_wc=cor(ehat,word_count))
    
    dat %>% group_by(LLM_rel) %>% mutate(resid=Yobs-LLM_score) %>%
        rename(Yhat=LLM_score, ehat=LLM_confidence) %>%
        summarise(n=n(), cor_Y_Yhat=cor(Yobs, Yhat),
                  cor_Y_ehat=cor(Yobs, ehat),
                  cor_Yhat_ehat=cor(Yhat, ehat),
                  cor_resid_Yhat=cor(resid,Yhat),
                  cor_resid_ehat=cor(resid,ehat),
                  cor_Yhat_wc=cor(Yhat,word_count),
                  cor_resid_wc=cor(resid,word_count),
                  cor_ehat_wc=cor(ehat,word_count))
}

# Reliability-weighted correlation
dat %>%
    filter(n_coders >= 3) %>%
    summarise(
        cor_all_3plus     = cor(Yobs, LLM_score),
        cor_weighted_inv  = cov.wt(cbind(Yobs, LLM_score),
                                   wt = 1 / pmax(human_sd^2, 0.1),cor=T)$cor[1, 2]
    )


dat_barbera = dat

save(dat_barbera, file="generated data/barbera-data-with-LLM-predictions-v3.RData")


library(ggplot2)

dat = dat %>% rename(Yhat=LLM_score, ehat=LLM_confidence) %>%
    mutate(resid=Yobs-Yhat, log_wc=log(word_count))

scatter_with_r2 <- function(data, xvar, yvar) {
    x <- data[[xvar]]
    y <- data[[yvar]]

    fit <- lm(y ~ x)
    r2  <- summary(fit)$r.squared

    ggplot(data, aes(x = x, y = y,col=as.factor(round(Yobs)))) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        annotate("text",
                 x = Inf, y = Inf,
                 label = paste0("R² = ", round(r2, 3)),
                 hjust = 1.1, vjust = 1.5, size = 5) +
        labs(x=xvar, y=yvar,col="Human score")+
        theme_minimal()
}

scatter_with_r2(dat, "Yhat", "Yobs")
scatter_with_r2(dat, "ehat", "Yobs")
scatter_with_r2(dat, "Yhat", "ehat")
scatter_with_r2(dat, "word_count", "ehat")
scatter_with_r2(dat, "log_wc", "ehat")

