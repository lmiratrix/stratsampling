

library(tidyverse)
library(lubridate)

# -----------------------------------------------------------------------------
# Read raw coder-level file ----
# Filter to relevant, tone-graded coder rows 
# -----------------------------------------------------------------------------

# One row per coder-article. Columns: articleid, headline, date, text,
# relevance ("yes"/"no"/missing), positivity (1-9 scale).
raw <- read.csv("Barbera 2020/data/5AC.csv",
                stringsAsFactors = FALSE,
                na.strings = c("", "NA"))

coder <- raw %>%
  filter(relevance == "yes") %>%
  filter(!is.na(positivity)) %>%
  mutate(positivity = as.numeric(positivity))


# -----------------------------------------------------------------------------
# Clean article text ----
# -----------------------------------------------------------------------------

# The raw `text` column uses </br></br> as a paragraph separator and carries
# a variety of smart quotes/dashes from the NYT scrape. Strip the HTML breaks,
# normalize whitespace, and drop curly punctuation so the text is ready for
# word-count features, embeddings, and LLM prompting.
clean_text <- function(x) {
    x %>%
        str_replace_all("</br>", " ") %>%
        str_replace_all("<br\\s*/?>", " ") %>%
        str_replace_all("[\u2018\u2019\u201A\u201B]", "'") %>%
        str_replace_all("[\u201C\u201D\u201E\u201F]", '"') %>%
        str_replace_all("[\u2013\u2014]", "-") %>%
        str_replace_all("[[:space:]]+", " ") %>%
        str_trim()
}

coder <- coder %>%
  mutate(text_clean = clean_text(text))


# -----------------------------------------------------------------------------
# Collapse to article level ----
# -----------------------------------------------------------------------------

article <- coder %>%
  group_by(articleid) %>%
  summarise(
    headline    = first(headline),
    date        = first(date),
    text        = first(text_clean),
    n_coders    = n(),
    human_mean  = mean(positivity),
    human_sd    = ifelse(n() > 1, sd(positivity), NA_real_),
    human_min   = min(positivity),
    human_max   = max(positivity),
    .groups = "drop"
  )


article <- article %>%
  mutate(
    date    = ymd(date),
    year    = year(date),
    month   = month(date),
    yrmon   = format(date, "%Y-%m"),
    quarter = paste0(year, "-Q", quarter(date)),
    decade  = paste0(floor(year / 10) * 10, "s")
  )


article <- article %>%
  mutate(
    word_count     = str_count(text, "\\S+"),
    char_count     = nchar(text),
    sentence_count = str_count(text, "[.!?]+")
  )


# -----------------------------------------------------------------------------
# Final assembly ----
# -----------------------------------------------------------------------------

# Standardize to the Yobs naming used elsewhere in the package so the
# downstream application scripts (stratification, estimation, simulation)
# can be applied without additional renaming.
dat_barbera <- article %>%
  rename(article_id = articleid) %>%
  mutate(
    ID   = row_number(),
    Yobs = human_mean
  ) %>%
  select(ID, article_id, date, year, month, yrmon, quarter, decade, recession,
         headline, text, word_count, char_count, sentence_count,
         n_coders, human_mean, human_sd, human_min, human_max, Yobs)


print(table(dat_barbera$decade))

print(summary(dat_barbera$Yobs))


# -----------------------------------------------------------------------------
# Save outputs ----
# -----------------------------------------------------------------------------

# Article-level finite population: one row per article, human outcome +
# metadata + cleaned text. This is the object downstream scripts will load.
save(dat_barbera, file = "generated data/barbera-article-data.RData")

# Coder-level file for reliability analyses
# One row per coder-article after filtering.
dat_barbera_coder <- coder %>%
  rename(article_id = articleid) %>%
  select(article_id, positivity, relevance) %>%
  mutate(coder_row = row_number())

save(dat_barbera_coder, file = "generated data/barbera-coder-data.RData")

