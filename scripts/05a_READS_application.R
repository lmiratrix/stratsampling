
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(doRNG)
library(kableExtra)
library(patchwork)
library(rpart)
library(gt)

#library(plyr)


# Load the data
load("generated data/MOREdata-with-LLM-predictions-v3.RData")

stats = dat %>% group_by( more ) %>% 
  summarise( mn = mean( score ),
             sd_Y = sd( score ), 
             n = n(), 
             sd_Yhat=sd(LLM_score))

avg_stats = stats %>% 
  summarise( sd_bar_Y = sqrt( weighted.mean( sd_Y^2, w=n ) ),
             sd_bar_Yhat=sqrt(weighted.mean(sd_Yhat^2,w=n)))


dat$Yobs_raw=dat$score
dat$Yhat_raw=dat$LLM_score
dat$score=dat$score/avg_stats$sd_bar_Y
dat$LLM_score=dat$LLM_score/avg_stats$sd_bar_Yhat

dat <- dat %>% 
  select(text, score, LLM_score, LLM_confidence, more, 
         s_id, t_id, sch_id, Yobs_raw,Yhat_raw) %>%
  rename(Yobs = score, 
         Yhat = LLM_score, 
         ehat=LLM_confidence,
         Z = more) %>%
  mutate(ID = row_number())


# Basic descriptives
cat("Sample size:", nrow(dat), "\n")
cat("Treatment:", sum(dat$Z), "Control:", sum(1-dat$Z), "\n")
cat("Mean outcome - Treatment:", mean(dat$Yobs[dat$Z==1]), "\n")
cat("Mean outcome - Control:", mean(dat$Yobs[dat$Z==0]), "\n")

# Overall Correlation and R²
overall_r = cor(dat$Yobs, dat$Yhat)
overall_r2 = cor(dat$Yobs, dat$Yhat)^2
SE_r = sqrt((1-overall_r^2)/(nrow(dat)-2))

cat(sprintf("Correlation: %.3f (SE=%.3f)\nOverall R² between human and LLM scores: %.3f\n", overall_r, SE_r, overall_r2))



# By treatment group
r2_by_treatment <- dat %>%
  group_by(Z) %>%
  summarise(
    R2 = cor(Yobs, Yhat)^2,
    cor = cor(Yobs, Yhat),
    N = n(),
    SE_cor = round(sqrt((1-cor^2)/(N-2)),2)
    
  )
print(r2_by_treatment)

# Calculate residuals
dat = dat %>% mutate(resid=Yobs-Yhat,
                     abs_resid=abs(resid),
                     sq_resid=resid^2)

# Summary statistics
residual_stats <- dat %>%
  group_by(Z) %>%
  summarise(
    mean_resid = mean(resid),
    sd_resid = sd(resid),
    mad_resid = mad(resid),
    min_resid = min(resid),
    max_resid = max(resid)
  )
print(residual_stats)

# Scatter plot of predictions vs. observed
p1 <- ggplot(dat, aes(x = Yhat, y = Yobs, color = factor(Z))) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~Z, labeller = labeller(Z = c("0" = "Control", "1" = "Treatment"))) +
  labs(title = "LLM Predictions vs. Human Scores",
       x = "LLM Prediction", y = "Human Score",
       color = "Treatment") +
  theme_bw()

# Residual distribution
p2 <- ggplot(dat, aes(x = resid, fill = factor(Z))) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Z, labeller = labeller(Z = c("0" = "Control", "1" = "Treatment"))) +
  labs(title = "Distribution of Residuals",
       x = "Residual (Human - LLM)", y = "Count") +
  theme_bw()

# Residuals vs predictions (check for heteroskedasticity)
p3 <- ggplot(dat, aes(x = Yhat, y = abs_resid, color = factor(Z))) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~Z) +
  labs(title = "Absolute Residuals vs. Predictions",
       x = "LLM Prediction", y = "|Residual|") +
  theme_bw()

gridExtra::grid.arrange(p1, p2, p3, nrow = 3)

#ggsave("figures/empirical_prediction_quality.pdf",  width = 10, height = 12)

# -----------------------------------------------------------------------------
#  Explore Potential Stratification Variables ----
# -----------------------------------------------------------------------------

## Step 2.1: Extract Text Features
# Calculate text features
dat <- dat %>%
  mutate(
    word_count = str_count(text, "\\S+"),
    char_count = nchar(text),
    sentence_count = str_count(text, "[.!?]+"),
    avg_word_length = char_count / word_count,
    avg_sentence_length = word_count / sentence_count
  )

# Create categorical versions
dat <- dat %>%
  mutate(
    length_cat = cut(word_count, 
                     breaks = c(0, 20, 40, 80, Inf),
                     labels = c("Very Short", "Short", "Medium", "Long")),
    pred_cat = cut(Yhat, 
                   breaks = c(-Inf, 2, 3, 4, Inf),
                   labels = c("Low", "Medium", "High", "Very High"))
  )

## Examine Residual Patterns by Features
# By word count categories
resid_by_length <- dat %>%
  group_by(Z, length_cat) %>%
  summarise(
    N = n(),
    mean_resid = mean(resid),
    var_resid = var(resid),
    sd_resid = sd(resid),
    .groups = "drop"
  )
print(resid_by_length)

# By prediction level
resid_by_pred <- dat %>%
  group_by(Z, pred_cat) %>%
  summarise(
    N = n(),
    mean_resid = mean(resid),
    var_resid = var(resid),
    sd_resid = sd(resid),
    .groups = "drop"
  )
print(resid_by_pred)

# Visualize
p_length <- ggplot(resid_by_length, aes(x = length_cat, y = var_resid, fill = factor(Z))) +
  geom_col(position = "dodge") +
  labs(title = "Residual Variance by Essay Length",
       x = "Essay Length", y = "Var(residual)", fill = "Treatment") +
  theme_bw()

p_pred <- ggplot(resid_by_pred, aes(x = pred_cat, y = var_resid, fill = factor(Z))) +
  geom_col(position = "dodge") +
  labs(title = "Residual Variance by Prediction Level",
       x = "LLM Prediction", y = "Var(residual)", fill = "Treatment") +
  theme_bw()

gridExtra::grid.arrange(p_length, p_pred, nrow = 2)

#ggsave("figures/empirical_residual_heterogeneity.pdf",
#  gridExtra::grid.arrange(p_length, p_pred, nrow = 2)#,
# width = 10, height = 8)

# -----------------------------------------------------------------------------
#  Create Candidate Stratifications ----
# -----------------------------------------------------------------------------
library(stratsampling)


# Create candidate stratifications
candidates <- create_stratification_candidates(dat, vars = c("Yhat",  "ehat", "word_count"),Z="Z")

# Evaluate all candidates
comparison_oracle <- compare_stratifications(dat, candidates, min_var_threshold = 0.01) # comparison using residuals (not feasible in practice)
comparison_precoding <- compare_stratifications_precoding(dat, candidates,Z="Z") # comparison based on Yhat only (can be done before coding has started)


# Select top candidates meeting criteria
top_candidates <- comparison_oracle$summary %>%
  filter(avg_size_ratio < 10, min_stratum_size > 100, min_var_resid > 0.01) %>%
  arrange(desc(mean_spread)) %>%
  head(5)

cat("\n=== Top 5 Candidate Stratifications ===\n")
cat("Based on oracle comparisons (using known residuals Yobs-Yhat):")
print(top_candidates)

cat("Based on pre-coding comparisons (using only predictions Yhat): ")
head(comparison_precoding,n=5)

# -----------------------------------------------------------------------------
#  Compare Candidates using Oracle Information ----
# -----------------------------------------------------------------------------

## Plot pre-coding proxy for stratification benefit vs. oracle
oracle_vs_precoding <- compare_oracle_vs_precoding(
    oracle_comparison = comparison_oracle,
    precoding_comparison = comparison_precoding,
    top_n = 18
)


# View scatter plot
print(oracle_vs_precoding$scatter_plot)
ggsave("figures/READS_precoding_strat_validation.pdf",oracle_vs_precoding$scatter_plot,  width = 8, height = 6)



# Compare variance decomposition oracle vs. precoding
results <- compare_decomposition_oracle(dat, candidates, n_code = floor(0.3*nrow(dat)),Z="Z")
results$summary  # View table
results$plot     # Display visualization
reads_decomp=plot_decomposition_comparison(results,show_components = F)
reads_decomp
ggsave("figures/READS_precoding_strat_variance-reduction.pdf",reads_decomp,  width = 8, height = 6)


## Select Final Stratification
# Apply best stratification
strat_name="Yhat_quintiles"
dat$strata <- candidates[[strat_name]] # choose one that doesn't use residuals (oracle)
dat$strata = factor(dat$strata, labels=c("0-1","2","3","4","5-7"))
plot_strata_residuals(dat, dat$strata)

compute_variance_decomposition(dat, strata=dat$strata)

strata_table <- dat %>%
  group_by(Z, strata) %>%
  summarise(
    N_zk = n(),
    mean_Yobs = mean(Yobs),
    sd_Yobs = sd(Yobs),
    mean_Yhat = mean(Yhat),
    mean_conf=mean(ehat),
    mean_resid = mean(resid),
    var_resid = var(resid),
    sd_resid = sd(resid),
    mean_word_count = mean(word_count), .groups = "drop") %>%
  mutate(
    k=rep(1:n_distinct(dat$strata),2), 
    Z = Z,
    `Predicted Quality` = gsub(" Yhat", "", strata),
    `Mean Y` = sprintf("%.2f", mean_Yobs),
    `Mean Yhat` = sprintf("%.2f", mean_Yhat),
    `Mean Resid.` = sprintf("%.2f", mean_resid),
    `Var(Resid.)` = sprintf("%.2f", var_resid),
  ) %>%
  select( Z,k, `Predicted Quality`,  N_zk,
          `Mean Resid.`, `Var(Resid.)`) %>%
  rename( `$Z$`      = Z,
          `$k$`      = k,
          `$N_{zk}$` = N_zk,
          `$\\hat{Y}$` = `Predicted Quality`
  )

print(strata_table)

kbl(strata_table, format = "latex", booktabs = TRUE,  
    caption = "Characteristics of final stratification for the MORE study determined by LLM-generated quality score (quintiles).",
    label="reads_strat_v2",linesep=c("", "","", "", "\\addlinespace"),escape=F) %>%   kable_styling(latex_options = "hold_position") %>% 
  write_file("tables/READS_empirical_strata_characteristics_v2.tex")


# -----------------------------------------------------------------------------
# Estimates output from a single run across all methods (Example of practical real-world use) ----
# -----------------------------------------------------------------------------
source("scripts/application_functions.R")  # From previous code
source("scripts/sim_functions.R")


# Single sample at 30% coding budget
set.seed(2025)
single_estimates <- compute_single_sample_estimates(dat, coding_fraction = 0.30, finite_pop=FALSE)

comparison_table <- create_comparison_table(single_estimates)
print(comparison_table)

comparison_table=comparison_table %>% 
  mutate( CI=sprintf("(%.2f, %.2f)", `95% CI Lower`, `95% CI Upper`), 
          `Var. Inflation`=sprintf("%.2f", `Var. Inflation`)) %>%
  select(Method, `n Coded`, Estimate, SE, CI, `Var. Inflation`)
comparison_table

kbl(comparison_table,  format = "latex", booktabs = TRUE,  
    caption = "Results for the MORE study estimating the impact of the MORE intervention on writing quality scores from a single iteration with coding budget of h=30\\%. Variance inflation is relative to full human coding (oracle).",
    label="reads_point_ests_v2",escape=F) %>%   kable_styling(latex_options = "hold_position") %>% 
  write_file("tables/READS-ests-single-iter.tex")


# Create forest plot
oracle <- est_oracle(dat)

forest_plot <- create_forest_plot(single_estimates, oracle_est = oracle$est)
print(forest_plot)
ggsave("figures/READS_emp_single_forest.pdf",forest_plot,  width = 8, height = 6)


# Multiple samples to show variability

multi_estimates <- compute_multiple_sample_estimates(dat, coding_fraction = 0.30, 
                                                     n_samples = 20, seed = 2026, finite_pop=F)
multi_summary <- create_multi_sample_summary(multi_estimates)
print(multi_summary)

# Multi-sample forest plot
forest_multi <- create_multi_sample_forest_plot(multi_estimates, oracle_est = oracle$est,jitter_seed=1234)
print(forest_multi)
ggsave("figures/READS_emp_multi_forest.pdf",forest_multi,  width = 8, height = 6)


scatter=oracle_vs_precoding$scatter_plot

reads_figs=list(decomp=reads_decomp, forest_multi=forest_multi, forest_plot=forest_plot, scatter=scatter)
names(reads_figs)=paste0("reads_",names(reads_figs))
save(reads_figs,file="results/reads_figs.RData")

# -----------------------------------------------------------------------------
# Power analysis based on READS data ----
# -----------------------------------------------------------------------------

cor( dat$Yhat, dat$Yobs )^2

R2s = estimate_R2( dat,  Yhat = "Yhat", Z = "Z", strata = "strata" )
R2s
# Overall explanatory power within strata
1 - (1-R2s$R2_strat)*(1-R2s$R2_ML)
sd( dat$Yobs[dat$Z==0] )
#sqrt( pp ) * sd( dat$Yobs[dat$Z==0] )

target_MDES = 0.15
pp <- plot_power_curve( nrow(dat), mean(dat$Z), 1,
                        R2_ML = R2s$R2_ML, R2_strat=R2s$R2_strat )+
    scale_y_continuous(breaks=seq(0,0.3,by=0.05))+expand_limits( y = 0 ) 


tbl <- attr( pp, "table" )
tbl
crosses = tbl %>%
  group_by( Method ) %>%
  filter( abs(MDES-target_MDES) == min(abs(MDES-target_MDES)))
crosses
crosses$h[[1]] / crosses$h[[2]]

pp <- pp +
  geom_hline( yintercept = target_MDES, lty=3 ) +
  geom_vline( xintercept = crosses$h )
pp
ggsave("figures/READS_power_curve.pdf", pp,  width = 6, height = 4)


# Sanity check: SE similar to what we actually found empirically?
ff <- filter( tbl, h == 0.3 )
ff
mm <- multi_estimates %>%
  group_by( method ) %>%
  summarise( mean_SE = sqrt( mean( SE^2 ) ) ) %>%
  mutate( pow_calc = c( NA, NA, ff$SE[[1]], NA, ff$SE[[2]] ) * sd( dat$Yobs[dat$Z==0] ),
          ratio = mean_SE / pow_calc )
mm





# -----------------------------------------------------------------------------
# Retrospective Simulation Study ----
# This is to see what would have happened in READS study if they had coded only a subsample using each method

# Note: this mirrors our simulations and is just a sanity check that everything works as expected
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Oracle Analysis (Using Full Data) ----
# -----------------------------------------------------------------------------

# Oracle estimator - using all data
oracle <- est_oracle(dat)
tau_oracle <- oracle$est
SE_oracle <- oracle$SE




# Simulation parameters
n_iter <- 1000
coding_budgets <- seq(0.1, 0.9, by = 0.1)

cores <- min(parallel::detectCores() - 2, 9)
cl <- makeCluster(cores)
registerDoParallel(cl)

## Step 4.2: Run Simulation Across Budgets
set.seed(12345)

system.time({
  results_empirical <- foreach(
    h = coding_budgets,
    .combine = 'rbind',
    .packages = c("dplyr", "stratsampling")
  ) %dorng% {
    
    n_code <- round(h * nrow(dat)/2)*2
    
    # Run 1000 iterations for this coding budget
    res <- plyr::rdply(n_iter, one_iter_all(dat, n_code = n_code, finite_pop = T))
    res$h <- h
    res$n_code <- n_code
    
    return(res)
  }
})

stopCluster(cl)

# Save results
saveRDS(results_empirical, "results/READS_empirical_results.rds")

# Summarize across iterations
performance_summary <- results_empirical %>%
  mutate(
    covers = (est - 1.96*SEhat <= tau_oracle) & (est + 1.96*SEhat >= tau_oracle)
  ) %>%
  group_by(method, h, n_code) %>%
  summarise(
    # Bias
    mean_est = mean(est),
    bias = mean(est) - tau_oracle,
    
    # Variance/SE
    emp_SE = sd(est),
    mean_SEhat = mean(SEhat),
    
    # MSE
    MSE = mean((est - tau_oracle)^2),
    
    # Coverage
    coverage = mean(covers),
    
    # Number of iterations
    n_iter = n(),
    
    .groups = "drop"
  )

print(performance_summary)

# Add oracle comparison
performance_summary <- performance_summary %>%
  mutate(
    # Variance inflation vs. oracle
    var_ratio_oracle = emp_SE^2 / SE_oracle^2,
    
    # Variance reduction vs. SRS
    var_ratio_srs = emp_SE^2 / emp_SE[method == "SRS"]^2,
    pct_reduction_srs = (1 - var_ratio_srs) * 100
  )

## Visualize Simulation Results
# Empirical SE by method and budget
p_emp_se <- performance_summary %>%
  ggplot(aes(x = h, y = emp_SE, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +    
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_hline(yintercept = SE_oracle, linetype = "dashed") +
  #annotate("text", x = 0.5, y = SE_oracle * 0.95, 
  #     label = "Oracle (Full Coding)", size = 3) +
  scale_x_continuous(breaks=seq(0.1,0.9,by=0.1),limits=c(0.08,0.92),labels = scales::percent) +
  scale_y_continuous(limits = c(0, 0.1), breaks=seq(0,0.1,by=0.02)) +
  
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Empirical Standard Errors: Simulation Results",
    x = "Proportion Coded",
    y = "Empirical SE"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Variance reduction vs SRS
p_var_reduction <- performance_summary %>%
  filter(!method %in% c("SRS","Subset")) %>%
  ggplot(aes(x = h, y = pct_reduction_srs/100, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks=seq(0.1,0.9,by=0.1), limits=c(0.08,0.92),labels = scales::percent) +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,by=0.20), labels = scales::percent) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Variance Reduction Relative to SRS",
    x = "Proportion Coded",
    y = "% Variance Reduction"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Coverage rates
p_coverage <- performance_summary %>%
  ggplot(aes(x = h, y = coverage, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(limits = c(0.90, 1.0), labels = scales::percent) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Coverage of 95% Confidence Intervals",
    x = "Proportion Coded",
    y = "Coverage Rate"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
gridExtra::grid.arrange(p_emp_se, p_var_reduction,p_coverage, nrow = 3)
# Combined plot
p_combined <- gridExtra::grid.arrange(p_emp_se, p_var_reduction, nrow = 2)
print(p_combined)
ggsave("figures/READS_empirical_simulation_frontiers.pdf", p_combined,
       width = 10, height = 12)

# -----------------------------------------------------------------------------
# Detailed Analysis at Target Budget (How do simulations perform specifically at 30% coding?) ----
# -----------------------------------------------------------------------------
# Filter to h = 0.30
c=unique(results_empirical$h)[3]
results_30 <- results_empirical %>% filter(h == c)

# Summary table
table_30 <- results_30 %>%
  group_by(method) %>%
  summarise(
    N_coded = unique(n_code),
    Mean_Est = mean(est),
    Emp_SE = sd(est),
    Mean_SEhat = mean(SEhat),
    Coverage = mean((est - 1.96*SEhat <= tau_oracle) & 
                      (est + 1.96*SEhat >= tau_oracle)),
    .groups = "drop"
  ) %>%
  mutate(
    Var_Ratio_Oracle = (Emp_SE / SE_oracle)^2,
    Var_Ratio_SRS = (Emp_SE / Emp_SE[method == "SRS"])^2,
    Pct_Reduction = (1 - Var_Ratio_SRS) * 100
  )

print(table_30)


