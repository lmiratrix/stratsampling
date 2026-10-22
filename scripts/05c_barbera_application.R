
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
##  Data cleaning and summary measures ----

# Load the data
load("generated data/barbera-data-with-all-LLM-predictions.RData")

all.pred$ehat.mean=rowMeans(all.pred %>% select(starts_with("Yhat.conf_")))

dat = all.pred %>% 
    #filter(LLM_relevance=="Yes") %>% 
    #filter(n_coders >= 3) %>% 
    rename(score = Yobs)



stats = dat %>%
    summarise( mn = mean( score ),
               sd_Y = sd( score ),
               n = n(),
               sd_Yhat = sd(Yhat.mean))

avg_stats = stats %>%
    summarise( sd_bar_Y = sqrt( weighted.mean( sd_Y^2, w=n ) ),
               sd_bar_Yhat = sqrt(weighted.mean(sd_Yhat^2, w=n)))


dat = dat %>%  mutate(Yhat_raw=Yhat.mean, 
           score = score/avg_stats$sd_bar_Y,
           Yhat = Yhat.mean / avg_stats$sd_bar_Yhat)


dat <- dat %>%
    select(score, Yhat_raw, Yhat, ehat.mean,
           article_id, year, word_count, n_coders, n_distinct) %>%
    rename(Yobs = score,
           ehat = ehat.mean,
           n_distinct_Yhat=n_distinct) %>%
    mutate(ID = row_number())


# Basic descriptives
cat("Sample size:", nrow(dat), "\n")
cat("Mean outcome:", mean(dat$Yobs), "\n")

# Overall R²
overall_r = cor(dat$Yobs, dat$Yhat)
overall_r2 = cor(dat$Yobs, dat$Yhat)^2
SE_r = sqrt((1-overall_r^2)/(nrow(dat)-2))

cat(sprintf("Correlation: %.3f (SE=%.3f)\nOverall R² between human and LLM scores: %.3f\n", overall_r, SE_r, overall_r2))

# Calculate residuals
dat = dat %>% mutate(resid=Yobs-Yhat,
                     abs_resid=abs(resid),
                     sq_resid=resid^2)

# Summary statistics
residual_stats <- dat %>%
    summarise(
        mean_resid = mean(resid),
        sd_resid = sd(resid),
        mad_resid = mad(resid),
        min_resid = min(resid),
        max_resid = max(resid)
    )
print(residual_stats)

# Scatter plot of predictions vs. observed
p1 <- ggplot(dat, aes(x = Yhat, y = Yobs)) +
    geom_point(alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = "LLM Predictions vs. Human Scores",
         x = "LLM Prediction", y = "Human Score") +
    theme_bw()

# Residual distribution
p2 <- ggplot(dat, aes(x = resid)) +
    geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Distribution of Residuals",
         x = "Residual (Human - LLM)", y = "Count") +
    theme_bw()

# Residuals vs predictions (check for heteroskedasticity)
p3 <- ggplot(dat, aes(x = Yhat, y = abs_resid)) +
    geom_point(alpha = 0.3) +
    labs(title = "Absolute Residuals vs. Predictions",
         x = "LLM Prediction", y = "|Residual|") +
    theme_bw()

gridExtra::grid.arrange(p1, p2, p3, nrow = 3)
#ggsave("figures/barbera_empirical_prediction_quality.pdf",  width = 10, height = 12)

#  Explore Potential Stratification Variables ----

# Create categorical versions. Yhat is on the (standardized) 1-9 Barbera
# positivity scale; confidence ratings remain on 0-100 from the LLM call.
dat <- dat %>%
    mutate(
        length_cat = cut(word_count,
                         breaks = quantile(word_count,breaks=5),
                         include.lowest=T,
                         labels = c("Very Short", "Short", "Medium", "Long")),
        pred_cat = cut(Yhat,
                       breaks = c(-Inf, 3, 4.5, 5.6, Inf)/avg_stats$sd_bar_Yhat,
                       labels = c("Negative", "Slightly Negative",
                                  "Neutral/Slightly Positive", "Positive")),
        conf_cat = cut(ehat,
                       breaks=c(0, 78,100 ),
                       include.lowest=T,
                       labels=c("Low Confidence", "High Confidence"))
    )

## Examine Residual Patterns by Features

# By word count categories
resid_by_length <- dat %>%
    group_by(length_cat) %>%
    summarise(
        N = n(),
        mean_resid = mean(resid),
        var_resid = var(resid),
        sd_resid = sd(resid),
        .groups = "drop"
    ) %>% mutate(dim="length") %>% rename(cat=length_cat)
print(resid_by_length)

# By prediction level
resid_by_pred <- dat %>%
    group_by(pred_cat) %>%
    summarise(
        N = n(),
        mean_resid = mean(resid),
        var_resid = var(resid),
        sd_resid = sd(resid),
        .groups = "drop"
    ) %>% mutate(dim="Yhat") %>% rename(cat=pred_cat)
print(resid_by_pred)

# By confidence level
resid_by_conf <- dat %>%
    group_by(conf_cat) %>%
    summarise(
        N = n(),
        mean_resid = mean(resid),
        var_resid = var(resid),
        sd_resid = sd(resid),
        .groups = "drop"
    ) %>% mutate(dim="confidence") %>% rename(cat=conf_cat)
print(resid_by_conf)

resids_by_group=rbind(resid_by_length,resid_by_pred,resid_by_conf)


# Visualize
p_length <- ggplot(resid_by_length, aes(y = cat, x = var_resid)) +
    geom_col(position = "dodge") +
    labs(title = "Article Length", x= "Var(residual)",y="") +
    theme_minimal()

p_pred <- ggplot(resid_by_pred, aes(y = cat, x = var_resid)) +
    geom_col(position = "dodge") +
    labs(title = "Predicted Tone", x = "Var(residual)",y="") +
    theme_minimal()

p_conf <- ggplot(resid_by_conf, aes(y = cat, x= var_resid)) +
    geom_col(position = "dodge") +
    labs(title = "Confidence Level",x = "Var(residual)",y="") +
    theme_minimal()

gridExtra::grid.arrange(p_length, p_pred, p_conf, ncol= 3)

#ggsave("figures/empirical_residual_heterogeneity.pdf",
#  gridExtra::grid.arrange(p_length, p_pred, nrow = 2)#,
# width = 10, height = 8)

#  Create Candidate Stratifications ----

library(stratsampling)
dat$Yhat_binned = cut(dat$Yhat_raw, breaks=c(0, 2,3,4, 5, 6, 9)) # group rare very low/high scores with closest values
dat$all_agree = ifelse(dat$n_distinct_Yhat==1,1,0)


candidates = create_stratification_candidates( dat,
    vars   = c("Yhat", "ehat", "word_count",
               "n_distinct_Yhat", "all_agree"),
    custom=list(Yhat_binned=dat$Yhat_binned,
                Yhat_binned_all_agree = interaction(dat$Yhat_binned,dat$all_agree)),
    cross_custom = F)
    

#  Compare Candidates using Oracle Information ----
# Evaluate all candidates
comparison_oracle <- compare_stratifications(dat, candidates, Z=NULL, min_var_threshold = 0.01) # comparison using residuals (not feasible in practice)
comparison_precoding <- compare_stratifications_precoding(dat, candidates,Z=NULL) # comparison based on Yhat only (can be done before coding has started)


# Select top candidates meeting criteria
top_candidates <- comparison_oracle$summary %>%
    filter( min_stratum_size > 10, min_var_resid > 0.01) %>%
    arrange(desc(mean_spread)) %>%
    head(8)

cat("\n=== Top 5 Candidate Stratifications ===\n")
cat("Based on oracle comparisons (using known residuals Yobs-Yhat):")
print(top_candidates)

cat("Based on pre-coding comparisons (using only predictions Yhat): ")
head(comparison_precoding,n=8)


## Plot pre-coding proxy for stratification benefit vs. oracle
oracle_vs_precoding <- compare_oracle_vs_precoding(
    oracle_comparison = comparison_oracle,
    precoding_comparison = comparison_precoding,
    top_n = 55
)


# View scatter plot
print(oracle_vs_precoding$scatter_plot)
ggsave("figures/supplemental/barbera_precoding_strat_validation.pdf",oracle_vs_precoding$scatter_plot,  width = 8, height = 6)


# Compare variance decomposition oracle vs. precoding
results <- compare_decomposition_oracle(dat, candidates, n_code = floor(0.3*nrow(dat)),Z=NULL)
results$summary  # View table

decomp=plot_decomposition_comparison(results,show_components = F,top_n=30,selected="Yhat_binned_all_agree")
decomp
ggsave("figures/supplemental/barbera_precoding_strat_variance-reduction.pdf",decomp,  width = 8, height = 6)


## Select Final Stratification
# Apply best stratification
strat_name="Yhat_binned_all_agree"
dat$strata <- candidates[[strat_name]] # choose one that doesn't use residuals (oracle)

plot_strata_residuals(dat, dat$strata)
compute_variance_decomposition(dat, strata=dat$strata,Z=NULL)


# Document characteristics
source("scripts/application_functions.R") 
dat$agreement_level=factor(dat$all_agree,levels=0:1,labels=c("Some or All Disagree","All Agree"))
dat$yhat_label=factor(dat$Yhat_binned,labels=c("\\$<=$2","2-3", "3-4","4-5","5-6",">6"))
strata_characteristics <- dat %>%
    group_by( strata,  agreement_level, yhat_label) %>%
    summarise(
        N_zk = n(),
        mean_Yobs = mean(Yobs),
        sd_Yobs = sd(Yobs),
        mean_Yhat = mean(Yhat),
        mean_resid = mean(resid),
        var_resid = var(resid),
        sd_resid = sd(resid),
        mean_word_count = mean(word_count),
        .groups = "drop"
    ) 

print(strata_characteristics)


## Output Stratification Characteristics Table
strata_table <- strata_characteristics %>%
    mutate(
       k= row_number(),
        #`Pred. Tone`=Yhat,
        N_k = N_zk,
        `LLM Agreement` = agreement_level,
       `Predicted Score ($\\hat{Y}$)` = yhat_label,
        `Mean $Y$` = sprintf("%.2f", mean_Yobs),
        `Mean $\\hat{Y}$` = sprintf("%.2f", mean_Yhat),
        `Mean Resid.` = sprintf("%.2f", mean_resid),
        `Var(Resid.)` = sprintf("%.2f", var_resid),
        `Avg. Length` = sprintf("%.0f", mean_word_count)
    ) %>%
    select( k, `LLM Agreement`, `Predicted Score ($\\hat{Y}$)`, N_k, 
           `Mean Resid.`, `Var(Resid.)`) %>%
    rename(
            `$k$`      = k,
            `$N_{k}$` = N_k
    )



kbl(strata_table, format = "latex", booktabs = TRUE,align='llllcc',
    caption = "Characteristics of final stratification for the Barbera et al. (2021) NYT economic-news corpus determined by LLM agreement and predicted tone scores.",
    linesep="",escape=FALSE) %>%  kable_styling(latex_options = "hold_position") %>% 
    label_below_caption("barbera_strat") %>%
    write_file("tables/barbera_empirical_strata.tex")


# Estimates output from a single run across all methods (Example of practical real-world use) ----

source("scripts/sim_functions.R")
coding_fraction=0.3
n_code=round(coding_fraction*nrow(dat))

# Single sample at 30% coding budget
single_estimates <- compute_single_sample_estimates(dat, Z=NULL, coding_fraction = coding_fraction, finite_pop = FALSE,seed=2025,
                                                    pilot_frac=1/5, pilot_select_frac=1/5,
                                                    candidates=candidates, min_cell=1, min_pilot_cell=1)


comparison_table <- create_comparison_table(single_estimates)
print(comparison_table)

comparison_table=comparison_table %>% mutate(CI=sprintf("(%.3f, %.3f)", `95% CI Lower`, `95% CI Upper`), `Var. Inflation`=sprintf("%.2f", `Var. Inflation`)) %>%
    select(Method, `n Coded`, Estimate, SE, CI, `Var. Inflation`)


# Create forest plot
oracle <- est_oracle(dat,Z=NULL)

# Multiple samples to show variability

#multi_estimates <- compute_multiple_sample_estimates(dat, Z=NULL, coding_fraction = 0.30, finite_pop=F,
   #                                                  n_samples = 20, seed = 2026) %>% filter(!is.na(SE))

# Multiple samples to show variability
n_samp=100
multi_estimates <- compute_multiple_sample_estimates(dat, Z=NULL, coding_fraction =coding_fraction, finite_pop=F,
                                                     pilot_frac=1/5, pilot_select_frac=1/5, candidates=candidates,
                                                     n_samples = n_samp, seed = 2025, min_cell=1, min_pilot_cell=1)
multi_summary <- create_multi_sample_summary(multi_estimates)
print(multi_summary)


# Combined table: single sample (left) and repeated samples (right)
include_methods=c("Oracle","Subset","SRS","Stratified","Stratified-Yhat", "Optimal","Stratified-Random","Optimal-Permuted","Optimal-Reversed",
                  "Pilot-Allocated-Discard", "Pilot-Allocated-Pooled",
                  "Pilot-Selected-Prop","Pilot-Selected-Optimal")

barbera_tab <- write_application_table(
    single_estimates %>% filter(method %in% include_methods),
    multi_estimates  %>% filter(method %in% include_methods),
    file    = "tables/barbera-ests.tex",
    caption = sprintf("Results for the Barbera corpus estimating the mean economic tone of articles with a coding budget of $h=30\\%%$ ($n=%s$ coded documents; full coding uses all $N=%s$). 
                      Left: estimates from a single sample. Right: averages across %s repeated samples; coverage is the share of 95\\%% confidence intervals containing the full-coding estimate, and variance inflation is the estimated variance relative to that of full human coding.
                      Pilot-based designs use $h^{\\star}$=1/5 of the coding budget for the pilot.",
                      n_code, nrow(dat), n_samp),
    label   = "barbera_ests")
print(barbera_tab)

# Multi-sample forest plot
forest_multi <- create_multi_sample_forest_plot(multi_estimates,  tx=FALSE,
                                                oracle_est = oracle$est,jitter_seed=1234,
                                                include_methods = c("Subset","SRS","Stratified","Stratified-Yhat", "Optimal"))
print(forest_multi)
ggsave("figures/supplemental/barbera_emp_multi_forest.pdf",forest_multi,  width = 8, height = 6)







scatter=oracle_vs_precoding$scatter_plot

barbera_figs=list(decomp=decomp, forest_multi=forest_multi, scatter=scatter)
names(barbera_figs)=paste0("barbera_",names(barbera_figs))
save(barbera_figs,file="results/barbera_figs.RData")


# Retrospective Simulation Study ----
# This is to see what would have happened in the Barbera NYT coding effort if
# only a subsample had been coded using each method

# Note: this mirrors our simulations and is just a sanity check that everything works as expected

oracle <- est_oracle(dat,Z=NULL)
tau_oracle <- oracle$est
SE_oracle <- oracle$SE

source("scripts/sim_functions.R")

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

        n_code <- round(h * nrow(dat))

        # Run 1000 iterations for this coding budget
        res <- plyr::rdply(n_iter, one_iter_all(dat, n_code = n_code, Z=NULL, finite_pop = T))
        res$h <- h
        res$n_code <- n_code

        return(res)
    }
})

stopCluster(cl)

# Save results
saveRDS(results_empirical, "results/BARBERA_empirical_results.rds")

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
    scale_y_continuous(limits = c(0, 0.06), breaks=seq(0,0.06,by=0.01)) +

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

# Combined plot
gridExtra::grid.arrange(p_emp_se, p_var_reduction,p_coverage, nrow = 3)
p_combined <- gridExtra::grid.arrange(p_emp_se, p_var_reduction, nrow = 2)
print(p_combined)
ggsave("figures/BARBERA_empirical_simulation_frontiers.pdf", p_combined,
       width = 10, height = 12)

# Detailed Analysis at Target Budget ----
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

