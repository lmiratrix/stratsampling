# ===============================
# Simulation results: figures & tables
# ===============================

# Packages
library(tidyverse)
library(janitor)
library(broom)
library(scales)
library(gt)

# -------------------------------
# 0) Load data and rebuild analysis_results
# -------------------------------
# Path to the CSV you uploaded
results1=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_balanced_exact_strata.rds")
results2=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_balanced_approx_strata.rds")
results3=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_unbalanced_strata.rds")

results=rbind(results1,results2,results3)

results$method[results$method=="Stratified-Opt"]="Optimal"



analysis_results <- results %>%
    filter(!is.na(est), !is.na(SEhat)) %>%  
    group_by(N, n_strat,  tau, bias, residual_var, R2_target, pop.sd, n_code, method) %>%
    mutate(
        CI_lower = est - 2 * SEhat,
        CI_upper = est + 2 * SEhat,
        covers   = if_else(CI_lower < 0 & CI_upper > 0, 1L, 0L)
    ) %>%
    summarise(
        mean_est   = mean(est),
        bias_est   = mean(est - tau),
        coverage   = mean(covers),
        emp_var    = var(est),
        emp_se     = sd(est),
        mean_se    = mean(SE),
        mean_sehat = mean(SEhat),
        mse        = mean((est - tau)^2),
        .groups = "drop"
    )

# Harmonize factor levels and derive useful fields
analysis_results <- analysis_results %>%
    mutate(
        coding_frac   = n_code / N,
        method        = fct_relevel(method, "Subset", "SRS", "Stratified", "Optimal"),
        residual_var  = fct_relevel(factor(residual_var), "homogeneous", "heterogeneous", "extreme-contrast"),
        bias          = fct_relevel(factor(bias), "none", "small", "moderate", "large", "extreme-contrast") #,
        #strata_config = fct_relevel(factor(strata_config), "balanced-exact", "balanced-approx", "unbalanced")
    )

# Keys that index a "cell" apart from method
cell_keys <- c("N","n_strat","tau","bias","residual_var","R2_target","pop.sd","n_code","coding_frac")

# Join SRS baseline (for % reduction vs SRS)
srs_baseline <- analysis_results %>%
    filter(method == "SRS") %>%
    select(all_of(cell_keys), emp_var_SRS = emp_var, emp_se_SRS = emp_se)

dat <- analysis_results %>%
    left_join(srs_baseline, by = cell_keys) %>%
    mutate(
        pct_var_red_vs_SRS = 100 * (1 - emp_var / emp_var_SRS),  # % reduction in variance
        re_vs_SRS          = emp_var / emp_var_SRS               # variance ratio to SRS
    )

# Join "Full coding" (if present) for variance inflation vs full coding

res.oracle=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_oracle.rds")


results_oracle = res.oracle %>%
    group_by(N, n_strat, tau, bias, residual_var, R2_target, pop.sd, method) %>%
    mutate(
        CI_lower = est - 2 * SE,
        CI_upper = est + 2 * SE,
        covers   = if_else(CI_lower < 0 & CI_upper > 0, 1L, 0L)
    ) %>%
    summarise(
        mean_est   = mean(est),
        bias_est   = mean(est - tau),
        coverage   = mean(covers),
        emp_var    = var(est),
        emp_se     = sd(est),
        mean_se    = mean(SE),
        mean_sehat = mean(SE),
        mse        = mean((est - tau)^2),
        .groups = "drop"
    )


results_oracle$method="Full coding"


    full_baseline <- results_oracle %>%
        select(all_of(setdiff(cell_keys, c("coding_frac", "n_code"))), emp_var_full = emp_var, emp_se_full = emp_se)

    dat <- dat %>%
        left_join(full_baseline, by = setdiff(cell_keys, c("coding_frac","n_code"))) %>%
        mutate(re_vs_fullcoding = emp_var / emp_var_full)  # 100% = equal to full coding

# Main-text filter: balanced-exact, K=4, N=1000, sd_y=3
main <- dat %>%
    filter(N == 1000, n_strat == 4, pop.sd == 3) %>%
    mutate(
        residual_var  = fct_relevel(residual_var, "homogeneous", "heterogeneous", "extreme-contrast"),
        bias          = fct_relevel(bias, "none","small", "moderate", "large", "extreme-contrast")
    )

source("scripts/plotting_functions.R")

# -------------------------------
# 1) Figure 1 — Efficiency frontiers (Empirical SE)
#    Two slices: R2 = 0.85 and R2 = 0.40; small-multiples by residual variance.
# -------------------------------

fig1_high <- make_frontier_plot(main, r2_val = 0.85, oracle_ref = F)+scale_y_continuous(breaks=seq(0.2,0.5,by=0.1),limits=c(0.19,0.51), oob=scales::oob_squish)
fig1_low  <- make_frontier_plot(main, r2_val = 0.40, oracle_ref=F)+scale_y_continuous(breaks=seq(0.2,0.5,by=0.1),limits=c(0.19,0.51), oob=scales::oob_squish)

ggsave("figures/fig1_frontier_highR2.pdf", fig1_high, width = 8.5, height = 6.5)
ggsave("figures/fig1_frontier_lowR2.pdf",  fig1_low,  width = 8.5, height = 6.5)

# -------------------------------
# 2) Figure 2 — Incremental value of stratification
#    ∆(% variance reduction vs SRS) across budgets and moderators
# -------------------------------

fig2_high <- make_delta_plot(main, r2_val = 0.85, bias_vals=c("none","small", "moderate","large"))
fig2_low  <- make_delta_plot(main, r2_val = 0.40, bias_vals=c("none","small", "moderate","large"))

ggsave("figures/fig2_delta_opt_over_strat_highR2.pdf", fig2_high, width = 7.2, height = 6.0)
ggsave("figures/fig2_delta_opt_over_strat_lowR2.pdf",  fig2_low,  width = 7.2, height = 6.0)


# Make the same figure but collapse over levels of bias (any vs. none) and residual variance (hetero vs. homogeneous)
d2 = results %>%
    filter(!is.na(est), !is.na(SEhat)) %>%  
    mutate(bias=ifelse(bias=="none",0,1), residual_var=ifelse(residual_var=="homogeneous",0,1), coding_frac   = n_code / N) %>%
    group_by(N, n_strat,  tau, bias, residual_var, R2_target, pop.sd, n_code, coding_frac, method) %>%
    summarise(        mean_est   = mean(est),
                      bias_est   = mean(est - tau),
                      emp_var    = var(est),
                      emp_se     = sd(est),
                      mean_se    = mean(SE),
                      mean_sehat = mean(SEhat),
                      mse        = mean((est - tau)^2),
                      .groups = "drop")


d2_srs_baseline <- d2 %>%
    filter(method == "SRS") %>%
    select(all_of(cell_keys), emp_var_SRS = emp_var, emp_se_SRS = emp_se)

dat2 <- d2 %>%
    left_join(d2_srs_baseline, by = cell_keys) %>%
    mutate(
        pct_var_red_vs_SRS = 100 * (1 - emp_var / emp_var_SRS),  # % reduction in variance
        re_vs_SRS          = emp_var / emp_var_SRS               # variance ratio to SRS
    )


base2 <- dat2 %>%
    filter(method %in% c("Stratified","Optimal")) %>%
    select(all_of(cell_keys), method, pct_var_red_vs_SRS) %>%
    mutate(bias = factor(bias, levels=c(0,1), labels=c("No Bias", "Biased")),
           residual_var=factor(residual_var, levels=c(0,1), labels=c("Homogenous Var.", "Heterogeneous Var.")))

fig2b=make_delta_plot(base2, facet=as.formula("R2_target ~ bias+residual_var")) + theme(panel.grid.minor = element_blank(), 
                                                                                        strip.text = element_text(size = 11),
                                                                                        legend.text = element_text(size=10))
print(fig2b)
ggsave("figures/fig2b_delta_opt_over_strat.pdf", fig2b, width = 8.5, height = 6.5)



# -------------------------------
# 4) Figure 3 — Variance inflation vs FULL CODING (if present)
#     y = 100 * Var(method) / Var(full-coding). 100% = as good as full coding.
# -------------------------------

    fig4_high <- make_full_ratio_plot(main, r2_val=0.85)
    fig4_low  <- make_full_ratio_plot(main, r2_val=0.40)
    ggsave("figures/fig4_vs_full_highR2.pdf", fig4_high, width = 8.5, height = 6.5)
    ggsave("figures/fig4_vs_full_lowR2.pdf",  fig4_low,  width = 8.5, height = 6.5)


## Figure 5

res<- analysis_results %>%
    mutate(coding_frac=n_code/N,
           any_bias=ifelse(bias=="none",0,1),
           any_hetero=ifelse(residual_var=="homogeneous",0,1)) %>%
    group_by(N, n_strat, tau, pop.sd, any_bias, any_hetero, R2_target, n_code, coding_frac, method) %>%
    summarize(emp_se = mean(emp_se, na.rm=TRUE),
              emp_var = mean(emp_var, na.rm=TRUE),
              .groups="drop") %>%
    mutate(any_bias = factor(any_bias, levels=c(0,1), labels=c("No Bias", "Biased")),
           any_hetero=factor(any_hetero, levels=c(0,1), labels=c("Homogenous Var.", "Heterogeneous Var.")))


results_oracle2 = res.oracle %>%
    mutate(any_bias=ifelse(bias=="none",0,1),
           any_hetero=ifelse(residual_var=="homogeneous",0,1)) %>%
    group_by(N, n_strat, tau, any_bias, any_hetero, R2_target, pop.sd, method) %>%
    summarise(
        mean_est   = mean(est),
        bias_est   = mean(est - tau),
        emp_var    = var(est),
        emp_se     = sd(est),
        mean_se    = mean(SE),
        mean_sehat = mean(SE),
        mse        = mean((est - tau)^2),
        .groups = "drop"
    ) %>% 
    mutate(any_bias = factor(any_bias, levels=c(0,1), labels=c("No Bias", "Biased")),
           any_hetero=factor(any_hetero, levels=c(0,1), labels=c("Homogenous Var.", "Heterogeneous Var.")))




results_oracle$method="Full coding"


full_baseline2 <- results_oracle2 %>%
    select(N, n_strat,  tau, any_bias, any_hetero, R2_target, pop.sd,emp_var_full = emp_var, emp_se_full = emp_se)

res <- res %>%
    left_join(full_baseline2, by = setdiff(names(res),c("method", "n_code","coding_frac","method","emp_se","emp_var"))) %>%
    mutate(re_vs_fullcoding = emp_var / emp_var_full)


res$bias.var.config = sprintf("%s, %s", res$any_bias, res$any_hetero)

fig5 =  make_frontier_plot(res,facet=as.formula(R2_target~any_bias+any_hetero),oracle_ref=T,smooth=F)+ 
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),legend.text = element_text(size=10))+
    scale_y_continuous(breaks=seq(0.2,0.5,by=0.1),limits=c(0.19,0.51), oob=scales::oob_squish)
print(fig5)
ggsave("figures/fig5_emp_SE_aggregated.pdf", fig5, width = 8.5, height = 6.5)

fig5b<- make_full_ratio_plot(res,facet=as.formula(R2_target~any_bias+any_hetero))+ theme(panel.grid.minor = element_blank())+
    scale_y_continuous(breaks=seq(2,10,by=2),limits=c(1,11), oob=scales::oob_squish)+   
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),legend.text = element_text(size=10))

print(fig5b)

ggsave("figures/fig5_VR_aggregated.pdf", fig5b, width = 8.5, height = 6.5)



