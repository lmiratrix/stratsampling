# Simulation results: figures & tables (main factorial grid, 01_run_full_sim.R)


# Packages
library(tidyverse)
library(janitor)
library(broom)
library(scales)
library(kableExtra)

source("scripts/plotting_functions.R")   # method registry + plot helpers (single source of truth)

# ------------------------------- #
# 0) Load data and rebuild analysis_results ----
# ------------------------------- #
results1=readRDS("results/sim_results_balanced_exact_all_fast_v2.rds")
results2=readRDS("results/sim_results_balanced_approx_all_fast_v2.rds")
results3=readRDS("results/sim_results_unbalanced_all_fast_v2.rds")

results=rbind(results1,results2,results3)
rm(results1,results2,results3)

cell_keys <- c("N","n_strat","strata_config","tau","bias","residual_var","R2_target","pop.sd","n_code","coding_frac")
PF_MAIN <- if ("pilot_frac" %in% names(results)) sort(unique(na.omit(results$pilot_frac))) else numeric(0)

analysis_results <- cell_summary(
    results %>% mutate(coding_frac = n_code / N),
    cell_keys = cell_keys
)




analysis_results <- analysis_results %>%
    mutate(
        method        = factor(method, levels = intersect(method_levels, unique(method))),
        residual_var  = fct_relevel(factor(residual_var), "homogeneous", "heterogeneous", "extreme-contrast"),
        bias          = fct_relevel(factor(bias), "none", "small", "moderate", "large", "extreme-contrast") #,
        #strata_config = fct_relevel(factor(strata_config), "balanced-exact", "balanced-approx", "unbalanced")
    )

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


# Join oracle for variance inflation vs full coding
res.oracle <- readRDS("results/sim_results_oracle_fast_v2.rds")

results_oracle = res.oracle %>%
    group_by(N, n_strat, strata_config, tau, bias, residual_var, R2_target, pop.sd, method) %>%   
    mutate(
        CI_lower = est - Z_CRIT * SE,
        CI_upper = est + Z_CRIT * SE,
        covers   = if_else(CI_lower < tau & CI_upper > tau, 1L, 0L)
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

# Main-text filter: K=4, N=1000, sd_y=3 (one row per cell, i.e. per strata configuration)
main <- dat %>%
    filter(N == 1000, n_strat == 4, pop.sd == 3) %>%
    mutate(
        residual_var  = fct_relevel(residual_var, "homogeneous", "heterogeneous", "extreme-contrast"),
        bias          = fct_relevel(bias, "none","small", "moderate", "large", "extreme-contrast")
    )

# `main` averaged over the three strata configurations (average the variances, then divide).  
# Used by the figures that draw one line per method within each bias x residual-var
# panel; with `main` itself those would now get three points per budget.
main_agg <- main %>%
    group_by(N, n_strat, tau, bias, residual_var, R2_target, pop.sd, n_code, coding_frac, method) %>%
    summarise(
        emp_var_SRS  = mean(emp_var_SRS),
        emp_var_full = mean(emp_var_full),
        bias_est     = mean(bias_est),
        coverage     = mean(coverage),
        mean_sehat   = mean(mean_sehat),
        emp_var      = mean(emp_var),
        .groups = "drop") %>%
    mutate(emp_se             = sqrt(emp_var),
           pct_var_red_vs_SRS = 100 * (1 - emp_var / emp_var_SRS),
           re_vs_SRS          = emp_var / emp_var_SRS,
           re_vs_fullcoding   = emp_var / emp_var_full)


# MAIN TEXT FIGURES ---- 

# ------------------------------- #
# 1) Figure 1 — Efficiency frontiers (Empirical SE) ----
#    Two slices: R2 = 0.85 and R2 = 0.40; small-multiples by residual variance.
# ------------------------------- #

res <- analysis_results %>%
    mutate(any_bias  = ifelse(bias == "none", 0, 1),
           any_hetero = ifelse(residual_var == "homogeneous", 0, 1)) %>%
    group_by(N, n_strat, tau, pop.sd, any_bias, any_hetero, R2_target, n_code, coding_frac, method) %>%
    summarize(emp_var = mean(emp_var, na.rm = TRUE), .groups = "drop") %>%
    mutate(emp_se = sqrt(emp_var),
           any_bias = factor(any_bias, levels=c(0,1), labels=c("No Bias", "Biased")),
           any_hetero=factor(any_hetero, levels=c(0,1), labels=c("Homogeneous Var.", "Heterogeneous Var.")))

results_oracle2 = res.oracle %>%
    mutate(any_bias=ifelse(bias=="none",0,1),
           any_hetero=ifelse(residual_var=="homogeneous",0,1)) %>%
    group_by(N, n_strat, strata_config, tau, bias, residual_var, any_bias, any_hetero, R2_target, pop.sd, method) %>%
    summarise(emp_var = var(est), mean_est = mean(est),
              bias_est = mean(est - tau), .groups = "drop") %>%
    group_by(N, n_strat, tau, any_bias, any_hetero, R2_target, pop.sd, method) %>%
    summarise(
        mean_est   = mean(mean_est),
        bias_est   = mean(bias_est),
        emp_var    = mean(emp_var),
        .groups = "drop"
    ) %>%
    mutate(emp_se = sqrt(emp_var),
           any_bias = factor(any_bias, levels=c(0,1), labels=c("No Bias", "Biased")),
           any_hetero=factor(any_hetero, levels=c(0,1), labels=c("Homogeneous Var.", "Heterogeneous Var.")))

results_oracle$method="Full coding"


full_baseline2 <- results_oracle2 %>%
    select(N, n_strat,  tau, any_bias, any_hetero, R2_target, pop.sd,emp_var_full = emp_var, emp_se_full = emp_se)

res <- res %>%
    left_join(full_baseline2, by = setdiff(names(res),c("method", "n_code","coding_frac","method","emp_se","emp_var"))) %>%
    mutate(re_vs_fullcoding = emp_var / emp_var_full)


res$bias.var.config = sprintf("%s, %s", res$any_bias, res$any_hetero)

fig1 =  make_frontier_plot(res,facet=as.formula(R2_target~any_bias+any_hetero),oracle_ref=T,smooth=F,
                           methods = c("Subset", "SRS", "Stratified", "Stratified-Yhat", "Optimal"))+ 
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),legend.text = element_text(size=10))+
    scale_y_continuous(breaks=seq(0.2,0.5,by=0.1),limits=c(0.19,0.51), oob=scales::oob_squish)+
    method_scales(use_labels = F) 
print(fig1)

ggsave("figures/fig1_emp_SE_aggregated.pdf", fig1, width = 8.5, height = 6.5)



# ------------------------------- #
# 2) Figure 2 — Incremental value of stratification ----
#    ∆(% variance reduction vs SRS) across budgets and moderators
# ------------------------------- #

# Collapse over levels of bias (any vs. none) and residual variance (hetero vs. homogeneous)
base2 <- analysis_results %>%
    mutate(any_bias = factor(bias != "none", levels = c(FALSE, TRUE),
                             labels = c("No Bias", "Biased")),
           any_het  = factor(residual_var != "homogeneous", levels = c(FALSE, TRUE),
                             labels = c("Homogeneous Var.", "Heterogeneous Var."))) %>%
    agg_vr(ref = "SRS",
           by  = c("any_bias", "any_het", "R2_target", "coding_frac")) %>%
    rename(bias = any_bias, residual_var = any_het,
           pct_var_red_vs_SRS = pct_red, re_vs_SRS = vr)

fig2=make_delta_plot(base2, facet=as.formula("R2_target ~ bias+residual_var"),
                          methods=c("Stratified","Stratified-Yhat","Stratified-Random", "Optimal","Optimal-Permuted","Optimal-Reversed"))+ 
    scale_y_continuous(breaks=seq(-160,60,by=40), limits=c(-120,60), oob=scales::oob_squish)+   
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11), legend.text = element_text(size=10))+
    method_scales(use_labels = F,nrow=1)

print(fig2)
ggsave("figures/fig2_delta_opt_over_strat.pdf", fig2, width = 8.5, height = 6.5)




# SUPPLEMENTAL FIGURES ---- 
# ------------------------------- #
# 1) Figure 1 — Efficiency frontiers (Empirical SE) ----
#    Two slices: R2 = 0.85 and R2 = 0.40; small-multiples by residual variance.
# ------------------------------- #

fig1_high <- make_frontier_plot(main_agg, r2_val = 0.85, oracle_ref = F,
                                methods = c("Subset", "SRS", "Stratified", "Stratified-Yhat", "Optimal"))+scale_y_continuous(breaks=seq(0.2,0.5,by=0.1),limits=c(0.19,0.51), oob=scales::oob_squish)+
    xlab("Proportion coded")
fig1_low  <- make_frontier_plot(main_agg, r2_val = 0.40, oracle_ref=F,
                                methods = c("Subset", "SRS", "Stratified", "Stratified-Yhat", "Optimal"))+scale_y_continuous(breaks=seq(0.2,0.5,by=0.1),limits=c(0.19,0.51), oob=scales::oob_squish)

ggsave("figures/supplemental/fig1_frontier_highR2.pdf", fig1_high, width = 8.5, height = 6.5)
ggsave("figures/supplemental/fig1_frontier_lowR2.pdf",  fig1_low,  width = 8.5, height = 6.5)

# ------------------------------- #
# 2) Figure 2 — Incremental value of stratification ----
#    ∆(% variance reduction vs SRS) across budgets and moderators
# ------------------------------- #

fig2_high <- make_delta_plot(main_agg, r2_val = 0.85, bias_vals=c("none","small", "moderate","large"),
                             methods=c("Stratified","Stratified-Yhat","Stratified-Random","Optimal","Optimal-Permuted","Optimal-Reversed"))+ 
    method_scales(use_labels = TRUE, nrow = 2) +scale_y_continuous(breaks=seq(-100,60,by=20),limits=c(-60,60),oob=scales::oob_squish)
fig2_low  <- make_delta_plot(main_agg, r2_val = 0.40, bias_vals=c("none","small", "moderate","large"),
                             methods=c("Stratified","Stratified-Yhat","Stratified-Random","Optimal","Optimal-Permuted","Optimal-Reversed"))+
    method_scales(use_labels = TRUE, nrow = 2) +scale_y_continuous(breaks=seq(-100,80,by=20),limits=c(-100,80),oob=scales::oob_squish)

ggsave("figures/supplemental/fig2_delta_opt_over_strat_highR2.pdf", fig2_high, width = 7.2, height = 6.0)
ggsave("figures/supplemental/fig2_delta_opt_over_strat_lowR2.pdf",  fig2_low,  width = 7.2, height = 6.0)



fig2b=make_delta_plot(base2, facet=as.formula("R2_target ~ bias+residual_var")) + 
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11), legend.text = element_text(size=10))+
    scale_y_continuous(breaks =c(-10,0,10,20,30,40,50),oob=scales::squish,
                       labels = percent_format(scale = 1),limits = c(-10, 60))

print(fig2b)
ggsave("figures/supplemental/fig2_orig_delta_opt_over_strat.pdf", fig2b, width = 8.5, height = 6.5)




# ------------------------------- #
# 3) Figure 3 — Variance inflation vs FULL CODING  ----
#     y = 100 * Var(method) / Var(full-coding). 100% = as good as full coding.
# ------------------------------- #

fig3_high <- make_full_ratio_plot(main_agg, r2_val=0.85,methods = c("SRS", "Stratified", "Stratified-Yhat", "Optimal"))+
    scale_x_log10(n.breaks=8)+xlab("Proportion coded (log scale)")
fig3_low  <- make_full_ratio_plot(main_agg, r2_val=0.40,methods = c( "SRS", "Stratified", "Stratified-Yhat", "Optimal"))+
    scale_x_log10(n.breaks=8)+xlab("Proportion coded (log scale)")

ggsave("figures/supplemental/fig3_vs_full_highR2.pdf", fig3_high, width = 8.5, height = 6.5)
ggsave("figures/supplemental/fig3_vs_full_lowR2.pdf",  fig3_low,  width = 8.5, height = 6.5)


# 4) Figure 4 - Variance inflation vs. full coding ----

fig4<- make_full_ratio_plot(res,facet=as.formula(R2_target~any_bias+any_hetero),
                            methods = c("Subset", "SRS", "Stratified", "Stratified-Yhat", "Optimal")   )+ theme(panel.grid.minor = element_blank())+
    scale_y_continuous(breaks=seq(2,10,by=2),limits=c(1,10), oob=scales::oob_squish)+   
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),legend.text = element_text(size=10))+
    scale_x_log10()+xlab("Proportion coded (log scale)")


print(fig4)
ggsave("figures/supplemental/fig4_var_ratio_aggregated.pdf", fig4, width = 8.5, height = 6.5)



fig5<- make_full_ratio_plot(res,facet=as.formula(R2_target~any_bias+any_hetero),
                                 methods=c("SRS","Stratified","Stratified-Yhat","Optimal","Optimal-Reversed","Stratified-Random"))+ theme(panel.grid.minor = element_blank())+
    scale_y_continuous(breaks=seq(1,6,by=1),limits=c(0.9,8), oob=scales::oob_squish)+   
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),legend.text = element_text(size=10))+
    method_scales(use_labels = TRUE, nrow = 2) 

print(fig5)

ggsave("figures/supplemental/fig4b_VR_aggregated_with_wrong_optimal.pdf", fig5, width = 8.5, height = 6.5)



# ------------------------------- #
# 6) Pilot-informed designs on the main grid ----
#    All five pilot arms are computed on the same populations as the
#    baselines, at the single pilot share PF_MAIN (share of the coding
#    budget; the pilot counts toward the budget).  Case 2 arms use the
#    KNOWN strata and a stratified pilot to estimate variances; Case 1
#    arms draw an SRS pilot, SELECT the stratification of Yhat from it
#    (tertiles / quartiles / quintiles cut within arm), estimate the
#    variances on the selected strata (raw or EB-shrunk), and discard
#    the pilot before drawing the estimation sample.
#    The main grid is the only simulation with an SRS arm, so this is
#    where "is a pilot-selected stratification better than no
#    stratification?" is answered; 03c covers the pilot-share tradeoff.
# ------------------------------- #
pf_txt <- sprintf("%.0f%%", 100 * PF_MAIN[1])
HEADLINE_H <- 0.2
h_near <- function(x, target, tol = 1e-6) abs(x - target) < tol

pilot_show <- c("SRS", "Stratified", "Pilot-Selected-Prop",
                "Optimal", "Pilot-Selected-Optimal", "Pilot-Selected-Optimal-Shrink",
                "Pilot-Allocated-Pooled","Pilot-Allocated-Discard","Pilot-Selected-Optimal")

# 6a) Efficiency frontier (empirical SE), aggregated over bias / heterogeneity
#     classes as in Figure 5, with the full-coding reference line
fig6a <- suppressMessages(
    make_frontier_plot(res, facet = as.formula(R2_target ~ any_bias + any_hetero),
                       oracle_ref = TRUE, smooth = FALSE,
                       methods = pilot_show) +
        method_scales(use_labels = TRUE, nrow = 3) +
        labs(caption = sprintf("Pilot arms use a pilot of %s of the coding budget. Dashed grey: full coding.", pf_txt)) +
        theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),
              legend.text = element_text(size = 9)))
print(fig6a)
ggsave("figures/supplemental/fig6_pilot_frontier_aggregated.pdf", fig6a, width = 9, height = 7)

# 6b) Incremental value vs SRS (% variance reduction), aggregated classes; pilot-based allocation
fig6b = make_pilot_plot(base2, refs = c("Stratified", "Optimal"),
                pilot = c("Pilot-Allocated-Pooled", "Pilot-Allocated-Discard"),
                facet = R2_target ~ bias + residual_var)
print(fig6b)

ggsave("figures/fig6b_pilot_delta_vs_SRS_allocation.pdf", fig6b, width = 9, height = 6.5)


# 6c: pilot-based selection (same bias x variance classes)
fig_sel <- make_pilot_plot(base2, refs = c("Stratified", "Stratified-Yhat"),
                     pilot = c("Pilot-Selected-Prop", "Pilot-Selected-Optimal"),
                     facet = R2_target ~ bias + residual_var)
print(fig_sel)
ggsave("figures/fig6b_pilot_delta_vs_SRS_selection.pdf", fig_sel, width = 9, height = 6.5)



# 7) Validity: Monte Carlo bias and coverage of every arm (main grid) ----
validity6 <- main %>%
    group_by(Design = method) %>%
    summarise(`Mean bias` = mean(bias_est), `Max |bias|` = max(abs(bias_est)),
              `Mean coverage` = mean(coverage), `Min coverage` = min(coverage),
              `Max coverage` = max(coverage), 
              `Mean SE` = sqrt(mean(emp_var)),             # RMS SE (agg rule 2026-09-16)
              `Mean SEhat` = mean(mean_sehat), 
              `SE-hat / SE` = mean(mean_sehat / emp_se),   # per-cell calibration; averaging ratios is intended here
              .groups = "drop") %>%
    mutate(Design = method_labels[as.character(Design)])
print(validity6)
kbl(validity6, format = "latex", booktabs = TRUE, digits = 4,
    caption = sprintf("Monte Carlo bias and 95\\%% CI coverage by design, aggregated over all scenarios and budgets of the main simulation grid (three strata configurations pooled). Pilot arms use a pilot of %s of the coding budget.", pf_txt),
    label = "all_sims_validity_main") %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    readr::write_file("tables/all_sims_validity_main_grid.tex")



coverage_groups <- list(
    "core" = c(
        "Subset"     = "Subset",
        "SRS"        = "SRS",
        "Stratified-Yhat" = "Stratified ($\\hat{Y}_i$ quantile strata)",
        "Stratified" = "Stratified (proportional allocation)",
        "Optimal"    = "Optimal (oracle variances)"),
    "Misinformed designs" = c(
        "Stratified-Random" = "Stratified (random strata)",
        "Optimal-Permuted"  = "Optimal (permuted variances)",
        "Optimal-Reversed"  = "Optimal (reversed variances)"),
    "Pilot-based allocation" = c(
        "Pilot-Allocated-Pooled" = "Pilot-allocated (pilot pooled)",
        "Pilot-Allocated-Discard"    = "Pilot-allocated (pilot discarded)"),
    "Pilot-based selection" = c(
        "Pilot-Selected-Prop" = "Pilot-selected (proportional allocation)",
        "Pilot-Selected-Optimal"    = "Pilot-selected (optimal allocation)")
)

coverage_tab <- main %>%
    mutate(method = as.character(method)) %>%
    filter(method %in% unlist(lapply(coverage_groups, names))) %>%
    group_by(method) %>%
    summarise(cov_mean  = mean(coverage),
              cov_min   = min(coverage),
              cov_max   = max(coverage),
              bias_mean = mean(bias_est),
              bias_max  = max(abs(bias_est)),
              .groups = "drop")
print(coverage_tab)

pct_tex <- function(x) paste0(formatC(100 * x, format = "f", digits = 1), "\\%")
num_tex <- function(x) formatC(x, format = "f", digits = 4)

coverage_rows <- function(group_name) {
    labs <- coverage_groups[[group_name]]
    d    <- coverage_tab[match(names(labs), coverage_tab$method), ]
    stopifnot(!anyNA(d$method))   # every listed arm must be in the results
    rows <- sprintf("%s & %s & %s & %s & %s & %s \\\\ ",
                    labs, pct_tex(d$cov_mean), pct_tex(d$cov_min), pct_tex(d$cov_max),
                    num_tex(d$bias_mean), num_tex(d$bias_max))
    if (group_name == "core") return(rows)
    c("\\addlinespace",
      sprintf("\\multicolumn{6}{l}{\\underline{\\textit{%s}}} \\\\ ", group_name),
      rows)
}

c("\\begin{table}[h]",
  sprintf("\\caption{Coverage and Monte Carlo bias by estimation method (aggregated across all simulation configurations). Pilot-based designs use a pilot of %s\\%% of the coding budget.}",
          sub("%", "", pf_txt, fixed = TRUE)),
  "\\label{tab:coverage_bias_summary}",
  "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrrrr}",
  "\\toprule",
  " & \\multicolumn{3}{c}{95\\% CI Coverage} & \\multicolumn{2}{c}{Monte Carlo Bias} \\\\ ",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-6}",
  "Method & Mean & Min & Max & Mean & Max $|$Bias$|$ \\\\ ",
  "\\midrule\\addlinespace[2.5pt]",
  unlist(lapply(names(coverage_groups), coverage_rows)),
  "\\bottomrule",
  "\\end{tabular*}",
  "\\end{table}") %>%
    writeLines("tables/supplemental/coverage_main_sims.tex")

# 6d) Key numbers at the headline budget: every arm vs SRS and vs full coding

key6_base <- dat %>%
    filter(N == 1000, n_strat == 4, pop.sd == 3, h_near(coding_frac, HEADLINE_H))
key6 <- key6_base %>%
    group_by(Design = method, `Residual var.` = residual_var) %>%
    summarise(`Emp. SE`                   = sqrt(mean(emp_var)),
              `Var. ratio vs SRS`         = mean(emp_var) / mean(emp_var_SRS),
              `Var. ratio vs full coding` = mean(emp_var) / mean(emp_var_full),
              Coverage = mean(coverage), .groups = "drop") %>%
    mutate(Design = method_labels[as.character(Design)]) %>%
    arrange(`Residual var.`, `Var. ratio vs SRS`)
print(key6, n = 60)

if (FALSE){
    kbl(key6, format = "latex", booktabs = TRUE, digits = 3,
        caption = sprintf("All designs at a %.0f\\%% coding budget, by residual-variance pattern (ratio of the mean empirical variance to the mean empirical variance of the reference, averaged over bias patterns and surrogate accuracy; three strata configurations pooled). Pilot arms use a pilot of %s of the budget.", 100 * HEADLINE_H, pf_txt),
        label = "pilot_arms_h30_main") %>%
        kable_styling(latex_options = c("hold_position", "scale_down")) %>%
        readr::write_file("tables/supplemental/pilot_arms_h30_main_grid.tex")
}



# 6e) By strata configuration: does the price of learning the strata
#     depend on how the true strata are shaped?
by_config <- results %>%
    filter(!is.na(est), !is.na(SEhat), N == 1000, n_strat == 4, pop.sd == 3,
           method %in% pilot_show) %>%
    mutate(coding_frac = n_code / N) %>%
    group_by(strata_config, R2_target, residual_var, bias, coding_frac, method) %>%
    summarise(emp_var = var(est), .groups = "drop") %>%
    group_by(strata_config, R2_target, residual_var, coding_frac, method) %>%
    summarise(emp_se = sqrt(mean(emp_var)), .groups = "drop") %>%   # (agg rule 2026-09-16)
    mutate(method = factor(method, levels = method_levels))

fig6e <- ggplot(by_config %>% filter(R2_target == 0.4),
                aes(coding_frac, emp_se, color = method, shape = method, linetype = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_grid(strata_config ~ residual_var, labeller = labs_all) +
    method_scales(use_labels = TRUE, nrow = 3) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(1)) +
    labs(x = "Proportion coded (total budget, pilot included)", y = "Empirical Std. Error (avg. over bias patterns)",
         caption = sprintf("R\u00b2 target = 0.4. Pilot arms use a pilot of %s of the budget.", pf_txt)) +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), legend.text = element_text(size = 9))+
    guides(color=guide_legend(nrow=2))
print(fig6e)
ggsave("figures/supplemental/fig6g_pilot_by_strata_config_lowR2.pdf", fig6e, width = 10, height = 8)

fig6e_high <- fig6e %+% (by_config %>% filter(R2_target == 0.85)) +
    labs(caption = sprintf("R\u00b2 target = 0.85. Pilot arms use a pilot of %s of the budget.", pf_txt))
ggsave("figures/supplemental/fig6g_pilot_by_strata_config_highR2.pdf", fig6e_high, width = 10, height = 8)



## BAD DESIGNS PROPERTIES ACROSS SCENARIOS

source( "scripts/application_functions.R" )
cell_keys <- c("strata_config", "bias", "residual_var", "R2_target", "n_code")
designs <- c("Stratified", "Optimal", "Stratified-Yhat", "Stratified-Random", "Optimal-Permuted", "Optimal-Reversed", pilot_methods )

# ---- distribution of per-scenario variance ratios vs SRS ------------
cells <- cell_summary(results %>% filter(method %in% c("SRS", designs)), cell_keys)
srs   <- cells %>% filter(method == "SRS") %>% select(all_of(cell_keys), v_srs = emp_var)
risk0  <- cells %>% filter(method != "SRS") %>% left_join(srs, by = cell_keys) %>%
    mutate(vr = emp_var / v_srs,
           any_bias=ifelse(bias=="none",0,1), any_hetero=ifelse(residual_var=="homogeneous",0,1)) %>% group_by(method) %>%
    summarise(mean=mean(vr),
              median = median(vr),
              Q05=quantile(vr,0.05),  
              Q95 = quantile(vr,0.95), 
              n_scenarios = n(), .groups = "drop")

method_order = c("Oracle", "Subset", "SRS", "Stratified-Yhat", "Stratified", "Optimal","Stratified-Random",
                 "Optimal-Permuted","Optimal-Reversed","Optimal-Random","Pilot-Allocated-Pooled","Pilot-Allocated-Discard",
                 "Pilot-Selected-Prop","Pilot-Selected-Optimal")

risk = risk0 %>%
    mutate(method = factor(method,levels=method_order)) %>%
    filter(method %in% unlist(lapply(coverage_groups, names))) %>%
    group_by(Design=method)%>%
    mutate(Design = method_labels[as.character(Design)])%>%
    arrange(method)%>%
    select(Design, mean, median, Q05, Q95)%>%
    rename(Method=Design, Mean=mean, Median=median)


kbl(risk, format = "latex", booktabs = TRUE, digits = 3,
    caption = sprintf("Variance ratio relative to SRS across all simulation scenarios. Q05 and Q95 are 5th and 95th percentile values across simulation scenarios.")) %>%   
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    label_below_caption(label="sim_results_bad_designs") %>% 
    readr::write_file("tables/supplemental/bad_design_var_ratios.tex")



