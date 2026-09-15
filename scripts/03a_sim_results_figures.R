# ===============================
# Simulation results: figures & tables (main factorial grid, 01_run_full_sim.R)
#
# Sections 1-5 are the manuscript's main-grid figures (unchanged file
# names).  Section 6 is the pilot-informed arms that one_iter_all()
# adds when 01 is run with pilot_frac set --
#   Case 2 (strata known, pilot estimates variances):
#       Optimal-Pilot-Pooled, Optimal-Pilot-Ind
#   Case 1 (SRS pilot selects the strata, pilot discarded):
#       Stratified-PilotSelect, Optimal-PilotSelect, Optimal-PilotSelect-Shrink
# all at the single pilot share recorded in the results (column pilot_frac).
#
# Run from the stratsampling/ root.
# ===============================

# Packages
library(tidyverse)
library(janitor)
library(broom)
library(scales)
library(kableExtra)

source("scripts/plotting_functions.R")   # method registry + plot helpers (single source of truth)


# -------------------------------
# 0) Load data and rebuild analysis_results
# -------------------------------
results1=readRDS("results/sim_results_balanced_exact_all_fast.rds")
results2=readRDS("results/sim_results_balanced_approx_all_fast.rds")
results3=readRDS("results/sim_results_unbalanced_all_fast.rds")

results=rbind(results1,results2,results3)
rm(results1,results2,results3)

unknown <- setdiff(unique(results$method), method_levels)
if (length(unknown)) stop("Unknown method labels in results: ", paste(unknown, collapse = ", "))

have_pilot_arms <- all(pilot_methods %in% unique(results$method))
PF_MAIN <- if ("pilot_frac" %in% names(results)) sort(unique(na.omit(results$pilot_frac))) else numeric(0)
if (have_pilot_arms) {
    message("Pilot-informed arms present at pilot share(s): ", paste(PF_MAIN, collapse = ", "))
    if (length(PF_MAIN) != 1) warning("Expected a single pilot share in the main grid; found ", length(PF_MAIN))
} else {
    message("Pilot-informed arms not present (01 was run without pilot_frac); Section 6 will be skipped.")
}


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
        method        = factor(method, levels = intersect(method_levels, unique(method))),
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

res.oracle <- readRDS("results/sim_results_oracle_fast.rds")


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

# -------------------------------
# 1) Figure 1 — Efficiency frontiers (Empirical SE)
#    Two slices: R2 = 0.85 and R2 = 0.40; small-multiples by residual variance.
# -------------------------------

fig1_high <- make_frontier_plot(main, r2_val = 0.85, oracle_ref = F)+scale_y_continuous(breaks=seq(0.2,0.5,by=0.1),limits=c(0.19,0.51), oob=scales::oob_squish)
fig1_low  <- make_frontier_plot(main, r2_val = 0.40, oracle_ref=F)+scale_y_continuous(breaks=seq(0.2,0.5,by=0.1),limits=c(0.19,0.51), oob=scales::oob_squish)

ggsave("figures/supplemental/fig1_frontier_highR2.pdf", fig1_high, width = 8.5, height = 6.5)
ggsave("figures/supplemental/fig1_frontier_lowR2.pdf",  fig1_low,  width = 8.5, height = 6.5)

# -------------------------------
# 2) Figure 2 — Incremental value of stratification
#    ∆(% variance reduction vs SRS) across budgets and moderators
# -------------------------------

fig2_high <- make_delta_plot(main, r2_val = 0.85, bias_vals=c("none","small", "moderate","large"),methods=c("Stratified","Stratified-Random","Optimal","Optimal-Permuted","Optimal-Reversed"))
fig2_low  <- make_delta_plot(main, r2_val = 0.40, bias_vals=c("none","small", "moderate","large"),methods=c("Stratified","Stratified-Random","Optimal","Optimal-Permuted","Optimal-Reversed"))

ggsave("figures/supplemental/fig2_delta_opt_over_strat_highR2.pdf", fig2_high, width = 7.2, height = 6.0)
ggsave("figures/supplemental/fig2_delta_opt_over_strat_lowR2.pdf",  fig2_low,  width = 7.2, height = 6.0)




make_delta_plot(main, r2_val = 0.85, bias_vals=c("none","small", "moderate","large"),methods=c("Stratified","Stratified-Random","Optimal","Optimal-Permuted","Optimal-Random"))

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
    #filter(method %in% c("Stratified","Optimal")) %>%
    select(all_of(cell_keys), method, pct_var_red_vs_SRS) %>%
    mutate(bias = factor(bias, levels=c(0,1), labels=c("No Bias", "Biased")),
           residual_var=factor(residual_var, levels=c(0,1), labels=c("Homogenous Var.", "Heterogeneous Var.")))

fig2b=make_delta_plot(base2, facet=as.formula("R2_target ~ bias+residual_var")) + 
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11), legend.text = element_text(size=10))+
    scale_y_continuous(breaks = seq(-10, 100, by = 20),oob=scales::squish,
                       labels = percent_format(scale = 1),limits = c(-10, 85))

print(fig2b)
ggsave("figures/fig2_orig_delta_opt_over_strat.pdf", fig2b, width = 8.5, height = 6.5)


fig2b.new=make_delta_plot(base2, facet=as.formula("R2_target ~ bias+residual_var"),
                          methods=c("Stratified","Stratified-Random", "Optimal","Optimal-Random","Optimal-Reversed"))+ 
    scale_y_continuous(breaks=seq(-80,80,by=20), limits=c(-80,80),oob=scales::oob_squish)+
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11), legend.text = element_text(size=10))


print(fig2b.new)
ggsave("figures/fig2b_delta_opt_over_strat.pdf", fig2b.new, width = 8.5, height = 6.5)



# -------------------------------
# 4) Figure 3 — Variance inflation vs FULL CODING (if present)
#     y = 100 * Var(method) / Var(full-coding). 100% = as good as full coding.
# -------------------------------

fig4_high <- make_full_ratio_plot(main, r2_val=0.85)
fig4_low  <- make_full_ratio_plot(main, r2_val=0.40)
ggsave("figures/supplemental/fig4_vs_full_highR2.pdf", fig4_high, width = 8.5, height = 6.5)
ggsave("figures/supplemental/fig4_vs_full_lowR2.pdf",  fig4_low,  width = 8.5, height = 6.5)


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

ggsave("figures/supplemental/fig5_var_ratio_aggregated.pdf", fig5b, width = 8.5, height = 6.5)


fig5b.new<- make_full_ratio_plot(res,facet=as.formula(R2_target~any_bias+any_hetero),
                                 methods=c("SRS","Stratified","Optimal","Optimal-Reversed","Stratified-Random"))+ theme(panel.grid.minor = element_blank())+
    scale_y_continuous(breaks=seq(1,6,by=1),limits=c(0.9,7), oob=scales::oob_squish)+   
    theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),legend.text = element_text(size=10))

print(fig5b.new)

ggsave("figures/supplemental/fig5_VR_aggregated_with_wrong_optimal.pdf", fig5b.new, width = 8.5, height = 6.5)



# -------------------------------
# 6) Pilot-informed designs on the main grid
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
# -------------------------------

    pf_txt <- sprintf("%.0f%%", 100 * PF_MAIN[1])
    HEADLINE_H <- 0.3
    h_near <- function(x, target, tol = 1e-6) abs(x - target) < tol
    
    pilot_show <- c("SRS", "Stratified", "Stratified-PilotSelect",
                    "Optimal", "Optimal-PilotSelect", "Optimal-PilotSelect-Shrink",
                    "Optimal-Pilot-Pooled")
    
    # 6a) Efficiency frontier (empirical SE), aggregated over bias / heterogeneity
    #     classes as in Figure 5, with the full-coding reference line
    fig6a <- suppressMessages(
        make_frontier_plot(res, facet = as.formula(R2_target ~ any_bias + any_hetero),
                           oracle_ref = TRUE, smooth = FALSE,
                           methods = setdiff(pilot_show, "Optimal-PilotSelect")) +
            method_scales(use_labels = TRUE, nrow = 2) +
            labs(caption = sprintf("Pilot arms use a pilot of %s of the coding budget. Dashed grey: full coding.", pf_txt)) +
            theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),
                  legend.text = element_text(size = 9)))
    print(fig6a)
    ggsave("figures/fig6_pilot_frontier_aggregated.pdf", fig6a, width = 9, height = 7)
    
    # 6b) Incremental value vs SRS (% variance reduction), aggregated classes
    fig6b <- suppressMessages(
        make_delta_plot(base2, facet = as.formula("R2_target ~ bias+residual_var"),
                        methods = setdiff(pilot_show, "SRS")) +
            method_scales(use_labels = TRUE, nrow = 3) +
            scale_y_continuous(breaks = seq(-100, 100, by = 20), oob = scales::squish,
                               labels = percent_format(scale = 1), limits = c(-60, 85)) +
            ggtitle(NULL) +
            labs(caption = sprintf("Pilot arms use a pilot of %s of the coding budget; Case 1 arms discard it before estimation.", pf_txt)) +
            theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11),
                  legend.text = element_text(size = 9)))
    print(fig6b)
    ggsave("figures/fig6b_pilot_delta_vs_SRS.pdf", fig6b, width = 9, height = 7)
    
    # 6c) Paired variance ratios: what does each pilot design cost or gain
    #     relative to its natural reference?  Per simulation cell, then
    #     geometric mean over bias patterns.
    cellv <- main %>% select(all_of(cell_keys), method, emp_var)
    ratio_pair <- function(num, den, label) {
        cellv %>% filter(method == num) %>%
            inner_join(cellv %>% filter(method == den) %>%
                           select(all_of(cell_keys), emp_var_den = emp_var), by = cell_keys) %>%
            transmute(across(all_of(cell_keys)), comparison = label,
                      log_vr = log(emp_var / emp_var_den))
    }
    pairs6 <- bind_rows(
        ratio_pair("Stratified-PilotSelect", "SRS",
                   "Pilot-selected strata vs SRS (proportional)"),
        ratio_pair("Optimal-PilotSelect-Shrink", "SRS",
                   "Pilot-selected strata vs SRS (optimal, shrunk var.)"),
        ratio_pair("Stratified-PilotSelect", "Stratified",
                   "Proportional: pilot-selected vs known strata"),
        ratio_pair("Optimal-PilotSelect-Shrink", "Optimal",
                   "Optimal: pilot-selected (shrunk) vs oracle on known strata"),
        ratio_pair("Optimal-PilotSelect", "Optimal-PilotSelect-Shrink",
                   "Shrinkage: raw vs shrunk pilot variances"),
        ratio_pair("Optimal-PilotSelect-Shrink", "Optimal-Pilot-Pooled",
                   "Case 1 vs Case 2: selected (shrunk) vs pooled pilot on known strata")
    ) %>%
        mutate(comparison = factor(comparison, levels = unique(comparison)))
    
    pairs6_disp <- pairs6 %>%
        group_by(comparison, R2_target, residual_var, coding_frac) %>%
        summarise(vr = exp(mean(log_vr)), .groups = "drop") %>%
        mutate(R2_lab = factor(paste0("R\u00b2 = ", R2_target)))
    
    fig6c <- ggplot(pairs6_disp, aes(coding_frac, vr, color = R2_lab, group = R2_lab)) +
        geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
        geom_line(linewidth = 0.8) +
        geom_point(size = 1.6) +
        facet_grid(comparison ~ residual_var, scales = "free_y",
                   labeller = labeller(residual_var = c("homogeneous" = "Residual var: homogeneous",
                                                        "heterogeneous" = "Residual var: heterogeneous",
                                                        "extreme-contrast" = "Residual var: extreme-contrast"),
                                       comparison = label_wrap_gen(26))) +
        scale_color_manual(name = NULL, values = c("#1b9e77", "#7570b3")) +
        scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(1)) +
        labs(x = "Proportion coded (total budget, pilot included)",
             y = "Variance ratio (geometric mean over bias patterns)",
             caption = sprintf("Values below 1 favour the first-named design. Pilot share = %s of the budget.", pf_txt)) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank(),
              strip.text.y.right = element_text(angle = 0, size = 8, hjust = 0))
    print(fig6c)
    ggsave("figures/new-sims/fig6c_pilot_paired_ratios.pdf", fig6c, width = 10, height = 9.5)
    readr::write_csv(pairs6_disp, "tables/pilot_paired_ratios_main_grid.csv")
    
    # 6d) Validity: Monte Carlo bias and coverage of every arm (main grid)
    validity6 <- main %>%
        group_by(Design = method) %>%
        summarise(`Mean bias` = mean(bias_est), `Max |bias|` = max(abs(bias_est)),
                  `Mean coverage` = mean(coverage), `Min coverage` = min(coverage),
                  `Max coverage` = max(coverage), 
                  `Mean SE` = mean(emp_se),
                  `Mean SEhat` = mean(mean_sehat), 
                  `SE-hat / SE` = mean(mean_sehat / emp_se),
                  .groups = "drop") %>%
        mutate(Design = method_labels[as.character(Design)])
    print(validity6)
    kbl(validity6, format = "latex", booktabs = TRUE, digits = 4,
        caption = sprintf("Monte Carlo bias and 95\\%% CI coverage by design, aggregated over all scenarios and budgets of the main simulation grid (three strata configurations pooled). Pilot arms use a pilot of %s of the coding budget.", pf_txt),
        label = "all_sims_validity_main") %>%
        kable_styling(latex_options = c("hold_position", "scale_down")) %>%
        readr::write_file("tables/all_sims_validity_main_grid.tex")
    
    fig6d_low  <- suppressMessages(make_coverage_plot(main, r2_val = 0.40, methods = c("Stratified", pilot_methods)) +
                                       method_scales(use_labels = TRUE, nrow = 2) + ggtitle(NULL))
    fig6d_high <- suppressMessages(make_coverage_plot(main, r2_val = 0.85, methods = c("Stratified", pilot_methods)) +
                                       method_scales(use_labels = TRUE, nrow = 2) + ggtitle(NULL))
    ggsave("figures/pilot-sims/fig6d_pilot_coverage_lowR2.pdf",  fig6d_low,  width = 9, height = 8)
    ggsave("figures/pilot-sims/fig6d_pilot_coverage_highR2.pdf", fig6d_high, width = 9, height = 8)
    
    # 6e) Key numbers at the headline budget: every arm vs SRS and vs full coding
    key6 <- dat %>%
        filter(N == 1000, n_strat == 4, pop.sd == 3, h_near(coding_frac, HEADLINE_H)) %>%
        group_by(Design = method, `Residual var.` = residual_var) %>%
        summarise(`Emp. SE` = mean(emp_se),
                  `Var. ratio vs SRS` = exp(mean(log(re_vs_SRS))),
                  `Var. ratio vs full coding` = exp(mean(log(re_vs_fullcoding))),
                  Coverage = mean(coverage), .groups = "drop") %>%
        mutate(Design = method_labels[as.character(Design)]) %>%
        arrange(`Residual var.`, `Var. ratio vs SRS`)
    print(key6, n = 60)
    kbl(key6, format = "latex", booktabs = TRUE, digits = 3,
        caption = sprintf("All designs at a %.0f\\%% coding budget, by residual-variance pattern (geometric mean of variance ratios over bias patterns and surrogate accuracy; three strata configurations pooled). Pilot arms use a pilot of %s of the budget.", 100 * HEADLINE_H, pf_txt),
        label = "pilot_arms_h30_main") %>%
        kable_styling(latex_options = c("hold_position", "scale_down")) %>%
        readr::write_file("tables/supplemental/pilot_arms_h30_main_grid.tex")
    
    # 6f) Which stratification of Yhat does the selection pilot choose?
    if ("strat_choice" %in% names(results)) {
        choice6 <- results %>%
            filter(method == "Optimal-PilotSelect", !is.na(strat_choice)) %>%
            mutate(coding_frac = n_code / N,
                   strat_choice = factor(strat_choice,
                                         levels = c("Yhat_tertiles", "Yhat_quartiles", "Yhat_quintiles"),
                                         labels = c("Tertiles", "Quartiles", "Quintiles"))) %>%
            count(strata_config, R2_target, residual_var, bias, coding_frac, strat_choice) %>%
            group_by(strata_config, R2_target, residual_var, bias, coding_frac) %>%
            mutate(share = n / sum(n)) %>% ungroup()
        
        fig6f <- choice6 %>%
            group_by(strata_config, R2_target, coding_frac, strat_choice) %>%
            summarise(share = mean(share), .groups = "drop") %>%
            ggplot(aes(coding_frac, share, fill = strat_choice)) +
            geom_col(width = 0.08, color = "grey30", linewidth = 0.2) +
            facet_grid(R2_target ~ strata_config,
                       labeller = labeller(R2_target = c("0.4" = "R\u00b2 target = 0.4", "0.85" = "R\u00b2 target = 0.85"))) +
            scale_fill_manual(name = "Selected stratification of Yhat",
                              values = c("Tertiles" = "#c7e9c0", "Quartiles" = "#74c476", "Quintiles" = "#238b45")) +
            scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(1)) +
            scale_y_continuous(labels = percent_format(1)) +
            labs(x = "Proportion coded (total budget, pilot included)", y = "Share of iterations",
                 caption = sprintf("Selection pilot = %s of the budget; averaged over bias and residual-variance patterns.", pf_txt)) +
            theme_minimal() +
            theme(legend.position = "bottom", panel.grid.minor = element_blank())
        print(fig6f)
        ggsave("figures/pilot-sims/fig6f_pilot_select_choice.pdf", fig6f, width = 9, height = 5.5)
        
       
    }
    
    # 6g) By strata configuration: does the price of learning the strata
    #     depend on how the true strata are shaped?
    by_config <- results %>%
        filter(!is.na(est), !is.na(SEhat), N == 1000, n_strat == 4, pop.sd == 3,
               method %in% pilot_show) %>%
        mutate(coding_frac = n_code / N) %>%
        group_by(strata_config, R2_target, residual_var, bias, coding_frac, method) %>%
        summarise(emp_se = sd(est), .groups = "drop") %>%
        group_by(strata_config, R2_target, residual_var, coding_frac, method) %>%
        summarise(emp_se = mean(emp_se), .groups = "drop") %>%
        mutate(method = factor(method, levels = method_levels))
    
    fig6g <- ggplot(by_config %>% filter(R2_target == 0.4),
                    aes(coding_frac, emp_se, color = method, shape = method, linetype = method)) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 1.8) +
        facet_grid(strata_config ~ residual_var, labeller = labs_all) +
        method_scales(use_labels = TRUE, nrow = 3) +
        scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(1)) +
        labs(x = "Proportion coded (total budget, pilot included)", y = "Empirical Std. Error (avg. over bias patterns)",
             caption = sprintf("R\u00b2 target = 0.4. Pilot arms use a pilot of %s of the budget.", pf_txt)) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank(), legend.text = element_text(size = 9))
    print(fig6g)
    ggsave("figures/pilot-sims/fig6g_pilot_by_strata_config_lowR2.pdf", fig6g, width = 10, height = 8)
    
    fig6g_high <- fig6g %+% (by_config %>% filter(R2_target == 0.85)) +
        labs(caption = sprintf("R\u00b2 target = 0.85. Pilot arms use a pilot of %s of the budget.", pf_txt))
    ggsave("figures/pilot-sims/fig6g_pilot_by_strata_config_highR2.pdf", fig6g_high, width = 10, height = 8)
    

