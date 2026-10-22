# ===================================================================== #
# 03c_sim_pilot_results.R
#
# Figures and tables for the pilot-based design simulation
# (results/sim_results_pilot_balanced_exact_fast.rds, produced by
# 01b_run_sim_pilot_fast.R; identical to 01b_run_sim_pilot.R).
#
# Arms in the results file (one population per iteration, all arms on it):
#   Baselines (no pilot):  Stratified (proportional, known strata),
#                          Optimal (oracle Neyman), Optimal-Random,
#                          Optimal-Permuted, Optimal-Reversed
#   Case 2 (strata known, pilot estimates variances), per pilot share:
#                          Pilot-Allocated-Pooled, Pilot-Allocated-Discard
#   Case 1 (SRS pilot SELECTS the strata), per share:
#                          Pilot-Selected-Prop        (pilot discarded, proportional wave 2)
#                          Pilot-Selected-Optimal           (pilot discarded, Neyman wave 2 on
#                                                         EB-shrunk pilot variances)
#
# Sections 1-5 reproduce the outputs the manuscript already references
# (same file names).  Section 6 regenerates the fixed-share decision
# table.  Section 7 covers the Case 1 (pilot-selected stratification)
# arms, including the two-phase estimators and the calibration of their
# standard errors.
#
# Run from the stratsampling/ root.
# ===================================================================== #

library(tidyverse)
library(scales)
library(kableExtra)

source("scripts/plotting_functions.R")   # method registry (single source of truth)

# ---------------------------------------------------------------------#
# Settings ----
# ---------------------------------------------------------------------#
RESULTS_FILE   <- "results/sim_results_pilot_balanced_exact_fast_v2.rds"
HEADLINE_PF    <- 1/5      # pre-specified pilot share for headline tables
HEADLINE_H     <- 0.2      # pre-specified coding budget for headline tables
PF_SHOW        <- c(0.10, 0.25, 1/3, 0.50, 0.75)   # shares shown where 11 are too many
HH_SHOW        <- c(0.3, 0.5, 0.7)             # budgets shown in tradeoff figures
Z_CRIT         <- 2        # CI half-width multiplier used for coverage (as in prior runs)
MIN_PILOT_CELL <- 1        # per-cell pilot floor quoted in the text (annotation only)

pf_near <- function(x, target, tol = 1e-6) abs(x - target) < tol

# ---------------------------------------------------------------------#
# Load and sanity-check
# ---------------------------------------------------------------------#
results_pilot <- readRDS(RESULTS_FILE)
pilot_small=readRDS("results/sim_results_pilot_small.rds")
results_pilot=rbind(results_pilot, pilot_small)

PILOT_FRACS <- sort(unique(na.omit(results_pilot$pilot_frac)))
N_SIM   <- unique(results_pilot$N)[1]
N_CELLS <- unique(results_pilot$n_strat)[1] * 2          # true strata x arms

method_levels_pilot <- c("Stratified", "Optimal",  "Optimal-Permuted","Stratified-Yhat",
                         "Optimal-Reversed",
                         "Pilot-Allocated-Pooled", "Pilot-Allocated-Discard",
                         "Pilot-Selected-Prop", "Pilot-Selected-Optimal")

missing_arms <- setdiff(method_levels_pilot, unique(results_pilot$method))
if (length(missing_arms)) {
    message("Arms not present in the results file (figures using them are skipped): ",
            paste(missing_arms, collapse = ", "))
}
unknown_arms <- setdiff(unique(results_pilot$method), method_levels_pilot)
if (length(unknown_arms)) stop("Unknown method labels: ", paste(unknown_arms, collapse = ", "))
have_case1 <- all(pilot_case1_methods %in% unique(results_pilot$method))

# The two Neyman pilot-select rows (Pilot-Selected-Optimal and its two-phase
# version) are NA (not dropped) when the selection pilot was too thin to
# yield variance parameters (an arm with <= 1 pilot unit in every cell
# of the chosen stratification); the proportional rows are still valid.
# Count them before filtering.
thin_pilot <- results_pilot %>%
    filter(method == "Pilot-Selected-Optimal") %>%
    group_by(n_code, pilot_frac) %>%
    summarise(n_iter = n(), n_na = sum(is.na(est)), .groups = "drop") %>%
    filter(n_na > 0)
if (nrow(thin_pilot)) {
    message("Pilot-select iterations returned NA (pilot too thin):")
    print(thin_pilot)
}

# Palette for the pilot figures (subset of the registry, same colours)
pal_pilot    <- method_pal[method_levels_pilot]
shapes_pilot <- method_shapes[method_levels_pilot]
lty_pilot    <- method_linetypes[method_levels_pilot]

pf_label <- function(pf) {
    factor(sprintf("Pilot = %.0f%% of budget", 100 * pf),
           levels = sprintf("Pilot = %.0f%% of budget", 100 * sort(unique(pf))))
}
pf_pct <- function(pf) factor(percent(pf, 1), levels = percent(sort(unique(pf)), 1))

lab_rv   <- c("homogeneous"      = "Homogeneous residual var.",
              "heterogeneous"    = "Heterogeneous residual var.",
              "extreme-contrast" = "Extreme-contrast residual var.")
lab_bias <- c("none" = "Bias: none", "small" = "Bias: small", "moderate" = "Bias: moderate",
              "large" = "Bias: large", "extreme-contrast" = "Bias: extreme-contrast")
lab_R2   <- c("0.4" = "R\u00b2 target = 0.4", "0.85" = "R\u00b2 target = 0.85")
rv_lab <- labeller(residual_var = lab_rv, bias = lab_bias, R2_target = lab_R2,
                   any_bias = label_value, any_residual_var = label_value)

# ---------------------------------------------------------------------#
# Full per-cell summary (all methods, all shares) ----
# ---------------------------------------------------------------------#
pilot_summary0 <- results_pilot %>%
    filter(!is.na(est), !is.na(SEhat)) %>%
    mutate(coding_frac = n_code / N,
           CI_lower = est - Z_CRIT * SEhat,
           CI_upper = est + Z_CRIT * SEhat,
           covers = if_else(CI_lower < tau & CI_upper > tau, 1L, 0L)) %>%
    group_by(method, pilot_frac, bias, residual_var, R2_target,
             n_code, coding_frac) %>%
    summarise(
        mean_est   = mean(est),
        bias_est   = mean(est - tau),
        coverage   = mean(covers),
        emp_var    = var(est),
        emp_se     = sd(est),
        mean_se    = mean(SE),
        mean_sehat = mean(SEhat),
        mse        = mean((est - tau)^2),
        n_iter     = n(),
        .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = method_levels_pilot),
           bias = fct_relevel(factor(bias), "none", "small", "moderate", "large", "extreme-contrast"),
           residual_var = fct_relevel(factor(residual_var), "homogeneous",
                                      "heterogeneous", "extreme-contrast"),
           se_logvar = sqrt(2 / (n_iter - 1)),           # MC SE of log(emp_var)
           any_bias = factor(ifelse(bias == "none", 0, 1), levels = 0:1,
                             labels = c("No Bias", "Biased")),
           any_residual_var = factor(ifelse(residual_var == "homogeneous", 0, 1), levels = 0:1,
                                     labels = c("Homogeneous Var.", "Heterogeneous Var.")))

# Headline summary for the pilot simulations
pilot_summary <- pilot_summary0 %>%
    filter(!method %in% c( "Optimal-Random",
                           "Optimal-Permuted", "Optimal-Reversed"))

# Oracle Optimal reference per cell
opt_ref <- pilot_summary0 %>%
    filter(method == "Optimal") %>%
    select(bias, residual_var, R2_target, n_code, coding_frac, emp_var_opt = emp_var)

# Proportional (known strata) reference per cell
prop_ref <- pilot_summary0 %>%
    filter(method == "Stratified") %>%
    select(bias, residual_var, R2_target, n_code, coding_frac,
           emp_var_prop = emp_var, se_lv_prop = se_logvar)

# Per-cell log variance ratio of every arm vs the proportional (known
# strata) and oracle references, with MC intervals (independent-cell
# approximation; conservative for arms computed on the same population).
ratios_all <- pilot_summary0 %>%
    left_join(prop_ref, by = c("bias", "residual_var", "R2_target", "n_code", "coding_frac")) %>%
    left_join(opt_ref,  by = c("bias", "residual_var", "R2_target", "n_code", "coding_frac")) %>%
    mutate(log_vr_prop = log(emp_var / emp_var_prop),
           log_vr_opt  = log(emp_var / emp_var_opt),
           se_log      = sqrt(se_logvar^2 + se_lv_prop^2),
           pct_red     = 100 * (1 - exp(log_vr_prop)),
           pct_red_lo  = 100 * (1 - exp(log_vr_prop + 1.96 * se_log)),
           pct_red_hi  = 100 * (1 - exp(log_vr_prop - 1.96 * se_log)))

# ---------------------------------------------------------------------#
# 1) Case 2 frontier ---- 
# (aggregated over scenarios), pilot share 50%;
#    baselines replicated across facets 
# ---------------------------------------------------------------------#
plot_df <- bind_rows(
    pilot_summary %>% filter(is.na(pilot_frac)) %>% select(-pilot_frac) %>%
        crossing(pilot_frac = PILOT_FRACS),
    pilot_summary %>% filter(!is.na(pilot_frac))
) %>%
    mutate(pf_lab = pf_label(pilot_frac))

fig1_df <- plot_df %>%
    group_by(method, pilot_frac, R2_target, coding_frac, any_bias, any_residual_var) %>%
    summarise(emp_se = mean(emp_se), .groups = "drop")

fig1 <- ggplot(fig1_df %>% filter(is.na(pilot_frac) | pf_near(pilot_frac, 0.5)),
               aes(coding_frac, emp_se, color = method, shape = method, linetype = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_grid(R2_target ~ any_bias + any_residual_var, labeller = rv_lab) +
    scale_color_manual(name = "Method", values = pal_pilot) +
    scale_shape_manual(name = "Method", values = shapes_pilot) +
    scale_linetype_manual(name = "Method", values = lty_pilot) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(accuracy = 1)) +
    labs(x = "Proportion coded (total budget, pilot included)",
         y = "Empirical Std. Error (avg. across scenarios)") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())+
    method_scales(use_labels = TRUE)+
  guides(color = guide_legend(nrow = 3))

print(fig1)
ggsave("figures/pilot-sims/fig_pilot_frontier.pdf", fig1, width = 10, height = 6.5)

# ---------------------------------------------------------------------#
# 2) Is the pilot worthwhile? ----
# Pooled pilot and proportional relative to
#    oracle Optimal at the same budget
# ---------------------------------------------------------------------#
fig2_df <- pilot_summary %>%
    filter(method %in% c("Pilot-Allocated-Pooled", "Pilot-Allocated-Discard","Stratified")) %>%
    left_join(opt_ref, by = c("bias", "residual_var", "R2_target", "n_code", "coding_frac")) %>%
    mutate(vr_vs_opt = emp_var / emp_var_opt) %>%
    group_by(method, pilot_frac, any_bias, any_residual_var, R2_target, coding_frac) %>%
    summarise(vr_vs_opt = mean(vr_vs_opt), .groups = "drop") %>%
    mutate(pf_lab = as.factor(ifelse(is.na(pilot_frac), 0, round(pilot_frac, 3))))

fig2 <- ggplot(fig2_df %>% filter(is.na(pilot_frac) | pf_near(pilot_frac, 1/4) | pf_near(pilot_frac, 1/2)),
               aes(coding_frac, vr_vs_opt, color = pf_lab, linetype = method,
                   group = interaction(pf_lab, method))) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_grid(R2_target ~ any_bias + any_residual_var, labeller = rv_lab) +
    scale_linetype_manual(name = "Method",
                          values = c("Pilot-Allocated-Pooled" = "dashed","Pilot-Allocated-Discard" = "dotted" , "Stratified" = "solid")) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(accuracy = 1)) +
    labs(x = "Proportion coded", y = "Empirical variance relative to oracle Optimal",
         col = "Pilot fraction (0 = stratified)") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(fig2)
ggsave("figures/pilot-sims/fig_is_pilot_worthwhile.pdf", fig2, width = 10, height = 6.5)

# ---------------------------------------------------------------------#
# 3) Value of pooling ----
# Var(Ind) / Var(Pooled) (> 1 = pooling more
#    efficient; expect roughly 1 / (1 - pilot share))
# ---------------------------------------------------------------------#
fig3_df <- pilot_summary0 %>%
    filter(method %in% pilot_case2_methods) %>%
    select(method, pilot_frac, bias, residual_var, R2_target, coding_frac, emp_var) %>%
    pivot_wider(names_from = method, values_from = emp_var) %>%
    mutate(ratio_ind_pool = `Pilot-Allocated-Discard` / `Pilot-Allocated-Pooled`) %>%
    group_by(pilot_frac, R2_target, coding_frac) %>%
    summarise(ratio_ind_pool = mean(ratio_ind_pool), .groups = "drop") %>%
    mutate(pf_lab = pf_label(pilot_frac))

fig3 <- ggplot(fig3_df, aes(coding_frac, ratio_ind_pool, color = pf_lab)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ R2_target, labeller = rv_lab) +
    scale_color_viridis_d(name = "Pilot share", end = 0.9) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(accuracy = 1)) +
    labs(x = "Proportion coded", y = "Var(pilot discarded) / Var(pilot pooled)") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    guides(color = guide_legend(nrow = 2))
print(fig3)
ggsave("figures/pilot-sims/fig_pilot_pool_vs_ind.pdf", fig3, width = 9, height = 4.5)

# ---------------------------------------------------------------------#
# 4) Information spectrum ----
# ---------------------------------------------------------------------#
pilot_arm_labs <- sprintf("Optimal, pilot variances (%.0f%% of budget)", 100 * PILOT_FRACS)
arm_levels <- c("Optimal, reversed variances (anti-informative)",
                "Optimal, random variances (zero info)",
                "Proportional (no variance info used)", pilot_arm_labs)
arm_cols <- c("Optimal, reversed variances (anti-informative)" = "goldenrod",
              "Optimal, random variances (zero info)" = "firebrick",
              "Proportional (no variance info used)"  = "black",
              setNames(colorRampPalette(c("#6baed6", "#08519c"))(length(PILOT_FRACS)),
                       pilot_arm_labs))

info_df <- bind_rows(
    pilot_summary0 %>% filter(method == "Stratified") %>%
        mutate(arm = "Proportional (no variance info used)"),
    pilot_summary0 %>% filter(method == "Optimal-Random") %>%
        mutate(arm = "Optimal, random variances (zero info)"),
    pilot_summary0 %>% filter(method == "Optimal-Reversed") %>%
        mutate(arm = "Optimal, reversed variances (anti-informative)"),
    pilot_summary0 %>% filter(method == "Pilot-Allocated-Pooled") %>%
        mutate(arm = sprintf("Optimal, pilot variances (%.0f%% of budget)", 100 * pilot_frac))
) %>%
    left_join(opt_ref, by = c("bias", "residual_var", "R2_target", "n_code", "coding_frac")) %>%
    mutate(vr_vs_opt = emp_var / emp_var_opt) %>%
    group_by(arm, residual_var, R2_target, coding_frac) %>%
    summarise(vr_vs_opt = mean(vr_vs_opt), .groups = "drop") %>%
    mutate(arm = factor(arm, levels = arm_levels))

fig4 <- ggplot(info_df, aes(coding_frac, vr_vs_opt, color = arm, group = arm)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_grid(R2_target ~ residual_var, labeller = rv_lab) +
    scale_color_manual(name = NULL, values = arm_cols) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(accuracy = 1)) +
    labs(x = "Proportion coded", y = "Empirical variance relative to oracle Optimal") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 4, byrow = TRUE))
print(fig4)
ggsave("figures/pilot-sims/fig_pilot_information_spectrum.pdf", fig4, width = 10, height = 6.5)

info_tab <- pilot_summary0 %>%
    filter(method %in% c("Optimal", "Optimal-Random", "Optimal-Reversed", "Stratified",
                         "Pilot-Allocated-Pooled"),
           pf_near(coding_frac, HEADLINE_H)) %>%
    left_join(opt_ref, by = c("bias", "residual_var", "R2_target", "n_code", "coding_frac")) %>%
    mutate(vr_vs_opt = emp_var / emp_var_opt,
           Allocation = case_when(
               method == "Optimal"              ~ "Oracle variances",
               method == "Pilot-Allocated-Pooled" ~ sprintf("Pilot variances (%.0f%% of budget)", 100 * pilot_frac),
               method == "Optimal-Random"       ~ "Random variances",
               method == "Optimal-Reversed"     ~ "Reversed variances (anti-informative)",
               method == "Stratified"           ~ "Proportional (no variance info)")) %>%
    group_by(Allocation, residual_var) %>%
    summarise(`Emp. SE` = mean(emp_se),
              `Var. ratio vs oracle` = mean(vr_vs_opt),
              `Coverage` = mean(coverage), .groups = "drop") %>%
    arrange(residual_var, `Var. ratio vs oracle`)
print(info_tab, n = 40)

kbl(info_tab, format = "latex", booktabs = TRUE, digits = 3,
    caption = sprintf("The information spectrum for optimal allocation at a %.0f\\%% coding budget (balanced-exact design): oracle variance parameters vs. pilot estimates (pooled estimator, by pilot share of budget) vs. uninformative random draws vs. rank-reversed variances (the anti-informative floor), with proportional allocation as the no-variance-information reference. Coverage remains nominal throughout; only efficiency varies.", 100 * HEADLINE_H),
    label = "pilot_information_spectrum") %>%
    kable_styling(latex_options = "hold_position") %>%
    readr::write_file("tables/supplemental/pilot_information_spectrum.tex")

# ---------------------------------------------------------------------#
# 5) Validity table ----
# bias and coverage by method (ALL arms, incl. Case 1)
# ---------------------------------------------------------------------#
validity_tab <- pilot_summary0 %>%
    group_by(method) %>%
    summarise(`Mean bias` = mean(bias_est),
              `Max |bias|` = max(abs(bias_est)),
              `Mean coverage` = mean(coverage),
              `Min coverage` = min(coverage),
              `Max coverage` = max(coverage),
              `Mean SE` = mean(emp_se),
              `Mean SEhat` = mean(mean_sehat),
              `SE-hat / SE` = mean(mean_sehat / emp_se),
              .groups = "drop")
print(validity_tab)

kbl(validity_tab, format = "latex", booktabs = TRUE, digits = 4,
    caption = "Monte Carlo bias and 95\\% CI coverage by estimator, aggregated across all balanced-exact scenarios, coding budgets, and pilot budget fractions. All allocation regimes, including randomly drawn and pilot-estimated variance parameters and pilot-selected stratifications, remain unbiased with nominal coverage; only efficiency differs. The last column is the average ratio of the mean estimated standard error to the empirical standard error.",
    label = "pilot_sims_validity") %>%
    kable_styling(latex_options = "hold_position") %>%
    readr::write_file("tables/pilot_sims_validity.tex")

# ---------------------------------------------------------------------#
# 5b) Case 2 focused figures: pilot vs proportional (worth it?), the
#     pilot-share tradeoff, and gain capture (same outputs as before,
#     now computed from the same loaded results object)
# ---------------------------------------------------------------------#
pilot_r  <- ratios_all %>% filter(method == "Pilot-Allocated-Pooled")
oracle_r <- ratios_all %>% filter(method == "Optimal")

pf_cols <- setNames(viridisLite::viridis(length(PILOT_FRACS), end = 0.9), percent(PILOT_FRACS, 1))

# Aggregated (any bias x any heterogeneity) version, selected shares
ratios2 <- pilot_summary0 %>%
    filter(method %in% c("Stratified", "Optimal", "Pilot-Allocated-Pooled")) %>%
    group_by(method, pilot_frac, any_bias, any_residual_var, R2_target, coding_frac) %>%
    summarise(m = sum(n_iter), emp_var = mean(emp_var), .groups = "drop") %>%
    mutate(se_logvar = sqrt(2 / (m - 1)))
ratios2 <- ratios2 %>%
    filter(method != "Stratified") %>%
    left_join(ratios2 %>% filter(method == "Stratified") %>%
                  select(any_bias, any_residual_var, R2_target, coding_frac,
                         var_prop = emp_var, se_lv_prop = se_logvar),
              by = c("any_bias", "any_residual_var", "R2_target", "coding_frac")) %>%
    mutate(log_vr = log(emp_var / var_prop),
           se_log = sqrt(se_logvar^2 + se_lv_prop^2),
           pct_red = 100 * (1 - exp(log_vr)),
           pct_red_lo = 100 * (1 - exp(log_vr + 1.96 * se_log)),
           pct_red_hi = 100 * (1 - exp(log_vr - 1.96 * se_log)))

pilot_r2  <- ratios2 %>% filter(method == "Pilot-Allocated-Pooled",
                                pf_near(pilot_frac, 0.25) | pf_near(pilot_frac, 1/3) | pf_near(pilot_frac, 0.5))
oracle_r2 <- ratios2 %>% filter(method == "Optimal")

fig1b <- ggplot(pilot_r2, aes(coding_frac, pct_red, color = percent(pilot_frac, 1))) +
    geom_hline(yintercept = 0, linewidth = 0.45) +
    geom_line(linewidth = 0.95) +
    geom_point(size = 1.7) +
    geom_line(data = oracle_r2, aes(coding_frac, pct_red), inherit.aes = FALSE,
              linetype = "dashed", color = "grey25", linewidth = 0.8) +
    facet_grid(R2_target ~ any_residual_var + any_bias, labeller = rv_lab) +
    labs(col = "Pilot share of budget") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(1)) +
    scale_y_continuous(breaks = seq(-30, 30, by = 10), limits = c(-31, 30), oob = scales::oob_squish) +
    labs(x = "Total coding budget (pilot included)",
         y = "% variance reduction vs proportional stratified",
         caption = "Colored lines are pooled optimal estimator (pilot-estimated variances). Dashed grey: oracle Optimal (known variances), the ceiling on any pilot.\nCurves near 0 = piloting is incremental; curves near the dashed line = piloting captures most of what is achievable.") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(fig1b)
ggsave("figures/pilot-sims/fig_pilot_vs_prop_worthit.pdf", fig1b, width = 10, height = 6)

# Pilot-share tradeoff (Case 2): endpoints 0% and 100% are structural
hh_show <-  HH_SHOW[sapply(HH_SHOW, function(h) any(pf_near(pilot_r$coding_frac, h)))]

trade <- bind_rows(
    pilot_r %>% select(any_bias, any_residual_var, R2_target, coding_frac, pilot_frac,
                       pct_red, pct_red_lo, pct_red_hi),
    pilot_r %>% distinct(any_bias, any_residual_var, R2_target, coding_frac) %>%
        crossing(pilot_frac = c(0, 1)) %>%
        mutate(pct_red = 0, pct_red_lo = 0, pct_red_hi = 0)
) %>%
    filter(sapply(coding_frac, function(h) any(pf_near(hh_show, h)))) %>%
    group_by(any_bias, any_residual_var, R2_target, coding_frac, pilot_frac) %>%     # average the bias panels
    summarise(across(c(pct_red, pct_red_lo, pct_red_hi), mean), .groups = "drop") %>%
    mutate(h_lab = factor(percent(coding_frac, 1), levels = percent(sort(unique(coding_frac)), 1)),
           n_per_cell = ifelse(pilot_frac %in% c(0, 1), NA,
                               round(pilot_frac * coding_frac * N_SIM / N_CELLS)))

oracle_ref2 <- oracle_r %>%
    filter(sapply(coding_frac, function(h) any(pf_near(hh_show, h)))) %>%
    group_by(any_bias,any_residual_var, R2_target, coding_frac) %>%
    summarise(pct_red = mean(pct_red), .groups = "drop") %>%
    mutate(h_lab = factor(percent(coding_frac, 1), levels = levels(trade$h_lab)))

fig_trade <- ggplot(trade, aes(pilot_frac, pct_red, color = h_lab, group = h_lab)) +
    geom_hline(yintercept = 0, linewidth = 0.45) +
    geom_point()+
    #geom_hline(data = oracle_ref2, aes(yintercept = pct_red, color = h_lab),
    #          linetype = "dotted", show.legend = FALSE) +
    #geom_pointrange(aes(ymin = pct_red_lo, ymax = pct_red_hi),
    #                position = position_dodge(width = 0.02), size = 0.35) +
    geom_line(linewidth = 0.9) +
    #geom_text(aes(label = n_per_cell), nudge_y = 1.6, size = 2.7, show.legend = FALSE, na.rm = TRUE) +
    facet_grid(R2_target~ any_bias+any_residual_var , labeller = rv_lab) +
    scale_color_viridis_d(name = "Total budget", end = 0.85) +
    scale_x_continuous(breaks = seq(0,1,by=0.2), labels = percent_format(1)) +
    labs(x = "Share of the coding budget spent on the pilot",
         y = "% Variance Reduction vs Proportional Allocation\n (Pilot-allocated design, pilot pooled)") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(fig_trade)
ggsave("figures/fig_pilot_fraction_tradeoff.pdf", fig_trade, width = 9.5, height = 5)

# Gain capture bars (best share per budget; descriptive only)
cap <- pilot_r %>%
    group_by(bias, residual_var, R2_target, coding_frac) %>%
    slice_max(pct_red, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(bias, residual_var, R2_target, coding_frac, best_pf = pilot_frac, pilot_red = pct_red) %>%
    left_join(oracle_r %>% select(bias, residual_var, R2_target, coding_frac, oracle_red = pct_red),
              by = c("bias", "residual_var", "R2_target", "coding_frac")) %>%
    group_by(residual_var, coding_frac) %>%
    summarise(pilot_red = mean(pilot_red), oracle_red = mean(oracle_red),
              best_pf = best_pf[1], .groups = "drop") %>%
    mutate(captured = pmin(pmax(pilot_red, 0), pmax(oracle_red, 0)),
           unrealized = pmax(oracle_red - captured, 0)) %>%
    pivot_longer(c(captured, unrealized), names_to = "part", values_to = "gain") %>%
    mutate(part = factor(part, levels = c("unrealized", "captured")))

fig_cap <- ggplot(cap, aes(factor(percent(coding_frac, 1), levels = percent(sort(unique(coding_frac)), 1)),
                           gain, fill = part)) +
    geom_col(width = 0.7, color = "grey30", linewidth = 0.25) +
    geom_text(data = cap %>% filter(part == "captured"),
              aes(y = 0.4, label = percent(best_pf, 1)), size = 2.6, vjust = 0,
              color = "white", fontface = "bold") +
    facet_wrap(~ residual_var, labeller = rv_lab) +
    scale_fill_manual(name = NULL,
                      values = c("captured" = "#1f78b4", "unrealized" = "grey85"),
                      labels = c("captured" = "Captured by best pilot", "unrealized" = "Unrealized (oracle only)")) +
    labs(x = "Total coding budget", y = "% variance reduction vs proportional",
         caption = "Bar height = oracle Optimal's gain over proportional (the most any pilot could deliver).\nWhite label = the best-performing pilot share at that budget.") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
print(fig_cap)
ggsave("figures/pilot-sims/fig_pilot_gain_capture_bars.pdf", fig_cap, width = 9.5, height = 4.6)




# Aggregated (any bias x any heterogeneity) version, selected shares
ratios3 <- pilot_summary0 %>%
    filter(method %in% c("Stratified", "Optimal", pilot_methods)) %>%
    group_by(method, pilot_frac, any_bias, any_residual_var, R2_target, coding_frac) %>%
    summarise(m = sum(n_iter), emp_var = mean(emp_var), .groups = "drop") %>%
    mutate(se_logvar = sqrt(2 / (m - 1)))
ratios3 <- ratios3 %>%
    filter(method != "Stratified") %>%
    left_join(ratios3 %>% filter(method == "Stratified") %>%
                  select(any_bias, any_residual_var, R2_target, coding_frac,
                         var_prop = emp_var, se_lv_prop = se_logvar),
              by = c("any_bias", "any_residual_var", "R2_target", "coding_frac")) %>%
    mutate(log_vr = log(emp_var / var_prop),
           se_log = sqrt(se_logvar^2 + se_lv_prop^2),
           pct_red = 100 * (1 - exp(log_vr)),
           pct_red_lo = 100 * (1 - exp(log_vr + 1.96 * se_log)),
           pct_red_hi = 100 * (1 - exp(log_vr - 1.96 * se_log)))

pilot_r3  <- ratios3 %>% filter(pf_near(pilot_frac, 0.2) | is.na(pilot_frac))
oracle_r3 <- ratios3 %>% filter(method == "Optimal")

fig1c <- ggplot(pilot_r3, aes(coding_frac, pct_red, color = method,shape = method, linetype = method)) +
    geom_hline(yintercept = 0, linewidth = 0.45) +
    geom_line(linewidth = 0.95) +
    geom_point(size = 2) +
    facet_grid(R2_target ~ any_residual_var + any_bias, labeller = rv_lab) +
    labs(col = "Method") +
    labs(x = "Total coding budget (pilot included)",
         y = "% Variance Reduction vs Proportional Allocation")+    
    method_scales(use_labels = TRUE, nrow = 2)+
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())+
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), labels = percent_format(accuracy = 1),
                       limits = c(0.05, 1.05)) +
    scale_y_continuous(breaks = seq(-30, 30, by = 10), limits = c(-32, 32), oob = scales::squish,
                       labels = percent_format(scale = 1)) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(size = 11), legend.text = element_text(size = 9))
    
print(fig1c)
ggsave("figures/pilot-sims/fig_optimal_vs_prop_worthit.pdf", fig1c, width = 10, height = 6)

