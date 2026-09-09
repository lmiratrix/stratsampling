library(tidyverse)
library(janitor)
library(broom)
library(scales)
library(gt)

results_pilot <- readRDS("results/sim_results_pilot_balanced_exact_fast.rds")

PILOT_FRACS=sort(unique(results_pilot$pilot_frac))

#dir.create("figures/new-sims", recursive = TRUE, showWarnings = FALSE)

method_levels_pilot <- c("Stratified", "Optimal", "Optimal-Random",
                         "Optimal-Pilot-Pooled", "Optimal-Pilot-Ind")

pal6 <- c("Stratified"                = "black",
          "Optimal"                   = "forestgreen",
          "Optimal-Random"            = "firebrick",
          "Optimal-Pilot-Pooled" = "#1f77b4",
          "Optimal-Pilot-Ind"  = "#8c6bb1")

shapes6 <- c("Stratified" = 15, "Optimal" = 16, "Optimal-Random" = 17,
             "Optimal-Pilot-Pooled" = 18, "Optimal-Pilot-Ind" = 5)

lty6 <- c("Stratified" = "dashed", "Optimal" = "solid",
          "Optimal-Random" = "longdash",
          "Optimal-Pilot-Pooled" = "solid",
          "Optimal-Pilot-Ind" = "twodash")

pf_label <- function(pf) {
    factor(sprintf("Pilot = %.0f%% of budget", 100 * pf),
           levels = sprintf("Pilot = %.0f%% of budget", 100 * sort(unique(pf))))
}

# Full summary (all methods): source for fig2, fig3, fig4, validity
pilot_summary0 <- results_pilot %>%
    filter(!is.na(est), !is.na(SEhat)) %>%
    mutate(coding_frac = n_code / N,
           CI_lower = est - 2 * SEhat,
           CI_upper = est + 2 * SEhat,
           covers = if_else(CI_lower < tau & CI_upper > tau, 1L, 0L)) %>%
    group_by(method, pilot_frac, bias, residual_var, R2_target,
             n_code, coding_frac) %>%
    summarise(
        mean_est   = mean(est),
        bias_est   = mean(est - tau),
        coverage   = mean(covers),
        emp_var    = var(est),
        emp_se     = sd(est),
        mean_sehat = mean(SEhat),
        mse        = mean((est - tau)^2),
        n_iter     = n(),
        .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = method_levels_pilot),
           bias = fct_relevel(factor(bias), "none", "large"),
           residual_var = fct_relevel(factor(residual_var), "homogeneous",
                                      "heterogeneous"))

# Headline summary for the frontier figure only (drop the Ind arm)
pilot_summary <- pilot_summary0 %>% filter(method != "Optimal-Pilot-Ind")

# ---------------------------------------------------------------------
# Figure 1: efficiency frontiers (aggregated over scenarios), by pilot
# fraction; baselines (incl. Optimal-Random) replicated across facets
# ---------------------------------------------------------------------

plot_df <- bind_rows(
    pilot_summary %>%
        filter(is.na(pilot_frac)) %>%
        select(-pilot_frac) %>%
        crossing(pilot_frac = PILOT_FRACS),
    pilot_summary %>% filter(!is.na(pilot_frac))
) %>%
    mutate(pf_lab = pf_label(pilot_frac),
           any_bias=ifelse(bias=="none",0,1), 
           any_residual_var=ifelse(residual_var=="homogeneous",0,1))

fig1_df <- plot_df %>%
    group_by(method, pilot_frac, R2_target, coding_frac, any_bias, any_residual_var) %>%
    summarise(emp_se = mean(emp_se), .groups = "drop")

fig1_df = fig1_df %>% 
    mutate(any_bias = factor(any_bias, levels=c(0,1), labels=c("No Bias", "Biased")),
           any_residual_var=factor(any_residual_var, levels=c(0,1), labels=c("Homogenous Var.", "Heterogeneous Var.")))

fig1 <- ggplot(fig1_df %>% filter(is.na(pilot_frac) | pilot_frac==0.5), aes(coding_frac, emp_se, color = method,
                            shape = method, linetype = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_grid(R2_target ~ any_bias+any_residual_var) +
    scale_color_manual(name = "Method", values = pal6) +
    scale_shape_manual(name = "Method", values = shapes6) +
    scale_linetype_manual(name = "Method", values = lty6) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(accuracy = 1)) +
    #scale_y_continuous(breaks=seq(0.2,0.5,by=0.05),limits=c(0.1,0.6))+
    #labs(x = "Proportion coded (total budget, pilot included)",
      #   y = "Empirical Std. Error (avg. across scenarios)") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(fig1)
ggsave("figures/new-sims/fig_pilot_frontier.pdf", fig1, width = 10, height = 6.5)

# ---------------------------------------------------------------------
# Figure 2: variance of each pilot arm and of Optimal-Random relative
# to oracle Optimal at the same budget (Optimal-Random shown as pilot
# fraction 0 = zero information)
# ---------------------------------------------------------------------

opt_ref <- pilot_summary %>%
    filter(method == "Optimal") %>%
    select(bias, residual_var, R2_target, n_code, coding_frac,
           emp_var_opt = emp_var)

fig2_df <- pilot_summary %>%
    filter(method %in% c("Optimal-Pilot-Pooled","Stratified"
                        # , "Optimal-Pilot-Ind", "Optimal-Random"
                        )) %>%
    left_join(opt_ref, by = c("bias", "residual_var", "R2_target",
                              "n_code", "coding_frac")) %>%
    mutate(vr_vs_opt = emp_var / emp_var_opt) %>%
    mutate(pf_lab = pf_label(pilot_frac),
           any_bias=ifelse(bias=="none",0,1), 
           any_residual_var=ifelse(residual_var=="homogeneous",0,1))%>%
      mutate(any_bias = factor(any_bias, levels=c(0,1), labels=c("No Bias", "Biased")),
               any_residual_var=factor(any_residual_var, levels=c(0,1), labels=c("Homogenous Var.", "Heterogeneous Var."))) %>%
          group_by(method, pilot_frac, any_bias, any_residual_var, R2_target, coding_frac) %>%
    summarise(vr_vs_opt = mean(vr_vs_opt), .groups = "drop") %>%
    mutate(pf_lab = as.factor(ifelse(is.na(pilot_frac), 0, pilot_frac)))

fig2 <- ggplot(fig2_df %>% filter(pilot_frac%in%c(1/4,1/3,1/2) | is.na(pilot_frac)), aes(coding_frac, vr_vs_opt,
                            color = pf_lab, linetype = method,
                            group = interaction(pf_lab, method))) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_grid(R2_target ~ any_bias+any_residual_var)+
    scale_linetype_manual(name = "Method",
                          values = c("Optimal-Pilot-Pooled" = "solid",
                                     "Stratified" = "dashed")) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Proportion coded",
         y = "Empirical variance relative to oracle Optimal",
         col = "Pilot fraction (0 = stratified)") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(fig2)
ggsave("figures/new-sims/fig_is_pilot_worthwhile.pdf", fig2, width = 10, height = 6.5)

# ---------------------------------------------------------------------
# Figure 3: value of pooling (Ind variance / Pool variance; > 1 means
# pooling is more efficient; expect roughly 1/(1 - pilot share))
# ---------------------------------------------------------------------

fig3_df <- pilot_summary0 %>%
    filter(method %in% c("Optimal-Pilot-Pooled", "Optimal-Pilot-Ind")) %>%
    select(method, pilot_frac, bias, residual_var, R2_target,
           coding_frac, emp_var) %>%
    pivot_wider(names_from = method, values_from = emp_var) %>%
    mutate(ratio_ind_pool = `Optimal-Pilot-Ind` / `Optimal-Pilot-Pooled`) %>%
    group_by(pilot_frac, R2_target, coding_frac) %>%
    summarise(ratio_ind_pool = mean(ratio_ind_pool), .groups = "drop") %>%
    mutate(pf_lab = pf_label(pilot_frac))

fig3 <- ggplot(fig3_df, aes(coding_frac, ratio_ind_pool, color = pf_lab)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~ R2_target,
               labeller = labeller(R2_target = c("0.4" = "R\u00b2 target = 0.4",
                                                 "0.85" = "R\u00b2 target = 0.85"))) +
    scale_color_brewer(name = "Pilot share", palette = "Dark2") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Proportion coded",
         y = "Var(pilot discarded) / Var(pilot pooled)") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(fig3)
ggsave("figures/new-sims/fig_pilot_pool_vs_ind.pdf", fig3, width = 9, height = 4.5)

# ---------------------------------------------------------------------
# Figure 4 + table: the information spectrum. Oracle variances (ratio
# = 1 by construction) vs pilot-estimated variances (Pool, by pilot
# share) vs random variances (zero information), with proportional
# allocation as the ignore-the-variances reference.
# ---------------------------------------------------------------------

pilot_arm_labs <- sprintf("Optimal, pilot variances (%.0f%% of budget)",
                          100 * sort(PILOT_FRACS))
arm_levels <- c("Optimal, random variances (zero info)",
                "Proportional (no variance info used)",
                pilot_arm_labs)
arm_cols <- c("Optimal, random variances (zero info)" = "firebrick",
              "Proportional (no variance info used)"  = "black",
              setNames(colorRampPalette(c("#6baed6", "#08519c"))(length(PILOT_FRACS)),
                       pilot_arm_labs))

info_df <- bind_rows(
    pilot_summary0 %>%
        filter(method == "Stratified") %>%
        mutate(arm = "Proportional (no variance info used)"),
    pilot_summary0 %>%
        filter(method == "Optimal-Random") %>%
        mutate(arm = "Optimal, random variances (zero info)"),
    pilot_summary0 %>%
        filter(method == "Optimal-Pilot-Pooled") %>%
        mutate(arm = sprintf("Optimal, pilot variances (%.0f%% of budget)",
                             100 * pilot_frac))
) %>%
    left_join(opt_ref, by = c("bias", "residual_var", "R2_target",
                              "n_code", "coding_frac")) %>%
    mutate(vr_vs_opt = emp_var / emp_var_opt) %>%
    group_by(arm, residual_var, R2_target, coding_frac) %>%
    summarise(vr_vs_opt = mean(vr_vs_opt), .groups = "drop") %>%
    mutate(arm = factor(arm, levels = arm_levels))

fig4 <- ggplot(info_df, aes(coding_frac, vr_vs_opt, color = arm, group = arm)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_grid(R2_target ~ residual_var,
               labeller = labeller(
                   R2_target = c("0.4" = "R\u00b2 target = 0.4",
                                 "0.85" = "R\u00b2 target = 0.85"),
                   residual_var = c("homogeneous" = "Residual var: homogeneous",
                                    "heterogeneous" = "Residual var: heterogeneous"))) +
    scale_color_manual(name = NULL, values = arm_cols) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(accuracy = 1)) +
    ggplot2::labs(x = "Proportion coded",
         y = "Empirical variance relative to oracle Optimal") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    guides(color = guide_legend(nrow = 2))
print(fig4)
ggsave("figures/new-sims/fig_pilot_information_spectrum.pdf", fig4,
       width = 10, height = 6.5)

# Companion table at the 30% budget
info_tab <- bind_rows(
    pilot_summary0 %>%
        filter(method %in% c("Optimal", "Optimal-Random", "Stratified")),
    pilot_summary0 %>%
        filter(method == "Optimal-Pilot-Pooled")
) %>%
    filter(coding_frac == 0.3) %>%
    left_join(opt_ref, by = c("bias", "residual_var", "R2_target",
                              "n_code", "coding_frac")) %>%
    mutate(vr_vs_opt = emp_var / emp_var_opt,
           Allocation = case_when(
               method == "Optimal"                   ~ "Oracle variances",
               method == "Optimal-Pilot-Pooled" ~ sprintf("Pilot variances (%.0f%% of budget)",
                                                               100 * pilot_frac),
               method == "Optimal-Random"            ~ "Random variances",
               method == "Stratified"                ~ "Proportional (no variance info)")) %>%
    group_by(Allocation, residual_var) %>%
    summarise(`Emp. SE` = mean(emp_se),
              `Var. ratio vs oracle` = mean(vr_vs_opt),
              `Coverage` = mean(coverage),
              .groups = "drop") %>%
    arrange(residual_var, `Var. ratio vs oracle`)

print(info_tab, n = 20)

kbl(info_tab, format = "latex", booktabs = TRUE, digits = 3,
    caption = "The information spectrum for optimal allocation at a 30\\% coding budget (balanced-exact design): oracle variance parameters vs. pilot estimates (pooled estimator, by pilot share of budget) vs. uninformative random draws, with proportional allocation as the no-variance-information reference. Coverage remains nominal throughout; only efficiency varies.",
    label = "pilot_information_spectrum") %>%
    kable_styling(latex_options = "hold_position") %>%
    readr::write_file("results/pilot_information_spectrum.tex")

# ---------------------------------------------------------------------
# Validity table: bias and coverage by method (ALL methods, incl. Ind)
# ---------------------------------------------------------------------

validity_tab <- pilot_summary0 %>%
    group_by(method) %>%
    summarise(`Mean bias` = mean(bias_est),
              `Max |bias|` = max(abs(bias_est)),
              `Mean coverage` = mean(coverage),
              `Min coverage` = min(coverage),
              `Max coverage` = max(coverage),
              .groups = "drop")

print(validity_tab)

kbl(validity_tab, format = "latex", booktabs = TRUE, digits = 4,
    caption = "Monte Carlo bias and 95\\% CI coverage by estimator, aggregated across all balanced-exact scenarios, coding budgets, and pilot budget fractions. All allocation regimes, including randomly drawn and pilot-estimated variance parameters, remain unbiased with nominal coverage; only efficiency differs.",
    label = "pilot_sims_validity") %>%
    kable_styling(latex_options = "hold_position") %>%
    readr::write_file("results/pilot_sims_validity.tex")




# =====================================================================
# 03d_pilot_figures.R
#
# Focused figures for the pilot-based optimal allocation study:
#
#  Fig 1: Is piloting worthwhile, or incremental, relative to the
#         simpler proportional stratified design? (Oracle Optimal
#         shown as the reference ceiling.)
#  Fig 2: How does the answer vary with the PILOT FRACTION? The
#         tradeoff curve: too small a pilot -> noisy variance
#         estimates; too large -> nothing left to reallocate; both
#         endpoints (0% and 100% of budget) are structurally identical
#         to proportional, so any value lives in between.
#  Fig 3: Of the oracle's achievable gain, how much does the best
#         pilot capture at each budget? (Meaningful vs incremental,
#         in one bar per budget.)
# =====================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

dir.create("figures/new-sims", recursive = TRUE, showWarnings = FALSE)

results_pilot <- readRDS("results/new-sims/sim_results_pilot_balanced_exact.rds") 

PILOT_FRACS <- sort(unique(na.omit(results_pilot$pilot_frac)))
N_SIM   <- unique(results_pilot$N)[1]
N_CELLS <- unique(results_pilot$n_strat)[1] * 2   # strata x arms

rv_lab <- labeller(
  residual_var = c("homogeneous"   = "Homogeneous residual var.\n(nothing to exploit)",
                   "heterogeneous" = "Heterogeneous residual var.",
                   "extreme-contrast" = "Extreme-contrast residual var."),
  bias = c("none" = "Bias: none", "large" = "Bias: large", 
           "moderate" = "Bias: moderate", "small"="Bias: small",
           "extreme-contrast" = "Bias: extreme-contrast"))
# ---------------------------------------------------------------------
# Per-cell variance summaries and ratios vs the proportional baseline
# (MC error on the log scale: Var(log s^2) ~ 2/(m-1))
# ---------------------------------------------------------------------

cells <- results_pilot %>%
  filter(method %in% c("Stratified", "Optimal", "Optimal-Pilot-Pooled"),
         !is.na(est)) %>%
  mutate(h = n_code / N) %>%
  group_by(method, pilot_frac, bias, residual_var, R2_target, h) %>%
  summarise(m = n(), emp_var = var(est), .groups = "drop") %>%
  mutate(se_logvar = sqrt(2 / (m - 1)))

ratios <- cells %>%
  filter(method != "Stratified") %>%
  left_join(cells %>% filter(method == "Stratified") %>%
              select(bias, residual_var, R2_target, h,
                     var_prop = emp_var, se_lv_prop = se_logvar),
            by = c("bias", "residual_var", "R2_target", "h")) %>%
  mutate(log_vr     = log(emp_var / var_prop),
         se_log     = sqrt(se_logvar^2 + se_lv_prop^2),
         pct_red    = 100 * (1 - exp(log_vr)),
         pct_red_lo = 100 * (1 - exp(log_vr + 1.96 * se_log)),
         pct_red_hi = 100 * (1 - exp(log_vr - 1.96 * se_log)))

pilot_r  <- ratios %>% filter(method == "Optimal-Pilot-Pooled")
oracle_r <- ratios %>% filter(method == "Optimal")

pf_cols <- setNames(RColorBrewer::brewer.pal(max(3, length(PILOT_FRACS)),
                                             "Dark2")[seq_along(PILOT_FRACS)],
                    percent(PILOT_FRACS, 1))

# =====================================================================
# Figure 1: worthwhile or incremental? Pilot vs proportional (y = 0
# line), with oracle Optimal as the reference ceiling.
# =====================================================================

fig1 <- ggplot(pilot_r%>%filter(R2_target==0.4),
               aes(h, pct_red, color = percent(pilot_frac, 1),
                   fill = percent(pilot_frac, 1))) +
  geom_hline(yintercept = 0, linewidth = 0.45) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.7) +
  geom_line(data = oracle_r%>%filter(R2_target==0.4), aes(h, pct_red),
            inherit.aes = FALSE, linetype = "dashed",
            color = "grey25", linewidth = 0.8) +
  facet_grid(bias ~ residual_var, labeller = rv_lab) +
  scale_color_manual(name = "Pilot share of budget", values = pf_cols) +
  scale_fill_manual(name = "Pilot share of budget", values = pf_cols) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(1)) +
  labs(x = "Total coding budget (pilot included)",
       y = "% variance reduction vs proportional stratified",
       caption = "Dashed grey: oracle Optimal (known variances), the ceiling on any pilot. Ribbons: 95% MC intervals.\nRead: curves near 0 = piloting is incremental; curves near the dashed line = piloting captures most of what is achievable.") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(fig1)
ggsave("figures/new-sims/fig_pilot_vs_prop_worthit.pdf", fig1,
       width = 10, height = 6)


# =====================================================================
# Figure 1b: worthwhile or incremental? Aggregated by bias and variance structure
# =====================================================================

cells2 <- results_pilot %>%
  filter(method %in% c("Stratified", "Optimal", "Optimal-Pilot-Pooled"),
         !is.na(est)) %>%
  mutate(h = n_code / N,
         any_bias=ifelse(bias=="none",0,1),
         any_var=ifelse(residual_var=="homogeneous",0,1)) %>%
  mutate(any_bias=factor(any_bias,levels=0:1,labels=c("No Bias","Bias")),
         any_var=factor(any_var,levels=0:1,labels=c("Homogeneous Var.","Heterogeneous Var."))) %>%
  group_by(method, pilot_frac, any_bias, any_var, R2_target, h) %>%
  summarise(m = n(), emp_var = var(est), .groups = "drop") %>%
  mutate(se_logvar = sqrt(2 / (m - 1)))

ratios2 <- cells2 %>%
  filter(method != "Stratified") %>%
  left_join(cells2 %>% filter(method == "Stratified") %>%
              select(any_bias, any_var, R2_target, h,
                     var_prop = emp_var, se_lv_prop = se_logvar),
            by = c("any_bias", "any_var", "R2_target", "h")) %>%
  mutate(log_vr     = log(emp_var / var_prop),
         se_log     = sqrt(se_logvar^2 + se_lv_prop^2),
         pct_red    = 100 * (1 - exp(log_vr)),
         pct_red_lo = 100 * (1 - exp(log_vr + 1.96 * se_log)),
         pct_red_hi = 100 * (1 - exp(log_vr - 1.96 * se_log)))

pilot_r2  <- ratios2 %>% filter(method == "Optimal-Pilot-Pooled",pilot_frac%in%c(0.25,1/3,0.5))
oracle_r2 <- ratios2 %>% filter(method == "Optimal", pilot_frac%in%c(0.25,1/3,0.5,NA))


fig1b <- ggplot(pilot_r2,
                aes(h, pct_red, color = percent(pilot_frac, 1),
                    fill = percent(pilot_frac, 1))) +
  geom_hline(yintercept = 0, linewidth = 0.45) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.7) +
  geom_line(data = oracle_r2, aes(h, pct_red),
            inherit.aes = FALSE, linetype = "dashed",
            color = "grey25", linewidth = 0.8) +
  facet_grid(R2_target ~ any_var+any_bias, labeller = rv_lab) +
  scale_color_manual(name = "Pilot share of budget", values = pf_cols) +
  scale_fill_manual(name = "Pilot share of budget", values = pf_cols) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), labels = percent_format(1)) +
  scale_y_continuous(breaks=seq(-30,30,by=10),limits=c(-31,30))+
  labs(x = "Total coding budget (pilot included)",
       y = "% variance reduction vs proportional stratified",
       caption = "Dashed grey: oracle Optimal (known variances), the ceiling on any pilot. Ribbons: 95% MC intervals.\nRead: curves near 0 = piloting is incremental; curves near the dashed line = piloting captures most of what is achievable.") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(fig1b)
ggsave("figures/new-sims/fig_pilot_vs_prop_worthit.pdf", fig1b,
       width = 10, height = 6)


# =====================================================================
# Figure 2: the pilot-fraction tradeoff. Both endpoints are EXACTLY
# proportional: at f = 0 there is no variance information, at f = 1
# there is no budget left to reallocate. Point labels show the implied
# pilot size per stratum-by-arm cell, the mechanism behind the
# too-small regime.
# =====================================================================

hh_show <- intersect(c(0.1, 0.3, 0.5, 0.7), unique(pilot_r$h))

trade <- bind_rows(
  pilot_r %>%
    select(bias, residual_var, R2_target, h, pilot_frac,
           pct_red, pct_red_lo, pct_red_hi),
  pilot_r %>%
    distinct(bias, residual_var, R2_target, h) %>%
    crossing(pilot_frac = c(0, 1)) %>%
    mutate(pct_red = 0, pct_red_lo = 0, pct_red_hi = 0)
) %>%
  filter(h %in% hh_show) %>%
  group_by(residual_var, R2_target, h, pilot_frac) %>%   # average the bias panels
  summarise(across(c(pct_red, pct_red_lo, pct_red_hi), mean),
            .groups = "drop") %>%
  mutate(h_lab = factor(percent(h, 1), levels = percent(sort(unique(h)), 1)),
         n_per_cell = ifelse(pilot_frac %in% c(0, 1), NA,
                             round(pilot_frac * h * N_SIM / N_CELLS)))

oracle_ref2 <- oracle_r %>%
  filter(h %in% hh_show) %>%
  group_by(residual_var, R2_target, h) %>%
  summarise(pct_red = mean(pct_red), .groups = "drop") %>%
  mutate(h_lab = factor(percent(h, 1), levels = levels(trade$h_lab)))

fig2 <- ggplot(trade, aes(pilot_frac, pct_red, color = h_lab, group = h_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.45) +
  geom_hline(data = oracle_ref2, aes(yintercept = pct_red, color = h_lab),
             linetype = "dotted", show.legend = FALSE) +
  geom_pointrange(aes(ymin = pct_red_lo, ymax = pct_red_hi),
                  position = position_dodge(width = 0.02), size = 0.35) +
  geom_line(linewidth = 0.9) +
  geom_text(aes(label = n_per_cell), nudge_y = 1.6, size = 2.7,
            show.legend = FALSE, na.rm = TRUE) +
  facet_grid(R2_target~ residual_var, labeller = rv_lab) +
  scale_color_viridis_d(name = "Total budget", end = 0.85) +
  scale_x_continuous(breaks = sort(unique(c(0, PILOT_FRACS, 1))),
                     labels = percent_format(1)) +
  labs(x = "Share of the coding budget spent on the pilot",
       y = "% variance reduction vs proportional",
       caption = "Endpoints are structural: 0% and 100% pilot both reproduce the proportional design exactly.\nDotted lines: oracle ceiling per budget. Point labels: implied pilot documents per stratum-by-arm cell.") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
print(fig2)
ggsave("figures/new-sims/fig_pilot_fraction_tradeoff.pdf", fig2,
       width = 9.5, height = 5)

# =====================================================================
# Figure 3: meaningful vs incremental, one bar per budget. Bar height
# = oracle's achievable gain over proportional; filled part = share
# captured by the best pilot fraction; hollow top = unrealized.
# =====================================================================

cap <- pilot_r %>%
  group_by(bias, residual_var, R2_target, h) %>%
  slice_max(pct_red, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(bias, residual_var, R2_target, h,
         best_pf = pilot_frac, pilot_red = pct_red) %>%
  left_join(oracle_r %>%
              select(bias, residual_var, R2_target, h,
                     oracle_red = pct_red),
            by = c("bias", "residual_var", "R2_target", "h")) %>%
  group_by(residual_var, h) %>%                          # average bias panels
  summarise(pilot_red  = mean(pilot_red),
            oracle_red = mean(oracle_red),
            best_pf    = best_pf[1],
            .groups = "drop") %>%
  mutate(captured   = pmin(pmax(pilot_red, 0), pmax(oracle_red, 0)),
         unrealized = pmax(oracle_red - captured, 0)) %>%
  pivot_longer(c(captured, unrealized),
               names_to = "part", values_to = "gain") %>%
  mutate(part = factor(part, levels = c("unrealized", "captured")))

fig3 <- ggplot(cap, aes(factor(percent(h, 1),
                               levels = percent(sort(unique(h)), 1)),
                        gain, fill = part)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.25) +
  geom_text(data = cap %>% filter(part == "captured"),
            aes(y = 0.4, label = percent(best_pf, 1)),
            size = 2.6, vjust = 0, color = "white", fontface = "bold") +
  facet_wrap(~ residual_var, labeller = rv_lab) +
  scale_fill_manual(name = NULL,
                    values = c("captured" = "#1f78b4",
                               "unrealized" = "grey85"),
                    labels = c("captured" = "Captured by best pilot",
                               "unrealized" = "Unrealized (oracle only)")) +
  labs(x = "Total coding budget",
       y = "% variance reduction vs proportional",
       caption = "Bar height = oracle Optimal's gain over proportional (the most any pilot could deliver).\nWhite label = the best-performing pilot share at that budget.") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
print(fig3)
ggsave("figures/new-sims/fig_pilot_gain_capture_bars.pdf", fig3,
       width = 9.5, height = 4.6)