# ===================================================================== #
# 03d_pilot_selection_eval.R
# How often does the pilot-based diagnostic choose the TRUE strata, and
# under what scenarios does it work / fail?
# Input: results/sim_results_pilot_balanced_exact_fast_v2.rds (02_run_sim_pilot_fast.R)
# ===================================================================== #

library(tidyverse)
library(scales)
library(kableExtra)

RESULTS_FILE <- "results/sim_results_pilot_balanced_exact_fast_v2.rds"
HEADLINE_H   <- 0.2      # coding budget for Fig 1 and the table
HEADLINE_PF  <- 1/5      # pilot share for the table
near <- function(x, t, tol = 1e-6) abs(x - t) < tol

family_levels <- c("True strata", "Noisy strata", "Yhat cuts", "Random covariate")
family_pal    <- c("True strata" = "#1b9e77", "Noisy strata" = "#7fc97f",
                   "Yhat cuts" = "#7570b3", "Random covariate" = "#d95f02")
cand_family <- function(x) {
    factor(case_when(x == "strata_true"   ~ "True strata",
                     x == "strata_noisy"  ~ "Noisy strata",
                     grepl("^Yhat_", x)   ~ "Yhat cuts",
                     x == "strata_random" ~ "Random covariate"), levels = family_levels)
}
# Whether the true strata actually matter for estimation in this DGP
signal_class <- function(bias, residual_var) {
    factor(case_when(bias == "none" & residual_var == "homogeneous" ~ "No structure",
                     bias != "none" & residual_var == "homogeneous" ~ "Bias only",
                     bias == "none" & residual_var != "homogeneous" ~ "Variance only",
                     TRUE ~ "Bias + variance"),
           levels = c("No structure", "Bias only", "Variance only", "Bias + variance"))
}
lab_rv   <- c("homogeneous" = "Homogeneous", "heterogeneous" = "Heterogeneous", "extreme-contrast" = "Extreme-contrast")
lab_bias <- c("none" = "none", "small" = "small", "moderate" = "moderate", "large" = "large", "extreme-contrast" = "extreme")
lab_R2   <- c("0.4" = "R\u00b2 = 0.4", "0.85" = "R\u00b2 = 0.85")

# ---------------------------------------------------------------------#
# One row per iteration x scenario x pilot share: what was chosen
# ---------------------------------------------------------------------#
results_pilot <- readRDS(RESULTS_FILE)

pilot_small=readRDS("results/sim_results_pilot_small.rds")
results_pilot=rbind(results_pilot, pilot_small)


sel <- results_pilot %>%
    filter(method == "Pilot-Selected-Prop") %>%
    transmute(.n, bias, residual_var, R2_target, n_code, coding_frac = n_code / N,
              pilot_frac, n_pilot, family = cand_family(strat_choice),
              signal = signal_class(bias, residual_var))
stopifnot(!anyNA(sel$family))

# ---------------------------------------------------------------------#
# Fig 1: P(true strata chosen) over the scenario grid, headline budget ----
# Rows = pilot share, columns = bias level, panels = residual var x R2.
# Cell label gives the Yhat-cut share so the two "competitors" are visible.
# ---------------------------------------------------------------------#
grid <- sel %>%
    filter(near(coding_frac, HEADLINE_H)) %>%
    group_by(bias, residual_var, R2_target, signal, pilot_frac, n_pilot) %>%
    summarise(p_true = mean(family == "True strata"),
              p_yhat = mean(family == "Yhat cuts"),
              p_noisy = mean(family == "Noisy strata"),
              p_random = mean(family == "Random covariate"),
              n_iter = n(), .groups = "drop")

fig1 <- grid %>%
    mutate(bias = factor(lab_bias[bias], levels = lab_bias),
           residual_var = factor(lab_rv[residual_var], levels = lab_rv),
           pf = factor(percent(pilot_frac, 1), levels = percent(sort(unique(pilot_frac)), 1))) %>%
    ggplot(aes(x = bias, y = pf, fill = p_true)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.0f\n(%.0f)", 100 * p_true, 100 * p_yhat)), size = 2.4) +
    facet_grid(R2_target ~ residual_var, labeller = labeller(R2_target = lab_R2)) +
    scale_fill_gradient(name = "P(true strata chosen)", low = "red", high = "#1b9e77",
                        limits = c(0, 1), labels = percent) +
    labs(x = "Bias of surrogate across true strata",
         y = "Pilot share of the coding budget",
         title = sprintf("Pilot-based selection at a %.0f%% coding budget", 100 * HEADLINE_H),
         subtitle = "Cell: % of iterations choosing the true strata (% choosing a Yhat cut in parentheses); columns = residual-variance pattern") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid = element_blank())
fig1
ggsave("figures/pilot-sims/fig_select_true_heatmap.pdf", fig1, width = 10, height = 8)

# ---------------------------------------------------------------------#
# Fig 2: hit rate vs pilot size, by residual structure ----
# Pools all budgets and shares; x = realised pilot size.  "Works" = the
# true-strata curve climbs toward 1; "fails" = it stays flat and the
# Yhat / random curves do not fall.
# ---------------------------------------------------------------------#
curve <- sel %>%
    group_by(signal, R2_target, n_pilot) %>%
    summarise(`True strata` = mean(family == "True strata"),
              `Noisy strata` = mean(family == "Noisy strata"),
              `Yhat cuts` = mean(family == "Yhat cuts"),
              `Random covariate` = mean(family == "Random covariate"),
              n_iter = n(), .groups = "drop") %>%
    pivot_longer(all_of(family_levels), names_to = "family", values_to = "p") %>%
    mutate(family = factor(family, levels = family_levels))

fig2 <- curve %>%
    ggplot(aes(x = n_pilot, y = p, colour = family)) +
    geom_line() + geom_point(size = 0.9) +
    facet_grid(R2_target ~ signal, labeller = labeller(R2_target = lab_R2)) +
    #scale_x_log10() +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    scale_colour_manual(name = "Chosen stratification", values = family_pal) +
    labs(x = "Pilot size (coded units)", y = "Share of iterations",
         title = "When the pilot diagnostic works: selection rate by pilot size and residual structure") +
    theme_minimal() + theme(legend.position = "bottom")
fig2
ggsave("figures/pilot-sims/fig_select_learning_curve.pdf", fig2, width = 10, height =6)

# ---------------------------------------------------------------------#
# Table: headline design (budget HEADLINE_H, pilot share HEADLINE_PF) ----
# ---------------------------------------------------------------------#
tab <- grid %>%
    filter(near(pilot_frac, HEADLINE_PF)) %>%
    arrange(R2_target, signal, bias, residual_var) %>%
    transmute(R2 = R2_target, Structure = signal, Bias = bias, `Resid. var.` = residual_var,
              `True strata` = percent(p_true, 1), `Noisy strata` = percent(p_noisy, 1),
              `Yhat cuts` = percent(p_yhat, 1), `Random` = percent(p_random, 1))
print(tab, n = Inf)

tab %>%
    kbl(format = "latex", booktabs = TRUE, linesep = "", escape = FALSE,
        col.names = c("$R^2$", "Structure", "Bias", "Resid. var.",
                      "True strata", "Noisy strata", "Yhat cuts", "Random"),
        caption = sprintf(paste(
            "Share of iterations in which each candidate family wins the pre-coding diagnostic,",
            "at a %.0f\\%% coding budget with %.0f\\%% of it spent on the pilot (%d pilot units).",
            "`Structure' classifies whether the true strata carry any signal in the surrogate",
            "residuals (bias across strata, heterogeneous residual variance, both, or neither)."),
            100 * HEADLINE_H, 100 * HEADLINE_PF, unique(grid$n_pilot[near(grid$pilot_frac, HEADLINE_PF)])),
        label = "tab:pilot-selection") %>%
    add_header_above(c(" " = 4, "Share chosen" = 4)) %>%
    collapse_rows(columns = 1:2, latex_hline = "major") %>%
    writeLines("tables/tab_pilot_selection.tex")

