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



# Reproduce your summary object
analysis_results <- results %>%
    filter(!is.na(est)) %>%  # (you noted some NAs under unbalanced + small n_code)
    group_by(N, n_strat, strata_config, tau, bias, residual_var, R2_target, pop.sd, n_code, method) %>%
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
        bias          = fct_relevel(factor(bias), "none", "moderate", "large", "extreme-contrast"),
        strata_config = fct_relevel(factor(strata_config), "balanced-exact", "balanced-approx", "unbalanced")
    )

# Keys that index a "cell" apart from method
cell_keys <- c("N","n_strat","strata_config","tau","bias","residual_var","R2_target","pop.sd","n_code","coding_frac")

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
    group_by(N, n_strat, strata_config, tau, bias, residual_var, R2_target, pop.sd, method) %>%
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
    filter(N == 1000, n_strat == 4, strata_config == "balanced-exact", pop.sd == 3) %>%
    mutate(
        residual_var  = fct_relevel(residual_var, "homogeneous", "heterogeneous", "extreme-contrast"),
        bias          = fct_relevel(bias, "none", "moderate", "large", "extreme-contrast")
    )

# A consistent color palette for methods
method_pal <- c(
    "Subset"      = "#9467bd",
    "SRS"         = "#1f77b4",
    "Stratified"  = "#2ca02c",
    "Optimal"     = "#ff6961",
    "Full coding" = "black"
)

# Helper for clean labels
lab_resid <- labeller(residual_var = c(
    "homogeneous" = "Residual variance: homogeneous",
    "heterogeneous" = "Residual variance: heterogeneous",
    "extreme-contrast" = "Residual variance: extreme-contrast"
))
lab_bias <- labeller(bias = c(
    "none" = "Bias: none",
    "moderate" = "Bias: moderate",
    "large" = "Bias: large",
    "extreme-contrast" = "Bias: extreme-contrast"
))

labs = labeller(residual_var = c(
    "homogeneous" = "Residual variance: homogeneous",
    "heterogeneous" = "Residual variance: heterogeneous",
    "extreme-contrast" = "Residual variance: extreme-contrast"
),
bias = c(
    "none" = "Bias: none",
    "moderate" = "Bias: moderate",
    "large" = "Bias: large",
    "extreme-contrast" = "Bias: extreme-contrast"
)
)

# -------------------------------
# 1) Figure 1 — Efficiency frontiers (Empirical SE)
#    Two slices: R2 = 0.85 and R2 = 0.40; small-multiples by residual variance.
# -------------------------------
make_frontier_plot <- function(df, r2_val, bias_vals=c("none","moderate","large","extreme-contrast")) {
    d <- df %>%
        filter(R2_target == r2_val, method %in% c("Subset", "SRS", "Stratified", "Optimal"),
               bias %in% bias_vals)

    ref = df %>%
        filter(R2_target == r2_val, method=="SRS", bias %in% bias_vals) %>%
        group_by(bias, residual_var, R2_target) %>%
        summarise(full_coding_se=mean(emp_se_full))

        ggplot(d, aes(x = coding_frac, y = emp_se, color = method)) +
            geom_line() + geom_point(size = 1) +
            facet_grid(bias ~ residual_var, labeller = labs) +
            scale_color_manual(values = method_pal) +
            scale_x_continuous(breaks=seq(0.1,0.9,by=0.1), labels = percent_format(accuracy = 1)) +
            geom_hline(data=ref, aes(yintercept=full_coding_se), lty=2)+
            xlab("Proportion coded")+ylab("Empirical SE of impact estimate")+
            ggtitle(glue::glue("Efficiency frontiers (R² target = {r2_val})")) +
            theme_minimal() +
            theme(legend.position = "bottom", panel.grid.minor=element_blank() )
}

fig1_high <- make_frontier_plot(main, r2_val = 0.85)
fig1_low  <- make_frontier_plot(main, r2_val = 0.40)

ggsave("figures/fig1_frontier_highR2.pdf", fig1_high, width = 8.5, height = 6.5)
ggsave("figures/fig1_frontier_lowR2.pdf",  fig1_low,  width = 8.5, height = 6.5)

# -------------------------------
# 2) Figure 2 — Incremental value of stratification
#    ∆(% variance reduction vs SRS) across budgets and moderators
# -------------------------------
make_delta_plot <- function(df, r2_val, bias_vals=c("none","moderate","large","extreme-contrast")) {
    base <- df %>%
        filter(R2_target == r2_val, method %in% c("Stratified","Optimal"), bias%in%bias_vals) %>%
        select(all_of(cell_keys), method, pct_var_red_vs_SRS)

    ggplot(base, aes(x = coding_frac, y = pct_var_red_vs_SRS, color = method)) +
        facet_grid(bias~residual_var, labeller = labs)+
        geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        geom_line() + geom_point(size = 1) +
        scale_color_manual(values = c("Stratified"="#6baed6","Optimal"="#fb6a4a"), labels=c("Proportional Allocation", "Optimal Allocation")) +
        scale_x_continuous(breaks=seq(0.1,0.9,by=0.1), labels = percent_format(accuracy = 1)) +
        scale_y_continuous(breaks=seq(0,80,by=20), labels = percent_format(scale=1),limits=c(-12,85)) +
        xlab("Proportion coded")+ ylab("% Reduction in Empirical Variance vs. SRS")+
        ggtitle(glue::glue("Incremental value of stratified sampling vs SRS (R² target = {r2_val})") ) +
        theme_minimal() +
        theme(legend.position = "bottom")
}


fig2_high <- make_delta_plot(main, r2_val = 0.85, bias_vals=c("none","moderate","large"))
fig2_low  <- make_delta_plot(main, r2_val = 0.40, bias_vals=c("none","moderate","large"))

ggsave("figures/fig2_delta_opt_over_strat_highR2.pdf", fig2_high, width = 7.2, height = 6.0)
ggsave("figures/fig2_delta_opt_over_strat_lowR2.pdf",  fig2_low,  width = 7.2, height = 6.0)

# -------------------------------
# 3) Figure 3 — Heatmaps at fixed budgets (e.g., 20% and 40%)
#    (a) % reduction vs SRS for Stratified and Optimal
#    (b) Incremental gain (Optimal − Stratified)
# -------------------------------
fixed_budgets <- c(0.20, 0.40)  # modify as desired

pick_nearest <- function(x, target, tol = 0.02) {
    out <- x[which.min(abs(x - target))]
    if (abs(out - target) > tol) NA_real_ else out
}

make_heatmap = function(dat, r2_val, fixed_budget) {

    # Collapse over bias by averaging, but keep bias range for annotation if needed
    heat_df <- dat %>%
        filter(method %in% c("Stratified","Optimal"), R2_target==r2_val) %>%
        group_by(residual_var, bias, coding_frac, method) %>%
        summarise(pct_red = mean(pct_var_red_vs_SRS, na.rm=TRUE), .groups="drop") %>%
        mutate(bias = fct_relevel(bias, "none", "moderate", "large", "extreme-contrast"),
               residual_var=fct_relevel(residual_var, "homogeneous","heterogeneous", "extreme-contrast"))


    # Choose nearest available grid points to the requested budgets
    avail <- heat_df %>%
        distinct(coding_frac) %>%
        arrange(coding_frac) %>%
        pull(coding_frac)

    chosen <- map_dbl(fixed_budget, ~ pick_nearest(avail, .x))

    heat_df %>%
            filter(coding_frac == chosen) %>%
            ggplot(aes(x = residual_var, y = bias, fill = pct_red)) +
            geom_tile(color = "white") +
            facet_wrap(~ method) +
            scale_fill_distiller(palette = "Blues", direction = 1, oob = scales::squish, limits=c(0,80),
                                 name = "% var. reduction vs SRS", labels=percent_format(scale=1)) +
            xlab("Residual variance pattern")+ ylab("Bias pattern")+
            ggtitle(glue::glue("Heatmap of % variance reduction vs SRS at {percent(fixed_budget)} coded (R² target = {r2_val})")) +
            theme_minimal() + theme(legend.position = "bottom")

}

fig3_high <- make_heatmap(main, r2_val = 0.85, fixed_budget=0.2)
fig3_low <- make_heatmap(main, r2_val = 0.4, fixed_budget=0.2)


ggsave(glue::glue("figures/fig3_heatmap_20pct_highR2.pdf"), fig3_high, width = 8, height = 3.5)
ggsave(glue::glue("figures/fig3_heatmap_20pct_lowR2.pdf"), fig3_low, width = 8, height = 3.5)


# -------------------------------
# 4) Figure 4 — Variance inflation vs FULL CODING (if present)
#     y = 100 * Var(method) / Var(full-coding). 100% = as good as full coding.
# -------------------------------

    make_full_ratio_plot <- function(df, r2_val) {
        d <- df %>%
            filter(R2_target == r2_val, method %in% c("Subset","SRS","Stratified","Optimal","Full coding"))
        ggplot(d, aes(x = coding_frac, y = re_vs_fullcoding, color = method)) +
            geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
            geom_line() + geom_point(size = 1) +
            facet_grid(bias~ residual_var,  labeller = labs) +
            scale_color_manual(values = method_pal) +
            scale_x_continuous(breaks=seq(0.1,0.9,by=0.1), labels = percent_format(accuracy = 1)) +
            scale_y_continuous(limits=c(0.9,10), oob=scales::squish)+
            xlab("Proportion coded")+ylab("Variance Ratio")+
            ggtitle(glue::glue("Variance inflation vs full coding (R² target = {r2_val})"))+
            theme_minimal() +
            theme(legend.position = "bottom")
    }
    fig4_high <- make_full_ratio_plot(main, 0.85)
    fig4_low  <- make_full_ratio_plot(main, 0.40)
    ggsave("figures/fig4_vs_full_highR2.pdf", fig4_high, width = 8.5, height = 6.5)
    ggsave("figures/fig4_vs_full_lowR2.pdf",  fig4_low,  width = 8.5, height = 6.5)


## Figure 5

res<- analysis_results %>%
    mutate(pct_coded=n_code/N,
           any_bias=ifelse(bias=="none",0,1),
           any_hetero=ifelse(residual_var=="homogeneous",0,1)) %>%
    group_by(any_bias, any_hetero, R2_target, n_code, pct_coded, method) %>%
    summarize(emp_se = mean(emp_se, na.rm=TRUE), .groups="drop") %>%
    mutate(any_bias = factor(any_bias, levels=c(0,1), labels=c("No Bias", "Biased")),
           any_hetero=factor(any_hetero, levels=c(0,1), labels=c("Var: Homogenous", "Var: Heterogeneous")))


res$bias.var.config = sprintf("%s, %s Var", res$any_bias, res$any_hetero)

res.oracle=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_oracle.rds")
oracle_se=sd(res.oracle$est)



fig5 <-
        ggplot(res, aes(x = pct_coded, y = emp_se, color = method)) +
            geom_smooth(se=F) +
            geom_hline(yintercept=oracle_se, lty=2)+
            facet_grid(paste("Target R²=", R2_target)~any_bias+any_hetero) +
            scale_x_continuous(breaks=seq(0,1,by=0.2), labels = scales::percent_format(accuracy = 1),limits=c(0.1,0.95)) +
            scale_y_continuous(limits=c(0.175,0.5),oob=scales::squish) +
            scale_color_manual(values = method_pal) +

            ggplot2::labs(x="Coding budget", y="Empirical SE",
                          title="Average Efficiency Frontiers across Bias/Residual Structures",
                          color="Method") +
            theme_minimal()+theme(legend.position = "bottom",panel.grid.minor=element_blank())
ggsave("figures/fig5_avg_se.pdf", fig5, width = 8.5, height = 6.5)


# -------------------------------
# 5) Table 1 — Regression summary of main effects & interactions
#     Outcome: log(empirical SE)
# -------------------------------
reg_df <- analysis_results %>%
    filter(method %in% c("SRS","Stratified","Optimal")) %>%  # use SRS as baseline
    mutate(method = fct_relevel(method, "SRS"),
           residual_var=fct_relevel(residual_var, "homogeneous"),
           bias=fct_relevel(bias, "none"))

fit1 <- lm(log(emp_se) ~ coding_frac  +  as.factor(R2_target) + method + strata_config + bias + residual_var, data = reg_df)
tab1 <- broom::tidy(fit1, conf.int = TRUE) %>%
    mutate(
        pct_change = 100 * (exp(estimate) - 1),
        pct_LL     = 100 * (exp(conf.low) - 1),
        pct_UL     = 100 * (exp(conf.high) - 1)
    ) %>%
    select(term, estimate, std.error, pct_change, pct_LL, pct_UL, p.value) %>%
    arrange(desc(abs(pct_change))) %>%
    mutate(term = recode(term,
                         "(Intercept)" = "Intercept",
                         "coding_frac" = "Coding fraction",
                         "biasmoderate" = "Bias: Moderate (vs None)",
                         "biaslarge" = "Bias: Large (vs None)",
                         "biasextreme-contrast" = "Bias: Extreme contrast (vs None)",
                         "methodStratified" = "Method: Stratified (vs SRS)",
                         "methodOptimal" = "Method: Optimal (vs SRS)",
                         "residual_varheterogeneous" = "Residual variance: heterogeneous (vs homogeneous)",
                         "residual_varextreme-contrast" = "Residual variance: extreme-contrast (vs homogeneous)",
                         "strata_configbalanced-approx" = "Balanced approx strata (vs. balanced exact)",
                         "strata_configunbalanced" = "Unbalanced strata (vs. balanced exact)"
    )) %>% filter(term!="Intercept")

tab1_gt <- tab1 %>%
    mutate_at(vars(estimate, std.error), ~sprintf("%.2f", .)) %>%
    mutate(across(starts_with("pct"), ~ sprintf("%.1f%%", .))) %>%
    mutate(sig= cut(p.value, breaks=c(0, 0.001, 0.01, 0.05, 1), labels=c("<0.001", "<0.01", "<0.05","NS"),include.lowest=T)) %>%
    mutate(p.val = ifelse(p.value>0.05, sprintf("%.2f",p.value), as.character(sig))) %>%
    select(-p.value, -sig) %>%
    gt() %>%
    tab_header(title = "Regression on log(empirical SE): main effects",
               subtitle = "Coefficients shown as % change in SE (SRS as baseline)")

gtsave(tab1_gt, "figures/table1a_regression_main_effects.tex")






fit2 <- lm(log(emp_se) ~ coding_frac  +  as.factor(R2_target) + method * bias + method*residual_var + strata_config, data = reg_df)
tab2 <- broom::tidy(fit2, conf.int = TRUE) %>%
    mutate(
        pct_change = 100 * (exp(estimate) - 1),
        pct_LL     = 100 * (exp(conf.low) - 1),
        pct_UL     = 100 * (exp(conf.high) - 1)
    ) %>%
    select(term, estimate, std.error, pct_change, pct_LL, pct_UL, p.value) %>%
    arrange(desc(abs(pct_change))) %>%
    mutate(term = recode(term,
                         "(Intercept)" = "Intercept",
                         "coding_frac" = "Coding fraction",
                         "biasmoderate" = "Bias: Moderate (vs None)",
                         "biaslarge" = "Bias: Large (vs None)",
                         "biasextreme-contrast" = "Bias: Extreme contrast (vs None)",
                         "methodStratified" = "Method: Stratified (vs SRS)",
                         "methodOptimal" = "Method: Optimal (vs SRS)",
                         "residual_varheterogeneous" = "Residual variance: heterogeneous (vs homogeneous)",
                         "residual_varextreme-contrast" = "Residual variance: extreme-contrast (vs homogeneous)",
                         "strata_configbalanced-approx" = "Balanced approx strata (vs. balanced exact)",
                         "strata_configunbalanced" = "Unbalanced strata (vs. balanced exact)"
    )) %>% filter(term!="Intercept")

tab2_gt <- tab2 %>%
    mutate_at(vars(estimate, std.error), ~sprintf("%.2f", .)) %>%
    mutate(across(starts_with("pct"), ~ sprintf("%.1f%%", .))) %>%
    mutate(sig= cut(p.value, breaks=c(0, 0.001, 0.01, 0.05, 1), labels=c("<0.001", "<0.01", "<0.05","NS"),include.lowest=T)) %>%
    mutate(p.val = ifelse(p.value>0.05, sprintf("%.2f",p.value), as.character(sig))) %>%
    select(-p.value, -sig) %>%
    gt() %>%
    tab_header(title = "Regression on log(empirical SE): main effects and key interactions",
               subtitle = "Coefficients shown as % change in SE (SRS as baseline)")

gtsave(tab2_gt, "figures/table1b_regression_main_and_interaction_effects.tex")

# -------------------------------
# 6) Table 2 — Design-relevant contrasts at fixed budgets
#     (10%, 20%, 40%): % reduction vs SRS for Stratified and Optimal, and their difference
# -------------------------------
budget_targets <- c(0.10, 0.20, 0.40)

nearest <- function(v, target) v[which.min(abs(v - target))]
avail_budgets <- sort(unique(main$coding_frac))



get_contrast_table = function(dat, budget_targets, r2_val, bias_vals=c("none","moderate","large","extreme-contrast")){
    contrast_rows <- map_dfr(budget_targets, function(b) {
        b_star <- nearest(avail_budgets, b)
        dat %>%
            filter(coding_frac == b_star, method %in% c("Stratified","Optimal"), R2_target==r2_val) %>%
            group_by(bias,residual_var, method) %>%
            summarise(pct_red = mean(pct_var_red_vs_SRS, na.rm=TRUE), .groups="drop") %>%
            pivot_wider(names_from = method, values_from = pct_red) %>%
            mutate(
                budget = b_star,
                diff = Optimal - Stratified
            )
    })

    contrast_rows %>%
        mutate(
            budget = scales::percent(budget, accuracy = 1),
            across(c("Stratified","Optimal","diff"), ~ sprintf("%.1f", .))
        ) %>%
        arrange(budget, bias,residual_var) %>%
        select(bias, residual_var, budget, Stratified, Optimal, diff) %>%
        gt() %>%
        tab_header(title = "Design-relevant contrasts",
                   subtitle = glue::glue("% variance reduction vs SRS at selected coding budgets (balanced-exact, R² target = {r2_val})")) %>%
        cols_label(
            bias="Bias pattern",
            residual_var = "Residual variance pattern",
            budget = "Budget",
            Stratified = "Stratified",
            Optimal = "Optimal",
            diff = "Incremental"
        )

}

gt=get_contrast_table(main, budget_targets = c(0.1,0.2,0.4), r2_val=0.85,bias_vals=c("none","moderate","large"))
gtsave(gt, "figures/table2_design_contrasts.html")




# Fit linear model on log empirical SE for stability
lm_fit <- lm(log(emp_se) ~ bias + residual_var + as.factor(R2_target) + coding_frac + method,
             data = analysis_results[analysis_results$method!="Subset",])

importance_df <- broom::tidy(lm_fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(pct_change = 100*(exp(estimate)-1))

fig6=ggplot(importance_df, aes(x=reorder(term, abs(estimate)), y=estimate)) +
    geom_col() +
    coord_flip() +
    ggplot2::labs(x="Factor", y="Coefficient (log Emp SE scale)",
                  title="Relative Influence of Design Factors (Baseline: SRS with no bias & homogeneous variance)") +
    theme_minimal()+theme(plot.title=element_text(size=rel(0.9)),plot.title.position = "plot")
ggsave("figures/fig6_importance_plot.pdf", fig6, width = 6.5, height = 6)


