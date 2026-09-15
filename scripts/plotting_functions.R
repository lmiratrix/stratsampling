library(tidyverse)
library(scales)
# Canonical method ordering (single source of truth for factor levels)
method_levels <- c(
    "Subset",
    "SRS",
    "Stratified",
    "Stratified-Random",
    "Optimal",
    "Optimal-Permuted",
    "Optimal-Random",
    "Optimal-Reversed",
    "Optimal-Pilot-Pooled",          # Case 2: strata known, pilot estimates variances, pilot pooled
    "Optimal-Pilot-Ind",             # Case 2, pilot discarded (robustness arm)
    "Stratified-PilotSelect",        # Case 1: SRS pilot selects the strata; proportional wave 2
    "Optimal-PilotSelect",           # Case 1, Neyman wave 2 on EB-shrunk pilot variances (Sep 15 sim_functions.R)
    "Pilot+Stratified-PilotSelect",  # Case 1, proportional wave 2 precision-combined with the pilot estimate
    "Pilot+Optimal-PilotSelect",     # Case 1, Neyman wave 2 precision-combined with the pilot estimate
    "Optimal-PilotSelect-Shrink",    # legacy label (pre-Sep-15 results files only)
    "Full coding"
)

# Human-readable legend labels (use via method_scales(use_labels = TRUE);
# the existing plot helpers keep the raw method names so that figures
# already in the manuscript do not change).
method_labels <- c(
    "Subset"                     = "Subset (coded units only)",
    "SRS"                        = "SRS",
    "Stratified"                 = "Proportional allocation",
    "Stratified-Random"          = "Proportional (random strata)",
    "Optimal"                    = "Optimal (oracle variances)",
    "Optimal-Permuted"           = "Optimal (permuted variances)",
    "Optimal-Random"             = "Optimal (random variances)",
    "Optimal-Reversed"           = "Optimal (reversed variances)",
    "Optimal-Pilot-Pooled"       = "Optimal (pilot variances, pooled)",
    "Optimal-Pilot-Ind"          = "Optimal (pilot variances, discarded)",
    "Stratified-PilotSelect"     = "Proportional (pilot-selected strata)",
    "Optimal-PilotSelect"        = "Optimal (pilot-selected strata)",
    "Pilot+Stratified-PilotSelect" = "Proportional (pilot-selected strata) + pilot estimate",
    "Pilot+Optimal-PilotSelect"    = "Optimal (pilot-selected strata) + pilot estimate",
    "Optimal-PilotSelect-Shrink" = "Optimal (pilot-selected strata, shrunk var.; legacy)",
    "Full coding"                = "Full coding"
)

# Groups used to subset the registry in the new pilot figures
pilot_case2_methods <- c("Optimal-Pilot-Pooled", "Optimal-Pilot-Ind")
pilot_case1_methods <- c("Stratified-PilotSelect", "Optimal-PilotSelect",
                         "Pilot+Stratified-PilotSelect", "Pilot+Optimal-PilotSelect")
pilot_methods       <- c(pilot_case2_methods, pilot_case1_methods)

# A consistent color palette for methods
method_pal <- c(
    "Subset"      = "black", #"#9467bd",
    "SRS"         = "#228222",
    "Stratified"  = "#1f77b4",
    "Stratified-Random" = "purple3",
    "Optimal"     = "#ff6961",
    "Optimal-Permuted"  = "#e6842a",
    "Optimal-Random"  = "#8c1a11",
    "Optimal-Reversed" = "goldenrod",
    "Optimal-Pilot-Pooled"       = "#17becf",
    "Optimal-Pilot-Ind"          = "#7f7f7f",
    "Stratified-PilotSelect"     = "#08306b",
    "Optimal-PilotSelect"        = "#762a83",
    "Pilot+Stratified-PilotSelect" = "#4292c6",
    "Pilot+Optimal-PilotSelect"    = "#c51b7d",
    "Optimal-PilotSelect-Shrink" = "#762a83",
    "Full coding" = "black"
)

method_shapes <- c(
    "Subset"      = 4,
    "SRS"         = 16,
    "Stratified"  = 17,
    "Optimal"     = 15,
    "Stratified-Random" = 2,
    "Optimal-Permuted"  = 0,
    "Optimal-Random" = 7,
    "Optimal-Reversed"  = 7,
    "Optimal-Pilot-Pooled"       = 18,
    "Optimal-Pilot-Ind"          = 5,
    "Stratified-PilotSelect"     = 6,
    "Optimal-PilotSelect"        = 13,
    "Pilot+Stratified-PilotSelect" = 25,
    "Pilot+Optimal-PilotSelect"    = 3,
    "Optimal-PilotSelect-Shrink" = 13,
    "Full coding" = 8
)

# Line type palette (distinct line types for each method)
method_linetypes <- c(
    "Subset"      = "dotted",
    "SRS"         = "solid",
    "Stratified"  = "dashed",
    "Optimal"     = "longdash",
    "Stratified-Random" = "twodash",
    "Optimal-Permuted"  = "42",
    "Optimal-Random" = "4212",
    "Optimal-Reversed"  = "4212",
    "Optimal-Pilot-Pooled"       = "solid",
    "Optimal-Pilot-Ind"          = "twodash",
    "Stratified-PilotSelect"     = "1343",
    "Optimal-PilotSelect"        = "solid",
    "Pilot+Stratified-PilotSelect" = "1343",
    "Pilot+Optimal-PilotSelect"    = "3313",
    "Optimal-PilotSelect-Shrink" = "solid",
    "Full coding" = "dotdash"
)

# One call that adds the colour / shape / linetype scales for methods.
# use_labels = TRUE swaps the raw method names for method_labels in the
# legend (new pilot figures); the default keeps the legacy legends.
method_scales <- function(use_labels = FALSE, name = "Method", nrow = NULL) {
    lab <- if (use_labels) method_labels else waiver()
    list(
        scale_color_manual(name = name, values = method_pal, labels = lab, drop = TRUE),
        scale_shape_manual(name = name, values = method_shapes, labels = lab, drop = TRUE),
        scale_linetype_manual(name = name, values = method_linetypes, labels = lab, drop = TRUE),
        if (!is.null(nrow)) guides(color = guide_legend(nrow = nrow),
                                   shape = guide_legend(nrow = nrow),
                                   linetype = guide_legend(nrow = nrow))
    )
}



# Helper for clean labels.  Lookup labellers pass unmatched values
# through unchanged (a plain named vector turns them into NA, which is
# what happened to the aggregated "No Bias" / "Biased" facets).
passthrough <- function(map) {
    function(x) { x <- as.character(x); out <- unname(map[x]); ifelse(is.na(out), x, out) }
}
map_resid <- c("homogeneous" = "Residual variance: homogeneous",
               "heterogeneous" = "Residual variance: heterogeneous",
               "extreme-contrast" = "Residual variance: extreme-contrast")
map_bias  <- c("none" = "Bias: none", "small" = "Bias: small", "moderate" = "Bias: moderate",
               "large" = "Bias: large", "extreme-contrast" = "Bias: extreme-contrast")
map_R2    <- c("0.4" = "R\u00b2 target = 0.4", "0.85" = "R\u00b2 target = 0.85")
map_config <- c("balanced-exact" = "Strata: balanced-exact", "balanced-approx" = "Strata: balanced-approx",
                "unbalanced" = "Strata: unbalanced")

lab_resid <- labeller(residual_var = passthrough(map_resid))
lab_bias  <- labeller(bias = passthrough(map_bias))

# (named labs_all, not labs: a global called labs masks ggplot2::labs())
labs_all <- labeller(residual_var = passthrough(map_resid),
                     bias = passthrough(map_bias),
                     R2_target = passthrough(map_R2),
                     strata_config = passthrough(map_config))

# -------------------------------
# 1) Figure 1 — Efficiency frontiers (Empirical SE)
#    Two slices: R2 = 0.85 and R2 = 0.40; small-multiples by residual variance.
# -------------------------------
make_frontier_plot <- function(df, r2_val = NULL, facet = "bias~residual_var",
                               bias_vals = NULL, oracle_ref = FALSE, smooth = FALSE,
                               methods = c("Subset", "SRS", "Stratified", "Optimal")) {
    d <- df %>%
        filter(method %in% methods)
    
    if (!is.null(bias_vals)) {
        d <- d %>% filter(bias %in% bias_vals)
    }
    
    if (!is.null(r2_val)) {
        d <- d %>% filter(R2_target == r2_val)
    }
    
    if (oracle_ref) {
        ref <- df %>%
            filter(method == "SRS") %>%
            group_by(any_bias, any_hetero, bias.var.config, R2_target) %>%
            summarise(full_coding_se = mean(emp_se_full), .groups = "drop")
    }
    
    # Base aesthetic mapping with color, shape, and linetype
    base_aes <- aes(x = coding_frac, y = emp_se, color = method,
                    shape = method, linetype = method)
    
    if (smooth && oracle_ref) {
        out <- ggplot(d, base_aes) +
            geom_point(size = 2) +
            geom_smooth(se = FALSE) +
            facet_grid(facet, labeller = labs_all) +
            geom_hline(data = ref, aes(yintercept = full_coding_se),
                       lty = 2, col = "grey50")
    } else if (!smooth && oracle_ref) {
        out <- ggplot(d, base_aes) +
            geom_line(linewidth = 0.9) +
            geom_point(size = 2.5) +
            facet_grid(facet, labeller = labs_all) +
            geom_hline(data = ref, aes(yintercept = full_coding_se),
                       lty = 2, col = "grey50")
    } else if (smooth && !oracle_ref) {
        out <- ggplot(d, base_aes) +
            geom_point(size = 2) +
            geom_smooth(se = FALSE) +
            facet_grid(facet, labeller = labs_all)
    } else {
        # Default: lines with points (colorblind-friendly)
        out <- ggplot(d, base_aes) +
            geom_line(linewidth = 0.9) +
            geom_point(size = 2.5) +
            facet_grid(facet, labeller = labs_all)
    }
    
    # Apply colorblind-friendly scales (color + shape + linetype)
    out <- out +
        scale_color_manual(name = "Method", values = method_pal) +
        scale_shape_manual(name = "Method", values = method_shapes) +
        scale_linetype_manual(name = "Method", values = method_linetypes) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0.05, 0.95),
                           labels = percent_format(accuracy = 1)) +
        xlab("Proportion coded") +
        ylab("Empirical Std. Error") +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank())
    
    out
}



make_delta_plot <- function(df, facet = "bias~residual_var", r2_val = NULL,
                            bias_vals = NULL, smooth = FALSE,
                            methods=c("Stratified", "Optimal")) {
    base <- df %>%
        filter(method %in% methods) %>%
        #mutate(method = fct_relevel(method, "Stratified", "Optimal")) %>%
        select(all_of(cell_keys), method, pct_var_red_vs_SRS)
    
    if (!is.null(bias_vals)) {
        base <- base %>% filter(bias %in% bias_vals)
    }
    
    if (!is.null(r2_val)) {
        base <- base %>% filter(R2_target == r2_val)
    }
    
    # Base aesthetic with color, shape, and linetype
    base_aes <- aes(x = coding_frac, y = pct_var_red_vs_SRS,
                    color = method, shape = method, linetype = method)
    
    if (smooth) {
        g <- ggplot(base, base_aes) +
            facet_grid(facet, labeller = labs_all) +
            geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
            geom_point(size = 2) +
            geom_smooth(se = FALSE)
    } else {
        g <- ggplot(base, base_aes) +
            facet_grid(facet, labeller = labs_all) +
            geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
            geom_line(linewidth = 1.0) +
            geom_point(size = 2.5)
    }
    
    # Apply colorblind-friendly scales
    g <- g +
        scale_color_manual(name = "Method", values = method_pal) +
        scale_shape_manual(name = "Method", values = method_shapes) +
        scale_linetype_manual(name = "Method", values = method_linetypes) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                           labels = percent_format(accuracy = 1), limits = c(0.05, 0.95)) +
        scale_y_continuous(breaks = seq(-100, 100, by = 20),oob=scales::squish,
                           labels = percent_format(scale = 1)
                           , #limits = c(-10, 85)
        )+
        xlab("Proportion coded") +
        ylab("% Reduction in Empirical Variance vs. SRS") +
        ggtitle(glue::glue("Incremental value of stratified sampling vs SRS (R-squared target = {r2_val})")) +
        theme_minimal() +
        theme(legend.position = "bottom")
    
    return(g)
}




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



make_full_ratio_plot <- function(df, facet = "bias~residual_var", r2_val = NULL, bias_vals = NULL,
                                 methods=c("Subset", "SRS", "Stratified", "Optimal", "Full coding")) {
    d <- df %>%
        filter(method %in% methods)
    
    title <- "Variance inflation vs full coding"
    
    if (!is.null(r2_val)) {
        d <- d %>% filter(R2_target == r2_val)
        title <- paste0(title, glue::glue(" (R-squared target = {r2_val})"))
    }
    
    if (!is.null(bias_vals)) {
        d <- d %>% filter(bias %in% bias_vals)
    }
    
    max.y <- max(d$re_vs_fullcoding)
    y.ul <- ifelse(max.y > 10, 10, ceiling(max.y))
    
    # Base aesthetic with color, shape, and linetype
    ggplot(d, aes(x = coding_frac, y = re_vs_fullcoding,
                  color = method, shape = method, linetype = method)) +
        geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
        geom_line(linewidth = 1.0) +
        geom_point(size = 2.5) +
        facet_grid(facet, labeller = labs_all) +
        scale_color_manual(name = "Method", values = method_pal) +
        scale_shape_manual(name = "Method", values = method_shapes) +
        scale_linetype_manual(name = "Method", values = method_linetypes) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0.05, 0.95),
                           labels = percent_format(accuracy = 1)) +
        xlab("Proportion coded") +
        ylab("Variance Ratio") +
        ggplot2::ggtitle(title) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank())
}



make_coverage_plot <- function(df, facet = "bias~residual_var", r2_val = NULL, bias_vals = NULL,
                               methods= c("Subset", "SRS", "Stratified", "Optimal")) {
    d <- df %>%
        filter(method %in% methods)
    
    if (!is.null(bias_vals)) {
        d <- d %>% filter(bias %in% bias_vals)
    }
    
    if (!is.null(r2_val)) {
        d <- d %>% filter(R2_target == r2_val)
    }
    
    ggplot(d, aes(x = coding_frac, y = coverage, color = method)) +
        geom_hline(yintercept = 0.95, linetype = 2, color = "grey50") +
        geom_line() + geom_point(size = 1) +
        facet_grid(facet, labeller = labs_all) +
        scale_color_manual(name = "Method", values = method_pal) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2), 
                           labels = percent_format(accuracy = 1), 
                           limits = c(0.05, 0.95)) +
        scale_y_continuous(breaks = seq(0.85, 1.0, by = 0.05), 
                           labels = percent_format(accuracy = 1),
                           limits = c(0.85, 1.0)) +
        xlab("Proportion coded") + ylab("95% CI Coverage Rate") +
        ggtitle(glue::glue("Coverage rates across simulation scenarios{ifelse(!is.null(r2_val), paste0(' (R² target = ', r2_val, ')'), '')}")) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank())
}


make_bias_plot <- function(df, facet = "bias~residual_var", r2_val = NULL, bias_vals = NULL,
                           methods = c("Subset", "SRS", "Stratified", "Optimal")) {
    d <- df %>%
        filter(method %in% methods)
    
    if (!is.null(bias_vals)) {
        d <- d %>% filter(bias %in% bias_vals)
    }
    
    if (!is.null(r2_val)) {
        d <- d %>% filter(R2_target == r2_val)
    }
    
    ggplot(d, aes(x = coding_frac, y = bias_est, color = method)) +
        geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
        geom_line() + geom_point(size = 1) +
        facet_grid(facet, labeller = labs_all) +
        scale_color_manual(name = "Method", values = method_pal) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2), 
                           labels = percent_format(accuracy = 1), 
                           limits = c(0.05, 0.95)) +
        scale_y_continuous(limits = c(-0.1, 0.1)) +
        xlab("Proportion coded") + ylab("Monte Carlo Bias") +
        ggtitle(glue::glue("Estimator bias across simulation scenarios{ifelse(!is.null(r2_val), paste0(' (R² target = ', r2_val, ')'), '')}")) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank())
}