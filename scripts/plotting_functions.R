library(tidyverse)
library(scales)
# Canonical method ordering (single source of truth for factor levels)
method_levels <- c(
    "Subset",
    "SRS",
    "Stratified",
    "Stratified-Random",
    "Stratified-Yhat",
    "Optimal",
    "Optimal-Permuted",
    "Optimal-Random",
    "Optimal-Reversed",
    "Pilot-Allocated-Pooled",          # Case 2: strata known, pilot estimates variances, pilot pooled
    "Pilot-Allocated-Discard",             # Case 2, pilot discarded (robustness arm)
    "Pilot-Selected-Prop",        # Case 1: SRS pilot selects the strata; proportional wave 2
    "Pilot-Selected-Optimal",           # Case 1: SRS pilot selects strata; Neyman wave 2 on EB-shrunk pilot variances
    "Full coding",
    "Oracle"                         # application scripts' name for full coding (same colour)
)

# Human-readable legend labels (use via method_scales(use_labels = TRUE);
# the existing plot helpers keep the raw method names so that figures
# already in the manuscript do not change).
method_labels <- c(
    "Subset"                     = "Subset",
    "SRS"                        = "SRS",
    "Stratified"                 = "Stratified",
    "Stratified-Random"          = "Stratified (random strata)",
    "Optimal"                    = "Optimal (oracle variances)",
    "Stratified-Yhat"   = "Stratified (Yhat quantile strata)",  # method_labels
    "Optimal-Permuted"           = "Optimal (permuted variances)",
    "Optimal-Random"             = "Optimal (random variances)",
    "Optimal-Reversed"           = "Optimal (reversed variances)",
    "Pilot-Allocated-Pooled"       = "Pilot-allocated (pooled)",
    "Pilot-Allocated-Discard"          = "Pilot-allocated (discarded)",
    "Pilot-Selected-Prop"     = "Pilot-selected (proportional allocation)",
    "Pilot-Selected-Optimal"        = "Pilot-selected (optimal allocation)",
    "Full coding"                = "Full coding",
    "Oracle"                     = "Oracle"
)

# Groups used to subset the registry in the new pilot figures
pilot_case2_methods <- c("Pilot-Allocated-Pooled", "Pilot-Allocated-Discard")
pilot_case1_methods <- c("Pilot-Selected-Prop", "Pilot-Selected-Optimal")
pilot_methods       <- c(pilot_case2_methods, pilot_case1_methods)

# A consistent color palette for methods
method_pal <- c(
    "Subset"                  = "black",
    "SRS"                     = "#666666",
    "Stratified"              = "#1f77b4",   # true-strata family (blue)
    "Stratified-Random"       = "#1f77b4",
    "Optimal"                 = "#EE6677",   # oracle-optimal family (red)
    "Optimal-Permuted"        = "#EE6677",
    "Optimal-Random"          = "#EE6677",
    "Optimal-Reversed"        = "#8C1A11",   # darkest red: the anti-informative floor
    "Stratified-Yhat"         = "#228833",   # Yhat-only family (green)
    "Pilot-Allocated-Pooled"  = "#B8860B",   # pilot-allocated family (gold)
    "Pilot-Allocated-Discard" = "#B8860B",
    "Pilot-Selected-Prop"     = "#7B3294",   # pilot-selected family (purple)
    "Pilot-Selected-Optimal"  = "#7B3294",
    "Full coding" = "black", "Oracle" = "black")


method_shapes <- c(
    "Subset" = 4, "SRS" = 16,
    "Stratified" = 17, "Stratified-Random" = 2,
    "Optimal" = 15, "Optimal-Permuted" = 0, "Optimal-Random" = 1, "Optimal-Reversed" = 7,
    "Stratified-Yhat" = 18,
    "Pilot-Allocated-Pooled" = 16, "Pilot-Allocated-Discard" = 0,
    "Pilot-Selected-Prop" = 17, "Pilot-Selected-Optimal" = 5,
    "Full coding" = 8, "Oracle" = 8)
method_linetypes <- c(
    "Subset" = "dotted", "SRS" = "solid",
    "Stratified" = "solid", "Stratified-Random" = "dotted",
    "Optimal" = "longdash", "Optimal-Permuted" = "dotted", "Optimal-Random" = "4212", "Optimal-Reversed" = "42",
    "Stratified-Yhat" = "dotdash",
    "Pilot-Allocated-Pooled" = "solid", "Pilot-Allocated-Discard" = "22",
    "Pilot-Selected-Prop" = "solid", "Pilot-Selected-Optimal" = "22",
    "Full coding" = "dotdash", "Oracle" = "dotdash")

# One call that adds the colour / shape / linetype scales for methods.
# use_labels = TRUE swaps the raw method names for method_labels in the
# legend (new pilot figures); the default keeps the legacy legends.
method_scales <- function(use_labels = FALSE, name = "Method:", nrow = NULL) {
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
        scale_color_manual(name = "Method:", values = method_pal) +
        scale_shape_manual(name = "Method:", values = method_shapes) +
        scale_linetype_manual(name = "Method:", values = method_linetypes) +
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
        filter(method %in% methods)
    
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
        scale_color_manual(name = "Method:", values = method_pal) +
        scale_shape_manual(name = "Method:", values = method_shapes) +
        scale_linetype_manual(name = "Method:", values = method_linetypes) +
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
        theme(legend.position = "bottom", panel.grid.minor=element_blank())
    
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
        scale_color_manual(name = "Method:", values = method_pal) +
        scale_shape_manual(name = "Method:", values = method_shapes) +
        scale_linetype_manual(name = "Method:", values = method_linetypes) +
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
        scale_color_manual(name = "Method:", values = method_pal) +
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



make_pilot_plot <- function(dat, refs, pilot, facet, ylim = c(-30, 65)) {
    show <- c(refs, pilot)
    d <- dat %>%
        filter(method %in% show) %>%
        mutate(method = factor(as.character(method), levels = show))
    ggplot(d, aes(coding_frac, pct_var_red_vs_SRS,
                  colour = method, linetype = method, shape = method)) +
        geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
        geom_line(data = ~ filter(.x, method %in% pilot), linewidth = 0.8) +
        geom_point(data = ~ filter(.x, method %in% pilot), size = 2.1, fill = "white", stroke = 0.9) +  # pilot arms only
        geom_line(data = ~ filter(.x, method %in% refs), linewidth = 1.4, alpha = 0.75) +            # references: line only, no symbols
        facet_grid(facet, labeller = labs_all) +
        scale_colour_manual(values = method_pal[show], breaks = show,
                            labels = method_labels[show], name = NULL) +
        scale_linetype_manual(values = method_linetypes[show], breaks = show,
                              labels = method_labels[show], name = NULL) +
        scale_shape_manual(values = method_shapes[show], breaks = show,
                           labels = method_labels[show], name = NULL, na.value = NA) +
        scale_x_continuous(breaks = seq(0.2, 0.8, by = 0.2),
                           labels = percent_format(accuracy = 1), limits = c(0.05, 0.95)) +
        scale_y_continuous(breaks = seq(-20, 60, by = 20), limits = ylim,
                           oob = scales::squish, labels = percent_format(scale = 1)) +
        # legend keys mirror the plot: references are thick, semi-transparent
        # lines with no symbol; pilot arms are thin lines with their symbol
        guides(colour = guide_legend(
                                     override.aes = list(
                                         shape     = ifelse(show %in% refs, NA, unname(method_shapes[show])),
                                         linetype  = unname(method_linetypes[show]),
                                         linewidth = ifelse(show %in% refs, 1.4, 0.8),
                                         alpha     = ifelse(show %in% refs, 0.75, 1))),
               linetype = "none", shape = "none") +
        labs(x = "Proportion coded", y = "% reduction in empirical variance vs. SRS") +
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom", legend.key.width = unit(1.6, "lines"),
              legend.text = element_text(size = 8.5),
              panel.grid.minor = element_blank(), 
              strip.text = element_text(size = 9), panel.spacing.x = unit(0.7, "lines"))
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
        scale_color_manual(name = "Method:", values = method_pal) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2), 
                           labels = percent_format(accuracy = 1), 
                           limits = c(0.05, 0.95)) +
        scale_y_continuous(limits = c(-0.1, 0.1)) +
        xlab("Proportion coded") + ylab("Monte Carlo Bias") +
        ggtitle(glue::glue("Estimator bias across simulation scenarios{ifelse(!is.null(r2_val), paste0(' (R² target = ', r2_val, ')'), '')}")) +
        theme_minimal() +
        theme(legend.position = "bottom", panel.grid.minor = element_blank())
}




# CI half-width multiplier used to score coverage. 
Z_CRIT <-2

# Default number of paired bootstrap replicates and the seed, so that
# every table in the paper reproduces exactly.
BOOT_B    <- 1000
BOOT_SEED <- 20260916


# Fail loudly if a "cell" is not really one population: the iteration id
# must be unique within (cell x method).  Without this, stacking several
# results files and forgetting a key (e.g. strata_config) silently pools
# different populations into one variance.
check_unique_iterations <- function(raw, cell_keys, where = "") {
    dup <- raw %>%
        count(across(all_of(c(cell_keys, "method", ".n")))) %>%
        filter(n > 1)
    if (nrow(dup)) {
        stop(sprintf(paste0("%s: iteration id `.n` is not unique within cell x method ",
                            "(%d duplicated keys). A cell key is probably missing; ",
                            "keys used: %s"),
                     where, nrow(dup), paste(cell_keys, collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}



# ---------------------------------------------------------------------
# Per-cell summaries of simulation results
# ---------------------------------------------------------------------

#' Summarise raw per-iteration results to one row per (cell x method)
#'
#' @param raw long results as written by 01/02 (one row per iteration x
#'   method), with columns `.n`, `method`, `est`, `SE`, `SEhat`, `tau`.
#' @param cell_keys character vector naming the columns that define a
#'   simulation cell (everything except `method` and `.n`).
#' @return one row per cell x method with the Monte Carlo summaries.
#'
#' Coverage is scored against the cell's own `tau`
cell_summary <- function(raw, cell_keys) {
    stopifnot(all(c("method", "est", "SEhat", "tau") %in% names(raw)))
    raw %>%
        filter(!is.na(est), !is.na(SEhat)) %>%
        mutate(
            CI_lower = est - Z_CRIT * SEhat,
            CI_upper = est + Z_CRIT * SEhat,
            covers   = if_else(CI_lower < tau & CI_upper > tau, 1L, 0L)
        ) %>%
        group_by(across(all_of(c(cell_keys, "method")))) %>%
        summarise(
            n_iter     = n(),
            mean_est   = mean(est),
            bias_est   = mean(est - tau),
            coverage   = mean(covers),
            emp_var    = var(est),
            emp_se     = sd(est),
            mean_se    = mean(SE),
            mean_sehat = mean(SEhat),
            mse        = mean((est - tau)^2),
            .groups    = "drop"
        )
}


#' Aggregate cell-level variances to a variance ratio against a reference
#'
#'
#' @param cells output of cell_summary() (or anything with `method` and
#'   the column named by `value`).
#' @param ref name of the reference method (denominator).
#' @param by character vector of columns to aggregate WITHIN (the ratio
#'   is computed separately for each combination).  Empty = one number.
#' @param value the variance column to aggregate.
#' @return one row per `by` x method, with `mean_var`, `vr` and
#'   `pct_red` (= 100 * (1 - vr)).
agg_vr <- function(cells, ref = "SRS", by = character(0), value = "emp_var") {
    stopifnot(ref %in% unique(as.character(cells$method)))
    num <- cells %>%
        group_by(across(all_of(c(by, "method")))) %>%
        summarise(mean_var = mean(.data[[value]], na.rm = TRUE),
                  n_cells  = sum(!is.na(.data[[value]])),
                  .groups  = "drop")
    den <- num %>%
        filter(as.character(method) == ref) %>%
        select(all_of(by), mean_var_ref = mean_var)
    out <- if (length(by)) left_join(num, den, by = by) else
        mutate(num, mean_var_ref = den$mean_var_ref[1])
    out %>% mutate(vr = mean_var / mean_var_ref,
                   pct_red = 100 * (1 - vr))
}
