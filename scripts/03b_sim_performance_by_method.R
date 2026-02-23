library(tidyverse)
library(broom)
library(patchwork)


method_pal <- c(
    "Subset"      = "#9467bd",
    "SRS"         = "#2ca02c",
    "Stratified"  = "#1f77b4",
    "Optimal"     = "#ff6961",
    "Full coding" = "black"
)


# Load simulation results
results1=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_balanced_exact_strata.rds")
results2=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_balanced_approx_strata.rds")
results3=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_unbalanced_strata.rds")

results=rbind(results1,results2,results3)
results$method[results$method=="Stratified-Opt"]="Optimal"



results_oracle <- readRDS("results/sim_results_oracle.rds")



# Combine and prepare data
dat_sim <- results %>%
    filter(method %in% c("SRS", "Stratified", "Optimal")) %>%
    group_by(method, N, n_code, n_strat, strata_config, tau, bias, 
             residual_var, pop.sd, R2_target) %>%
    summarise(
        emp_SE = sd(est),
        mean_est = mean(est),
        mean_SEhat = mean(SEhat),
        .groups = "drop"
    ) %>%
    mutate(
        log_emp_SE = log(emp_SE),
        coding_frac = n_code / N,
        method = factor(method, levels = c("SRS", "Stratified", "Optimal")),
        bias = factor(bias, levels = c("none", "small", "moderate", "large", "extreme-contrast")),
        residual_var = factor(residual_var, levels = c("homogeneous", "heterogeneous", "extreme-contrast")),
        R2_target = factor(R2_target),
        strata_config = factor(strata_config, levels = c("balanced-exact", "balanced-approx","unbalanced"))
        
    )

# Check data
head(dat_sim)
table(dat_sim$method)


#' Fit regression model for a single method
#' 
#' @param data Data for one method
#' @param method_name Name of the method (for labeling)
fit_method_model <- function(data, method_name,interact=TRUE) {
    
    data = data %>% select(log_emp_SE, coding_frac, bias, residual_var, R2_target, strata_config)
    
    if (interact){
    # Fit model with main effects and key interactions
    model <- lm(log_emp_SE ~  .*., data = data)
    }
    
    else {
        model <- lm(log_emp_SE ~  ., data = data)
        
    }
    
    list(
        model = model,
        method = method_name,
        summary = summary(model),
        tidy = broom::tidy(model),
        glance = broom::glance(model)
    )
}

# Fit models for each method
models <- dat_sim %>%
    group_by(method) %>%
    group_split() %>%
    map(~fit_method_model(.x, unique(.x$method)))

names(models) <- map_chr(models, ~as.character(.x$method))

# Main effects only models
models_main <- dat_sim %>%
    group_by(method) %>%
    group_split() %>%
    map(~fit_method_model(.x, unique(.x$method), interact=F))

names(models_main) <- map_chr(models_main, ~as.character(.x$method))


# Print model summaries
walk(models, function(m) {
    cat("\n\n========================================\n")
    cat("Method:", m$method, "\n")
    cat("========================================\n")
    cat("R-squared:", round(m$glance$r.squared, 3), "\n")
    cat("Adj R-squared:", round(m$glance$adj.r.squared, 3), "\n\n")
    print(m$summary)
})


#' Calculate variable importance metrics
#' 
calc_variable_importance <- function(model_fit) {
    
    # Get model matrix and standardize predictors
    X <- model.matrix(model_fit$model)
    
    # Remove intercept column
    X_scaled <- scale(X[, -1, drop = FALSE])
    
    # Get coefficients and standard errors (excluding intercept)
    coef_summary <- summary(model_fit$model)$coefficients[-1, , drop = FALSE]
    coefs <- coef_summary[, "Estimate"]
    std_errors <- coef_summary[, "Std. Error"]
    
    # Calculate standardized coefficients
    # std_coef = coef * sd(X) / sd(Y)
    sd_y <- sd(model_fit$model$model$log_emp_SE)
    sd_x <- apply(X[, -1, drop = FALSE], 2, sd)
    
    std_coefs <- coefs * sd_x / sd_y
    
    # Calculate standard errors for standardized coefficients
    # SE(std_coef) = SE(coef) * sd(X) / sd(Y)
    std_se <- std_errors * sd_x / sd_y
    
    # Calculate 95% confidence intervals
    ci_lower <- std_coefs - 1.96 * std_se
    ci_upper <- std_coefs + 1.96 * std_se
    
    # Create importance data frame
    importance_df <- data.frame(
        term = names(std_coefs),
        coefficient = coefs,
        std_error = std_errors,
        std_coefficient = std_coefs,
        std_se = std_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        abs_std_coefficient = abs(std_coefs),
        method = model_fit$method
    ) %>%
        arrange(desc(abs_std_coefficient)) %>%
        mutate(
            term_clean = str_replace_all(term, ":", " × "),
            term_clean = str_replace_all(term_clean, "bias", "Bias: "),
            term_clean = str_replace_all(term_clean, "residual_var", "Var: "),
            term_clean = str_replace_all(term_clean, "R2_target", "R² = "),
            term_clean = str_replace_all(term_clean, "coding_frac", "Coding Fraction"),
            term_clean = str_replace_all(term_clean, "strata_configbalanced-approx", "Strata Config: Approx. Balance"),
            term_clean = str_replace_all(term_clean, "strata_configunbalanced", "Strata Config: Unbalanced")
            
            
        )
    
    importance_df
}

# Recalculate importance for each method with standard errors
importance_all <- map_dfr(models, calc_variable_importance)
importance_main <- map_dfr(models_main, calc_variable_importance)


# Print top 10 most important variables for each method
importance_all %>%
    group_by(method) %>%
    slice_max(abs_std_coefficient, n = 10) %>%
    select(method, term_clean, std_coefficient, abs_std_coefficient) %>%
    print(n = 30)


#' Create variable importance plot for one method
plot_importance_single <- function(importance_df, top_n = 15) {
    
    method_name <- unique(importance_df$method)
    
    plot_data <- importance_df %>%
        slice_max(abs_std_coefficient, n = top_n) %>%
        mutate(term_clean = fct_reorder(term_clean, abs_std_coefficient))
    
    ggplot(plot_data, aes(x = abs_std_coefficient, y = term_clean, 
                          fill = std_coefficient > 0)) +
        geom_col() +
        scale_fill_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#b2182b"),
                          labels = c("Decreases SE", "Increases SE"),
                          name = "Effect Direction") +
        labs(
            title = paste("Variable Importance:", method_name),
            subtitle = "Top 15 predictors of log(Empirical SE)",
            x = "|Standardized Coefficient|",
            y = NULL
        ) +
        theme_bw() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 12),
            axis.text.y = element_text(size = 9)
        )
}

# Create importance plots for each method
importance_plots <- importance_all %>%
    group_by(method) %>%
    group_split() %>%
    map(~plot_importance_single(.x, top_n = 15))

# Display plots
walk(importance_plots, print)

# Combined plot
p_combined <- importance_all %>%
    group_by(method) %>%
    slice_max(abs_std_coefficient, n = 15) %>%
    ungroup() %>%
    mutate(term_clean = fct_reorder(term_clean, abs_std_coefficient)) %>%
    ggplot(aes(x = abs_std_coefficient, y = term_clean, fill = std_coefficient > 0)) +
    geom_col() +
    facet_wrap(~method, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#b2182b"),
                      labels = c("Decreases SE", "Increases SE"),
                      name = "Effect Direction") +
    labs(
        title = "Variable Importance by Method",
        subtitle = "Top 15 predictors of log(Empirical SE) for each method",
        x = "|Standardized Coefficient|",
        y = NULL
    ) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 8)
    )

print(p_combined)

# Save plot
ggsave("figures/variable_importance_by_method.pdf", p_combined, 
       width = 8, height = 10)


#' Coefficient plot showing effect direction with confidence intervals
plot_coefficients_comparison <- function(importance_df, top_n = NULL) {
    
    # Get top variables across all methods
    top_vars <- importance_df %>%
        group_by(term_clean) %>%
        summarise(max_importance = max(abs_std_coefficient), .groups = "drop")
    
    if (!is.null(top_n)) {
        top_vars <- top_vars %>%
            slice_max(max_importance, n = top_n)
    }
    
    top_vars <- top_vars %>% pull(term_clean)
    
    # Exponentiate coefficients and CIs to get multiplicative effects
    # exp(coef) gives the multiplicative factor for SE
    # e.g., exp(0.1) = 1.105 means 10.5% increase in SE
    plot_data <- importance_df %>%
        filter(term_clean %in% top_vars) %>%
        mutate(
            term_clean = fct_reorder(term_clean, abs_std_coefficient, .fun = max),
            # Exponentiate to get multiplicative effect
            mult_effect = exp(std_coefficient),
            mult_ci_lower = exp(ci_lower),
            mult_ci_upper = exp(ci_upper)
        )
    
    # Create custom breaks for log scale (as percent change from 1.0)
    # Show breaks at meaningful percent changes
    #breaks <- c(0.5, 0.67, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0)
    breaks <- seq(0.25,2.0,by=0.25)# c(0.25, 0.5, 0.67, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0)
    
    labels <- scales::percent(breaks - 1, accuracy = 1)
    labels[breaks == 1.0] <- "0%"
    
    ggplot(plot_data, aes(x = mult_effect, y = term_clean, color = method)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
        geom_errorbarh(aes(xmin = mult_ci_lower, xmax = mult_ci_upper),
                       height = 0.3, alpha = 0.5, linewidth = 0.8,
                       position = position_dodge(width = 0.5)) +
        geom_point(size = 3, alpha = 0.8,
                   position = position_dodge(width = 0.5)) +
        scale_x_continuous(trans = "log",
                           breaks = breaks,
                           labels = labels) +
        scale_color_manual(name = "Method", values = method_pal) +
        labs(
            title = "Effect on Standard Error by Method",
            subtitle = "Top predictors of empirical SE across all methods (with 95% CI)",
            x = "% Change in SE (log scale)",
            y = NULL,
            caption = "Negative = Decreases SE (better), Positive = Increases SE (worse)\nError bars show 95% confidence intervals"
        ) +
        theme_bw() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 12),
            axis.text.y = element_text(size = 9)
        )
}

p_comparison <- plot_coefficients_comparison(importance_all,top_n=20)
print(p_comparison)

ggsave("figures/coefficient_comparison_with_CI.pdf", 
       p_comparison, width = 8, height = 8)


p_comparison_main <- plot_coefficients_comparison(importance_main)
p_comparison_main


ggsave("figures/coefficient_comparison_with_CI_main_effects_only.pdf", 
       p_comparison_main, width = 8, height = 8)


#' Alternative version: separate panels for clearer visualization
plot_coefficients_comparison_faceted <- function(importance_df,top_n=NULL) {
    
    # Get top variables across all methods
    top_vars <- importance_df %>%
        group_by(term_clean) %>%
        summarise(max_importance = max(abs_std_coefficient), .groups = "drop") 
    
    if (!is.null(top_n)){
        top_vars = top_vars %>% 
            slice_max(max_importance, n = top_n)
    }
    top_vars = top_vars %>% pull(term_clean)

    
    plot_data <- importance_df %>%
        filter(term_clean %in% top_vars) %>%
        mutate(term_clean = fct_reorder(term_clean, abs_std_coefficient, .fun = max),
               sig = ifelse(ci_lower * ci_upper > 0, "Significant", "Not Significant"))
    
    ggplot(plot_data, aes(x = std_coefficient, y = term_clean, color = sig)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                       height = 0.3, alpha = 0.5, linewidth = 0.8) +
        geom_point(size = 3, alpha = 0.8) +
        facet_wrap(~method, ncol = 1, scales = "free_x") +
        scale_color_manual(values = c("Significant" = "#d7191c", "Not Significant" = "grey60"),
                           name = "Statistical Significance") +
        labs(
            title = "Standardized Coefficients by Method",
            subtitle = "Top predictors (with 95% CI)",
            x = "Standardized Coefficient",
            y = NULL,
            caption = "Negative = Decreases SE (better), Positive = Increases SE (worse)"
        ) +
        theme_bw() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 12),
            strip.background = element_rect(fill = "grey90"),
            strip.text = element_text(face = "bold", size = 11),
            axis.text.y = element_text(size = 9)
        )
}

p_comparison_faceted <- plot_coefficients_comparison_faceted(importance_all,top_n=10)
print(p_comparison_faceted)

ggsave("figures/coefficient_comparison_faceted.pdf", p_comparison_faceted,
       width = 10, height = 12)


#' Create detailed coefficient table with CIs
create_coefficient_table <- function(model_results, top_n = 50) {
    
    ref_importance <- calc_variable_importance(model_results$SRS)
    top_terms <- ref_importance %>% 
        head(top_n) %>% 
        pull(term)
    
    coef_table <- map_dfr(model_results, function(res) {
        res$tidy %>%
            filter(term %in% top_terms) %>%
            mutate(
                method = res$method,
                ci_lower = estimate - 1.96 * std.error,
                ci_upper = estimate + 1.96 * std.error,
                ci = sprintf("[%.3f, %.3f]", ci_lower, ci_upper),
                est_coefficient = sprintf("%.3f", estimate),
                p_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
            ) %>%
            select(method, term, est_coefficient, ci, p_val)
    }) %>%
        pivot_wider(
            names_from = method,
            values_from = c(est_coefficient, ci, p_val),
            names_glue = "{method}_{.value}"
        ) %>%
        arrange(term)
    
    return(coef_table)
}

coef_table <- create_coefficient_table(models, top_n = 50)


# Export to CSV
write_csv(coef_table, "results/coefficient_table_with_CI.csv")

# Print formatted table
library(knitr)
library(kableExtra)

coef_table %>%
    kable(format = "latex", booktabs = TRUE, 
          caption = "Coefficient Estimates with 95% Confidence Intervals") %>%
    kable_styling(latex_options = c("scale_down", "hold_position")) %>%
    cat(file = "tables/coefficient_comparison_table.tex")

print(coef_table)
