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
results_oracle <- readRDS("results/sim_results_oracle.rds")

results1=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_balanced_exact_strata.rds")
results2=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_balanced_approx_strata.rds")
results3=readRDS("~/Dropbox/text and RCT sampling math notes/stratsampling/results/sim_results_unbalanced_strata.rds")

results=rbind(results1,results2,results3)
results$method[results$method=="Stratified-Opt"]="Optimal"



# =============================================================================
# Simplified model with binary bias and heterogeneity indicators
# =============================================================================

# Prepare data with binary recoding
dat_sim_binary <- results %>%
    filter(method %in% c("SRS", "Stratified", "Optimal"), bias!="extreme-contrast") %>%
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
        
        # Binary recoding: any_bias = 1 if bias is not "none"
        any_bias = factor(ifelse(bias == "none", 0, 1), 
                          levels = c(0, 1), 
                          labels = c("No Bias", "Bias")),
        bias2 = factor(ifelse(bias=="extreme-contrast","large",bias),levels=c("none","small", "moderate","large")),
        
        # Binary recoding: any_heterogeneity = 1 if residual_var is not "homogeneous"
        any_heterogeneity = factor(ifelse(residual_var == "homogeneous", 0, 1),
                                   levels = c(0, 1),
                                   labels = c("Homogeneous Var.", "Heterogeneous Var.")),
        
        R2_target = factor(R2_target),
        strata_config = factor(strata_config, levels = c("balanced-exact", "balanced-approx","unbalanced"))
        
    )



# Verify recoding
table(dat_sim_binary$any_bias)
table(dat_sim_binary$any_heterogeneity)

# Fit model with binary predictors
fit_method_model_binary <- function(data, method_name,interact=TRUE) {
    
    
    data = data %>% select(log_emp_SE, coding_frac, strata_config, any_bias, any_heterogeneity,R2_target)
    
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

# Fit models for each method with binary predictors
models_binary <- dat_sim_binary %>%
    group_by(method) %>%
    group_split() %>%
    map(~fit_method_model_binary(.x, unique(.x$method)))

names(models_binary) <- map_chr(models_binary, ~as.character(.x$method))

models_binary_main <- dat_sim_binary %>%
    group_by(method) %>%
    group_split() %>%
    map(~fit_method_model_binary(.x, unique(.x$method),interact=F))

names(models_binary_main) <- map_chr(models_binary_main, ~as.character(.x$method))


# Print model summaries
walk(models_binary, function(m) {
    cat("\n\n========================================\n")
    cat("Method:", m$method, "(Binary Predictors)\n")
    cat("========================================\n")
    cat("R-squared:", round(m$glance$r.squared, 3), "\n")
    cat("Adj R-squared:", round(m$glance$adj.r.squared, 3), "\n\n")
    print(m$summary)
})

# Calculate variable importance for binary models
calc_variable_importance_binary <- function(model_fit) {
    
    X <- model.matrix(model_fit$model)
    X_scaled <- scale(X[, -1, drop = FALSE])
    
    coef_summary <- summary(model_fit$model)$coefficients[-1, , drop = FALSE]
    coefs <- coef_summary[, "Estimate"]
    std_errors <- coef_summary[, "Std. Error"]
    
    sd_y <- sd(model_fit$model$model$log_emp_SE)
    sd_x <- apply(X[, -1, drop = FALSE], 2, sd)
    
    std_coefs <- coefs * sd_x / sd_y
    std_se <- std_errors * sd_x / sd_y
    
    ci_lower <- std_coefs - 1.96 * std_se
    ci_upper <- std_coefs + 1.96 * std_se
    
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
            term_clean = str_replace_all(term_clean, "any_biasBias", "Bias"),
            term_clean = str_replace_all(term_clean, "bias2small", "Bias (Small)"),
            term_clean = str_replace_all(term_clean, "bias2moderate", "Bias (Moderate)"),
            term_clean = str_replace_all(term_clean, "bias2large", "Bias (Large)"),
            term_clean = str_replace_all(term_clean, "any_heterogeneityHeterogeneous Var.", "Heterogeneous Var."),
            term_clean = str_replace_all(term_clean, "R2_target", "R² = "),
            term_clean = str_replace_all(term_clean, "coding_frac", "Coding Fraction"),
            term_clean = str_replace_all(term_clean, "strata_configbalanced-approx", "Strata Config: Approx. Balance"),
            term_clean = str_replace_all(term_clean, "strata_configunbalanced", "Strata Config: Unbalanced")
            
        )
    importance_df
}

# Calculate importance for all methods
importance_all_binary <- map_dfr(models_binary, calc_variable_importance_binary)
importance_all_binary_main <- map_dfr(models_binary_main, calc_variable_importance_binary)

# Print top predictors
importance_all_binary %>%
    group_by(method) %>%
    slice_max(abs_std_coefficient, n = 10) %>%
    select(method, term_clean, std_coefficient, abs_std_coefficient) %>%
    print(n = 30)

# Create coefficient comparison plot with binary predictors
plot_coefficients_comparison_binary <- function(importance_df,top_n=NULL) {
    
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
        mutate(term_clean = fct_reorder(term_clean, abs_std_coefficient, .fun = max))
    
    ggplot(plot_data, aes(x = std_coefficient, y = term_clean, color = method)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                       height = 0.3, alpha = 0.5, linewidth = 0.8,
                       position = position_dodge(width = 0.5)) +
        geom_point(size = 3, alpha = 0.8,
                   position = position_dodge(width = 0.5)) +
        scale_color_manual(name = "Method", values = method_pal) +
        labs(
            title = "Standardized Coefficients by Method (Binary Predictors)",
            subtitle = "Predictors of log(empirical SE) with 95% CI for each method",
            
            x = "Standardized Coefficient",
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

p_comparison_binary <- plot_coefficients_comparison_binary(importance_all_binary)
print(p_comparison_binary)

# Save the plot
ggsave("figures/coefficient_comparison_with_CI_binary.pdf", 
       p_comparison_binary, width = 8, height = 8)




# =============================================================================
# Simplified model with bias (large and extreme contrast collapsed) and heterogeneity indicators
# =============================================================================


# Fit model with binary predictors (preserving bias moderate vs. large; excluding bias:extreme-contrast)
fit_method_model_binary2 <- function(data, method_name, interact=TRUE) {
    
    data = data %>% select(log_emp_SE, coding_frac, strata_config, bias2, any_heterogeneity,R2_target)
    
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

# Fit models for each method with binary predictors
models_binary2 <- dat_sim_binary %>%
    group_by(method) %>%
    group_split() %>%
    map(~fit_method_model_binary2(.x, unique(.x$method)))

names(models_binary2) <- map_chr(models_binary2, ~as.character(.x$method))


# Calculate importance for all methods
importance_all_binary2 <- map_dfr(models_binary2, calc_variable_importance_binary)



p_comparison_binary2 <- plot_coefficients_comparison_binary(importance_all_binary2)
print(p_comparison_binary2)

# Save the plot
ggsave("figures/coefficient_comparison_with_CI_binary_v2.pdf", 
       p_comparison_binary, width = 8, height = 8)

