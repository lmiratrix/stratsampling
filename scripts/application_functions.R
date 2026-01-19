
# Function to calculate oracle variance (i.e., true variance of the estimator based on variance of residuals for all observations) for a given method and coding budget
calc_oracle_variance = function(dat, method = c("SRS", "Stratified-Prop", "Stratified-Opt"),
                                 h = 0.33) {
    method = match.arg(method)
    N = nrow(dat)
    n_total = round(h * N)
    
    # Experimental variance component (same for all methods)
    var_exp = sum(oracle_result$var_Y / oracle_result$N)
    
    if (method == "SRS") {
        # Simple random sampling variance component
        var_samp = dat %>%
            group_by(Z) %>%
            summarise(
                N_z = n(),
                n_z = round(n_total / 2),
                var_resid = var(resid),
                component = ((N_z - n_z) / N_z) * (var_resid / n_z),
                .groups = "drop"
            ) %>%
            pull(component) %>%
            sum()
        
    } else if (method == "Stratified-Prop") {
        # Proportional allocation
        var_samp = oracle_strata_stats %>%
            group_by(Z) %>%
            mutate(
                N_z = sum(N_zk),
                n_z = round(n_total / 2),
                n_zk = round(n_z * N_zk / N_z),
                component = (N_zk * (N_zk - n_zk)) / (N_z^2 * n_zk) * var_resid
            ) %>%
            pull(component) %>%
            sum()
        
    } else {  # Stratified-Opt
        # Optimal (Neyman) allocation
        var_samp = oracle_strata_stats %>%
            group_by(Z) %>%
            mutate(
                N_z = sum(N_zk),
                n_z = round(n_total / 2),
                weight = N_zk * sqrt(var_resid),
                n_zk = pmax(1, round(n_z * weight / sum(weight))),
                n_zk = pmin(n_zk, N_zk),  # Cap at stratum size
                component = (N_zk * (N_zk - n_zk)) / (N_z^2 * n_zk) * var_resid
            ) %>%
            pull(component) %>%
            sum()
    }
    
    total_var = var_exp + var_samp
    
    list(
        method = method,
        h = h,
        var_exp = var_exp,
        var_samp = var_samp,
        var_total = total_var,
        SE = sqrt(total_var)
    )
}

# -----------------------------------------------------------------------------
# Function to compute estimates from a single sample
# -----------------------------------------------------------------------------

#' Compute treatment effect estimates from a single stratified sample
#' 
#' @param dat Data frame with full data (ID, Yobs, Yhat, strata)
#' @param Z treatment indicator
#' @param coding_fraction Fraction of data to code
#' @param finite_pop Should finite population correction be applied?
#' @param seed Random seed for reproducibility
#' @return Data frame with estimates from all methods plus oracle
compute_single_sample_estimates = function(dat, Z="Z", coding_fraction = 0.30, finite_pop=TRUE, seed = NULL) {
    
    if (!is.null(seed)) set.seed(seed)

    N = nrow(dat)
    n_code = round(coding_fraction * N / 2) * 2  # Ensure even split
    
    # Run one iteration of all methods using the shared simulation function
    results = one_iter_all(dat, n_code = n_code, Z = Z, finite_pop = finite_pop)
    
    # Add oracle for comparison
    oracle_est = est_oracle(dat,Z=Z)
    
    # Combine results
    all_results = bind_rows(
        results %>% select(method, est, SE = SEhat),
        oracle_est
    ) %>%
        mutate(
            ci_lower = est - 1.96 * SE,
            ci_upper = est + 1.96 * SE,
            n_coded = ifelse(method == "Oracle", N, n_code),
            coding_pct = ifelse(method == "Oracle", 1.0, coding_fraction)
        )
    
    return(all_results)
}


#' Run multiple samples to show variability in estimates
#' 
#' @param dat Data frame with full data
#' @param Z treatment indicator
#' @param coding_fraction Fraction of data to code
#' @param n_samples Number of independent samples to draw
#' @param finite_pop Should finite population correction be applied?
#' @param seed Base random seed
#' @return Data frame with estimates from all samples
compute_multiple_sample_estimates = function(dat, Z="Z", coding_fraction = 0.30, 
                                              n_samples = 10, finite_pop=TRUE, seed=NULL) {
    
    if (!is.null(seed)) set.seed(seed)
    
    results = purrr::map_dfr(1:n_samples, function(i) {
        compute_single_sample_estimates(dat, Z=Z, coding_fraction=coding_fraction, finite_pop=finite_pop) %>%
            mutate(sample_id = i)
    })
    
    return(results)
}


# -----------------------------------------------------------------------------
# Create forest plot visualization
# -----------------------------------------------------------------------------

#' Create forest plot comparing estimates across methods
#' 
#' @param estimates Data frame from compute_single_sample_estimates()
#' @param oracle_est Oracle point estimate for reference line
#' @param tx Is the estimate a treatment effect or mean?
#' @return ggplot object
create_forest_plot = function(estimates, oracle_est = NULL, tx=TRUE) {
    
    xlab_name = ifelse(tx, "Estimated Treatment Effect (Effect Size Units)", "Estimated Mean")
    
    # Order methods logically
    method_order = c("Oracle", "Subset", "SRS", "Stratified", "Stratified-Opt")
    estimates = estimates %>%
        mutate(method = factor(method, levels = rev(method_order)))
    
    p = ggplot(estimates, aes(x = est, y = method)) +
        geom_point(size = 3, aes(color = method)) +
        geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, linewidth = 0.8) +
        labs(
            x = xlab_name,
            y = "",
            title = "Estimates by Method",
            subtitle = paste0("Coding budget: ", scales::percent(estimates$coding_pct[estimates$method != "Oracle"][1]))
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.text.y = element_text(size = 11),
            panel.grid.major.y = element_blank()
        )
    
    if (tx){
        p =   p + geom_vline(xintercept = 0, linetype = "dashed", color = "gray60")
    }
    else if (!tx & !is.null(oracle_est)){
        xl=round(oracle_est-5*max(estimates$SE),1)
        xu=round(oracle_est+5*max(estimates$SE),1)
        p = p + scale_x_continuous(limits=c(xl,xu))
    }
    # Add oracle reference line if provided
    if (!is.null(oracle_est)) {
        p = p + geom_vline(xintercept = oracle_est, linetype = "solid", 
                            color = "darkgreen", linewidth = 0.8, alpha = 0.7)
    }
    
    return(p)
}


#' Create forest plot showing multiple samples
#' 
#' @param multi_estimates Data frame from compute_multiple_sample_estimates()
#' @param oracle_est Oracle point estimate for reference line
#' @param tx Treatment effect or mean?
#' @param jitter_seed Optional seed for the jitter
#' @return ggplot object
create_multi_sample_forest_plot = function(multi_estimates, oracle_est = NULL, tx=TRUE, jitter_seed=NA) {
    
    xlab_name = ifelse(tx, "Estimated Treatment Effect (Effect Size Units)", "Estimated Mean")
    
    
    # Order methods logically
    method_order = c("Oracle", "Subset", "SRS", "Stratified", "Stratified-Opt")
    
    # Separate oracle (only need one) from other methods
    oracle_data = multi_estimates %>%
        filter(method == "Oracle") %>%
        slice(1) %>%
        mutate(method = factor(method, levels = rev(method_order)))
    
    other_data = multi_estimates %>%
        filter(method != "Oracle") %>%
        mutate(
            method = factor(method, levels = rev(method_order)),
            sample_id = factor(sample_id)
        )
    
    
    avg_data = other_data %>%
        group_by(method) %>%
        summarise(
            mean_est = mean(est),
            mean_SE = mean(SE),
            .groups = "drop"
        ) %>%
        mutate(
            ci_lower = mean_est - 1.96 * mean_SE,
            ci_upper = mean_est + 1.96 * mean_SE
        )
    
    # Combine with oracle
    plot_data = bind_rows(
        oracle_data %>% mutate(sample_id = factor(1)),
        other_data
    )
    
    p = ggplot(plot_data, aes(x = est, y = method)) +
        # Show average estimate with error bars for non-oracle methods
        geom_errorbarh(data = avg_data,
                       aes(x = mean_est, xmin = ci_lower, xmax = ci_upper),
                       height = 0.15, linewidth = 0.8, color = "black",alpha=0.6) +
        geom_point(data = avg_data,
                   aes(x = mean_est),
                   size = 4.5, color = "black", alpha=0.75) +
        # Show individual sample estimates with some jitter
        geom_point(data = other_data, 
                   aes(color = sample_id), 
                   position = position_jitter(height = 0.18, width = 0,seed=jitter_seed),
                   size = 3, alpha = 0.85) +
        # Show oracle as a distinct point
        geom_point(data = oracle_data, size = 4.5, color = "black") +
        geom_errorbarh(data = oracle_data, 
                       aes(xmin = ci_lower, xmax = ci_upper), 
                       height = 0.15, linewidth = 0.8, color = "black") +
        labs(
            x = xlab_name,
            y = "",
            title = "Treatment Effect Estimates Across Multiple Samples",
            subtitle = paste0("Each point = one random sample at ", 
                              scales::percent(other_data$coding_pct[1]), " coding budget"),
            color = "Sample"
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.text.y = element_text(size = 11),
            panel.grid.major.y = element_blank()
        )
    
    if (tx){
        p =   p + geom_vline(xintercept = 0, linetype = "dashed", color = "gray60")
    }
    else if (!tx & !is.null(oracle_est)){
        xl=round(oracle_est-5*max(multi_estimates$SE),1)
        xu=round(oracle_est+5*max(multi_estimates$SE),1)
        p = p + scale_x_continuous(limits=c(xl,xu))
    }
    
    # Add oracle reference line
    if (!is.null(oracle_est)) {
        p = p + geom_vline(xintercept = oracle_est, linetype = "solid", 
                            color = "darkgreen", linewidth = 0.5, alpha = 0.5)
    }
    
    return(p)
}


# -----------------------------------------------------------------------------
# Create comparison table
# -----------------------------------------------------------------------------

#' Create formatted table comparing methods
#' 
#' @param estimates Data frame from compute_single_sample_estimates()
#' @return Formatted data frame
create_comparison_table = function(estimates) {
    
    method_order = c("Oracle", "Subset", "SRS", "Stratified", "Stratified-Opt")
    
    oracle_se = estimates$SE[estimates$method == "Oracle"]
    
    estimates %>%
        mutate(
            method = factor(method, levels = method_order),
            # Variance inflation relative to oracle
            var_inflation = (SE / oracle_se)^2,
            # Significance
            significant = ci_lower > 0 | ci_upper < 0,
            sig_marker = ifelse(significant, "*", "")
        ) %>%
        arrange(method) %>%
        select(
            Method = method,
            `n Coded` = n_coded,
            Estimate = est,
            SE = SE,
            `95% CI Lower` = ci_lower,
            `95% CI Upper` = ci_upper,
            `Var. Inflation` = var_inflation,
            `Sig.` = sig_marker
        ) %>%
        mutate(across(where(is.numeric), ~round(., 3)))
}


#' Create summary table across multiple samples
#' 
#' @param multi_estimates Data frame from compute_multiple_sample_estimates()
#' @return Summary data frame
create_multi_sample_summary = function(multi_estimates) {
    
    method_order = c("Oracle", "Subset", "SRS", "Stratified", "Stratified-Opt")
    
    oracle_data = multi_estimates %>%
        filter(method == "Oracle") %>%
        slice(1)
    
    multi_estimates %>%
        group_by(method) %>%
        summarise(
            n_samples = n(),
            mean_est = mean(est),
            sd_est = sd(est),
            min_est = min(est),
            max_est = max(est),
            mean_SE = mean(SE),
            pct_significant = mean(ci_lower > 0 | ci_upper < 0) * 100,
            .groups = "drop"
        ) %>%
        mutate(
            var_inflation = sprintf("%.2f",(mean_SE / oracle_data$SE)^2),
            method = factor(method, levels = method_order),
            mean_est=sprintf("%.3f", mean_est),
            sd_est=sprintf("%.3f", sd_est),
            mean_SE=sprintf("%.3f", mean_SE),
            range=sprintf("[%.3f, %.3f]", min_est, max_est)
        ) %>%
        arrange(method) %>%
        select(
            Method = method,
            `Mean Est.` = mean_est,
            `SD of Est.` = sd_est,
            `Range` = range,
            `Mean SE` = mean_SE,
            `Var. Inflation` = var_inflation,
            `% Significant` = pct_significant
        ) %>%
        mutate(across(where(is.numeric), ~round(., 3)))
}

