#' Create Candidate Stratification Variables
#'
#' Generates a set of candidate stratification variables for use in stratified
#' model-assisted estimation. Creates single-variable stratifications (tertiles,
#' quartiles, quintiles) for each specified variable, as well as interaction
#' stratifications for all pairs of variables.
#'
#' @param dat Data frame containing the variables to stratify on.
#' @param vars Character vector specifying the names of columns in \code{dat}
#'   to use for creating stratifications.
#' @param Z Optional character string specifying the name of the treatment
#'   indicator variable in \code{dat}. If provided, stratifications are computed
#'   separately within each treatment arm, ensuring balanced strata across arms.
#'   Default is \code{NULL} (pooled stratification).
#' @param jitter_vars Character vector specifying which variables (subset of
#'   \code{vars}) should have small random noise added before computing quantile
#'   breaks. This helps avoid ties that can cause \code{cut()} to fail. Default
#'   is \code{NULL} (no jitter applied).
#' @param sd_jitter Numeric value specifying the standard deviation of the
#'   jitter noise. Default is \code{1e-3}.
#' @param seed Optional integer seed for reproducibility when jitter is applied.
#'   Default is \code{NULL}.
#'
#' @return A named list of factor vectors, each representing a candidate
#'   stratification. Names follow the convention:
#'   \describe{
#'     \item{Single-variable}{`{var}_tertiles`, `{var}_quartiles`, `{var}_quintiles`}
#'     \item{Interactions}{`{var1}_{var2}_2x2`, `{var1}_{var2}_2x3`, `{var1}_{var2}_3x2`, `{var1}_{var2}_3x3`}
#'   }
#'   Factor levels have meaningful labels (e.g., "Low Yhat", "Med Yhat / High word_count").
#'   Stratifications that fail due to non-unique quantile breaks are automatically
#'   excluded from the output.
#'
#' @details
#' For each variable in \code{vars}, the function creates three stratifications
#' based on quantiles:
#' \itemize{
#'   \item Tertiles (3 groups): Low, Med, High
#'   \item Quartiles (4 groups): Low, Med-Low, Med-High, High
#'   \item Quintiles (5 groups): Very Low, Low, Med, High, Very High
#' }
#'
#' For each unique pair of variables, the function creates four interaction
#' stratifications by crossing quantile groups:
#' \itemize{
#'   \item 2x2: median split on both variables (4 strata)
#'   \item 2x3: median split on first, tertiles on second (6 strata)
#'   \item 3x2: tertiles on first, median split on second (6 strata)
#'   \item 3x3: tertiles on both variables (9 strata)
#' }
#'
#' When \code{Z} is provided, quantile breaks are computed separately within
#' each treatment arm. This ensures that the same proportion of treated and
#' control units fall into each stratum level, which is important for
#' maintaining balance in stratified estimation.
#'
#' When a variable has many tied values, quantile breaks may not be unique,
#' causing \code{cut()} to fail. Use \code{jitter_vars} to add small random
#' noise (SD = 0.001) to such variables before computing breaks.
#'
#' @export
create_stratification_candidates <- function(dat, vars, Z = NULL, jitter_vars = NULL, sd_jitter = 1e-3, seed = NULL) {
    
    if (!is.null(seed)) set.seed(seed)
    
    n <- nrow(dat)
    
    # -------------------------------------------------------------------------
    # Helper functions
    # -------------------------------------------------------------------------
    
    # Labels for different numbers of groups
    get_labels <- function(n_groups, var_name) {
        base_labels <- switch(as.character(n_groups),
                              "2" = c("Low", "High"),
                              "3" = c("Low", "Med", "High"),
                              "4" = c("Low", "Med-Low", "Med-High", "High"),
                              "5" = c("Very Low", "Low", "Med", "High", "Very High")
        )
        paste(base_labels, var_name)
    }
    
    # Create quantile-based strata (returns NULL if breaks aren't unique)
    make_quantile_strata <- function(x, n_groups, var_name) {
        probs <- seq(0, 1, length.out = n_groups + 1)
        breaks <- quantile(x, probs = probs, na.rm = TRUE)
        if (dplyr::n_distinct(breaks) < n_groups + 1) return(NULL)
        labels <- get_labels(n_groups, var_name)
        cut(x, breaks = breaks, labels = labels, include.lowest = TRUE)
    }
    
    # Create quantile-based strata separately by treatment arm
    make_quantile_strata_by_arm <- function(x, z_vec, n_groups, var_name) {
        probs <- seq(0, 1, length.out = n_groups + 1)
        labels <- get_labels(n_groups, var_name)
        
        # Initialize result vector
        result <- factor(rep(NA_character_, length(x)), levels = labels)
        
        # Process each treatment arm separately
        for (z_val in unique(z_vec)) {
            idx <- which(z_vec == z_val)
            x_arm <- x[idx]
            
            breaks <- quantile(x_arm, probs = probs, na.rm = TRUE)
            if (dplyr::n_distinct(breaks) < n_groups + 1) return(NULL)
            
            result[idx] <- cut(x_arm, breaks = breaks, labels = labels, include.lowest = TRUE)
        }
        
        result
    }
    
    # Create interaction strata
    make_interaction_strata <- function(x1, n1, var_name1, x2, n2, var_name2, z_vec = NULL) {
        if (is.null(z_vec)) {
            s1 <- make_quantile_strata(x1, n1, var_name1)
            s2 <- make_quantile_strata(x2, n2, var_name2)
        } else {
            s1 <- make_quantile_strata_by_arm(x1, z_vec, n1, var_name1)
            s2 <- make_quantile_strata_by_arm(x2, z_vec, n2, var_name2)
        }
        if (is.null(s1) || is.null(s2)) return(NULL)
        
        # Create interaction with meaningful labels
        int_factor <- interaction(s1, s2, drop = TRUE, sep = " / ")
        int_factor
    }
    
    # Apply jitter to a variable
    add_jitter <- function(x) x + rnorm(length(x), 0, sd_jitter)
    
    # -------------------------------------------------------------------------
    # Validate and extract treatment indicator if provided
    # -------------------------------------------------------------------------
    
    z_vec <- NULL
    if (!is.null(Z)) {
        if (!Z %in% names(dat)) {
            stop(paste("Treatment variable", Z, "not found in dat"))
        }
        z_vec <- dat[[Z]]
        if (length(unique(z_vec)) < 2) {
            warning("Treatment variable Z has fewer than 2 unique values. Ignoring Z.")
            z_vec <- NULL
        }
    }
    
    # -------------------------------------------------------------------------
    # Prepare variables (with jitter for specified vars)
    # -------------------------------------------------------------------------
    
    if (is.null(jitter_vars)) jitter_vars <- character(0)
    
    var_data <- lapply(vars, function(v) {
        x <- dat[[v]]
        if (is.null(x)) stop(paste("Variable", v, "not found in dat"))
        if (v %in% jitter_vars) add_jitter(x) else x
    })
    names(var_data) <- vars
    
    # -------------------------------------------------------------------------
    # Generate single-variable stratifications (tertiles, quartiles, quintiles)
    # -------------------------------------------------------------------------
    
    single_strats <- list()
    
    for (v in vars) {
        if (is.null(z_vec)) {
            # Pooled stratification
            single_strats[[paste0(v, "_tertiles")]]  <- make_quantile_strata(var_data[[v]], 3, v)
            single_strats[[paste0(v, "_quartiles")]] <- make_quantile_strata(var_data[[v]], 4, v)
            single_strats[[paste0(v, "_quintiles")]] <- make_quantile_strata(var_data[[v]], 5, v)
        } else {
            # Stratification by treatment arm
            single_strats[[paste0(v, "_tertiles")]]  <- make_quantile_strata_by_arm(var_data[[v]], z_vec, 3, v)
            single_strats[[paste0(v, "_quartiles")]] <- make_quantile_strata_by_arm(var_data[[v]], z_vec, 4, v)
            single_strats[[paste0(v, "_quintiles")]] <- make_quantile_strata_by_arm(var_data[[v]], z_vec, 5, v)
        }
    }
    
    # -------------------------------------------------------------------------
    # Generate interaction stratifications for all pairs of variables
    # -------------------------------------------------------------------------
    
    interaction_strats <- list()
    
    if (length(vars) >= 2) {
        # All unique pairs
        var_pairs <- combn(vars, 2, simplify = FALSE)
        
        # Interaction configurations: (n_groups_var1, n_groups_var2)
        configs <- list(
            "2x2" = c(2, 2),
            "2x3" = c(2, 3),
            "3x2" = c(3, 2),
            "3x3" = c(3, 3)
        )
        
        for (pair in var_pairs) {
            v1 <- pair[1]
            v2 <- pair[2]
            
            for (config_name in names(configs)) {
                n1 <- configs[[config_name]][1]
                n2 <- configs[[config_name]][2]
                
                strat_name <- paste0(v1, "_", v2, "_", config_name)
                interaction_strats[[strat_name]] <- make_interaction_strata(
                    var_data[[v1]], n1, v1,
                    var_data[[v2]], n2, v2,
                    z_vec = z_vec
                )
            }
        }
    }

    candidates <- c(single_strats, interaction_strats)
    candidates <- Filter(Negate(is.null), candidates)
    
    return(candidates)
}


#' Evaluate a Candidate Stratification Variable
#'
#' Computes stratum-level residual statistics and heterogeneity metrics to 
#' assess the potential benefit of a stratification for model-assisted estimation.
#'
#' @param data Data frame containing a \code{resid} column (residuals = Yobs - Yhat)
#'   and the stratification variable specified by \code{strata}.
#' @param strata_var Character string specifying the name of the column to use as the stratification variable.
#' @param Z Character string specifying the name of the treatment indicator  or \code{NULL} to ignore treatment groups. 
#' @param min_var_threshold Numeric. Minimum within-stratum residual variance 
#'   required. If any stratum has variance below this threshold, the function 
#'   returns \code{NULL}. This guards against strata with near-zero variance 
#'   that could cause numerical issues. Default is \code{0.01}.
#'
#' @return A list with three elements, or \code{NULL} if the stratification is 
#'   invalid (e.g., empty data, missing values, or variance below threshold):
#'   \describe{
#'     \item{\code{strata_stats}}{Data frame with stratum-level statistics:
#'       \itemize{
#'         \item \code{N_zk} (or \code{N_k}): Sample size in each stratum
#'         \item \code{mean_resid}: Mean residual in each stratum
#'         \item \code{var_resid}: Variance of residuals in each stratum
#'         \item \code{sd_resid}: Standard deviation of residuals in each stratum
#'       }
#'       If \code{Z} is specified, includes treatment indicator \code{Z}.
#'     }
#'     \item{\code{heterogeneity}}{Data frame summarizing heterogeneity across strata:
#'       \itemize{
#'         \item \code{n_strata}: Number of strata
#'         \item \code{min_size}, \code{max_size}: Range of stratum sizes
#'         \item \code{imbalance}: Ratio of largest to smallest stratum size
#'         \item \code{min_var}, \code{max_var}: Range of residual variances
#'         \item \code{ratio_var}: Ratio of max to min residual variance (heteroskedasticity measure)
#'         \item \code{range_var}: Difference between max and min residual variance
#'         \item \code{var_mean_resid}: Variance of stratum mean residuals (proxy for BS)
#'       }
#'       If \code{Z} is specified, returns one row per treatment arm.
#'     }
#'     \item{\code{stratification}}{Character string echoing the stratification variable name.}
#'   }
#'
#' @details
#' This function is useful for comparing candidate stratification variables 
#' before selecting one for stratified sampling. Key metrics to examine:
#' 
#' \itemize{
#'   \item \code{var_mean_resid}: Higher values indicate greater between-strata 
#'     variance (BS), suggesting stratification will reduce sampling variance.
#'   \item \code{ratio_var}: Values greater than 2 indicate heteroskedasticity 
#'     across strata, suggesting optimal (Neyman) allocation will outperform 
#'     proportional allocation.
#'   \item \code{imbalance}: Values close to 1 indicate balanced strata; highly 

#'     unbalanced strata may reduce efficiency.
#' }
#' @export
evaluate_stratification <- function(data, strata_var, Z = "Z", 
                                    min_var_threshold = 0.01) {
    
    data$strata_eval <- data[[strata_var]]
    data <- data %>% filter(!is.na(strata_eval))
    if (nrow(dat) == 0) return(NULL)
    
    if (!is.null(Z)) {
        data$Z <- data[[Z]]
        
        strata_stats <- data %>%
            group_by(Z, strata_eval) %>%
            summarise(
                N_zk = n(),
                mean_resid = mean(resid),
                var_resid = var(resid),
                sd_resid = sd(resid),
                .groups = "drop"
            ) %>%
            filter(!is.na(var_resid))
        
        if (nrow(strata_stats) == 0) return(NULL)
        
        min_var <- min(strata_stats$var_resid, na.rm = TRUE)
        if (is.na(min_var) || min_var < min_var_threshold) return(NULL)
        
        heterogeneity <- strata_stats %>%
            group_by(Z) %>%
            summarise(
                n_strata = n(),
                min_size = min(N_zk),
                max_size = max(N_zk),
                imbalance = max(N_zk) / min(N_zk),
                min_var = min(var_resid),
                max_var = max(var_resid),
                ratio_var = max(var_resid) / min(var_resid),
                range_var = max(var_resid) - min(var_resid),
                var_mean_resid = var(mean_resid),
                .groups = "drop"
            )
        
    } else {
        strata_stats <- data %>%
            group_by(strata_eval) %>%
            summarise(
                N_k = n(),
                mean_resid = mean(resid),
                var_resid = var(resid),
                sd_resid = sd(resid),
                .groups = "drop"
            ) %>%
            filter(!is.na(var_resid))
        
        if (nrow(strata_stats) == 0) return(NULL)
        
        min_var <- min(strata_stats$var_resid, na.rm = TRUE)
        if (is.na(min_var) || min_var < min_var_threshold) return(NULL)
        
        heterogeneity <- strata_stats %>%
            summarise(
                n_strata = n(),
                min_size = min(N_k),
                max_size = max(N_k),
                imbalance = max(N_k) / min(N_k),
                min_var = min(var_resid),
                max_var = max(var_resid),
                ratio_var = max(var_resid) / min(var_resid),
                range_var = max(var_resid) - min(var_resid),
                var_mean_resid = var(mean_resid)
            )
    }
    
    list(
        strata_stats = strata_stats,
        heterogeneity = heterogeneity,
        stratification = strata_var
    )
}




#' Compare Multiple Candidate Stratification Variables
#'
#' Evaluates and compares multiple candidate stratification variables to help
#' select the best one for stratified model-assisted estimation. Each candidate
#' is assessed using \code{\link{evaluate_stratification}}, and results are
#' summarized in a ranked table.
#'
#' @param data Data frame containing a \code{resid} column (residuals = Yobs - Yhat),
#'   a treatment indicator column named \code{"Z"}, and any columns referenced
#'   by the candidate stratifications.
#' @param candidates Named list of candidate stratification vectors. Each element
#'   should be a vector of the same length as \code{nrow(data)} containing stratum
#'   assignments. Names are used to identify stratifications in the output.
#' @param min_var_threshold Numeric. Minimum within-stratum residual variance
#'   required. Stratifications with any stratum below this threshold are excluded
#'   from results. Default is \code{0.01}.
#'
#' @return A list of performance metrics for each stratification method
#' 
#' @export
compare_stratifications <- function(data, candidates, Z="Z", ...) {
    
    results <- list()
    
    for (cname in names(candidates)) {
        data[[cname]] <- candidates[[cname]]
        
        eval_result <- evaluate_stratification(data, cname, Z = Z, ...)
        
        if (!is.null(eval_result)) {
            results[[cname]] <- eval_result
        }
    }
    
    summary_df <- map_dfr(results, function(res) {
        het <- res$heterogeneity
        if ("Z" %in% names(het)) {
            het <- het %>%
                summarise(across(where(is.numeric), mean))
        }
        
        data.frame(
            stratification = res$stratification,
            n_strata = het$n_strata,
            min_size = het$min_size,
            max_size = het$max_size,
            imbalance_ratio = het$imbalance,
            min_var_resid = het$min_var,
            max_var_resid = het$max_var,
            ratio_var_resid = het$ratio_var,
            var_mean_resid = het$var_mean_resid
        )
    }) %>%
        arrange(desc(ratio_var_resid))
    
    list(
        summary = summary_df,
        details = results
    )
}



#' Evaluate a candidate stratification strategy using pre-coding proxies
#'
#' Since we don't have residuals, we use proxies:
#' - Variance of Yhat within strata (lower = more homogeneous predictions)
#' - Spread of stratum means (higher = better potential for bias correction)
#' - Balance of stratum sizes
#'
#' @param data Data frame with Yhat column
#' @param strata Stratification vector
#' @param Z Character string specifying the name of the treatment indicator  or \code{NULL} to ignore treatment groups. 
#' @return Data frame with proxy quality metrics
#' @export
evaluate_stratification_precoding <- function(data, strata=NULL, Z = "Z") {
    
    if (is.null(Z)){
        data$Z=1
    }
    
    
    # If strata provided as separate vector, add to data
    if (!is.null(strata)) {
        if (length(strata) == nrow(data)) {
            data$strata <- strata
        } else {
            stop("Length of strata vector must match number of rows in data")
        }
    }
    
    # Ensure strata column exists
    if (!"strata" %in% names(data)) {
        stop("No strata column found. Provide strata as argument or include in data.")
    }
    
    
    # Compute statistics by treatment arm and stratum
    strata_stats <- data %>%
        group_by(Z, strata) %>%
        summarise(
            N_zk = n(),
            mean_Yhat = mean(Yhat),
            var_Yhat = var(Yhat),
            .groups = "drop"
        )
    
    # Summarize by treatment arm
    arm_summary <- strata_stats %>% 
        group_by(Z) %>%
        mutate(N_z=sum(N_zk),
               w_k = N_zk / N_z,
               Yhat_bar_z=weighted.mean(mean_Yhat, w_k)) %>%
        summarise(
            N_z = sum(N_zk),
            n_strata = n(),
            
            # Stratum size balance (closer to 1 = more balanced)
            min_size = min(N_zk),
            max_size=max(N_zk),
            size_ratio = max(N_zk) / min(N_zk),
            
            # Spread of stratum means (higher = more separation)
            # This is a proxy for potential BS benefit
            var_means = sum(w_k * (mean_Yhat - Yhat_bar_z)^2),
            range_means = max(mean_Yhat) - min(mean_Yhat),
            
            # Average within-stratum variance of Yhat
            # Lower values suggest more homogeneous strata
            avg_within_var = weighted.mean(var_Yhat,w_k),
            sd_within_var=sqrt(sum(w_k * (var_Yhat - avg_within_var)^2)),
            .groups = "drop"
        )
    
    # Overall summary (single row)
    overall <- tibble(
        n_strata = first(arm_summary$n_strata),
        min_stratum_size = min(arm_summary$min_size),
        max_stratum_size=max(arm_summary$max_size),
        avg_size_ratio = mean(arm_summary$size_ratio),
        
        # Proxy for BS benefit: spread of Yhat means across strata
        # Higher values suggest LLM predictions differ systematically across strata
        # which often correlates with systematic bias in residuals
        mean_spread = mean(arm_summary$var_means),
        
        # Proxy for within-stratum homogeneity
        avg_within_var_Yhat = mean(arm_summary$avg_within_var),
        sd_within_var_Yhat=mean(arm_summary$sd_within_var),
        
        # Overall variance of Yhat (for computing R-squared proxy)
        total_var_Yhat = var(data$Yhat),
    ) 
    
    return(overall)
}


#' Compare candidate stratifications using pre-coding metrics
#'
#' @param data Data frame with Yhat and Z columns
#' @param candidates Named list of stratification vectors
#' @param Z Character string specifying the name of the treatment indicator  or \code{NULL} to ignore treatment groups. 

#' @return Data frame comparing stratifications
#' @export
compare_stratifications_precoding <- function(data, candidates, Z=NULL) {
    
    results <- purrr::map_dfr(names(candidates), function(cname) {
        
        strata <- candidates[[cname]]
        
        metrics <- evaluate_stratification_precoding(data, strata,Z=Z)
        metrics %>%
            mutate(stratification = cname, .before = 1)
    })
    
    if (nrow(results) == 0) {
        warning("No valid stratifications found")
        return(NULL)
    }
    
    # Rank by variance of strata means
    results %>%
        arrange(desc(mean_spread))
}


#' Compute variance decomposition (BS and WS components from Theorem 3.2)
#' 
#' @param dat Data frame with Yobs, Yhat
#' @param Z Treatment indicator (if applicable)
#' @param strata Stratification variable
#' @return Data frame with BS and WS components by treatment arm (if applicable) and overall
#' 
#' @export
compute_variance_decomposition <- function(dat, Z="Z", strata = NULL) {
    
    by_tx = !is.null(Z)
    
    if ( by_tx ) {
        stopifnot( Z %in% colnames(dat) )
        dat$Z = dat[[Z]]
    } else {
        dat$Z = 1
    }
    
    # If strata provided as separate vector, add to data
    if (!is.null(strata)) {
        if (length(strata) == nrow(dat)) {
            dat$strata <- strata
        } else {
            stop("Length of strata vector must match number of rows in dat")
        }
    }
    
    # Ensure strata column exists
    if (!"strata" %in% names(dat)) {
        stop("No strata column found. Provide strata as argument or include in dat.")
    }
    
    
    # Compute residuals if not present
    if (!"resid" %in% names(dat)) {
        dat <- dat %>%
            mutate(resid = Yobs - Yhat)
    }
    
    # -------------------------------------------------------------------------
    # Overall (unstratified) statistics by treatment arm
    # -------------------------------------------------------------------------
    overall_stats <- dat %>%
        group_by(Z) %>%
        summarise(
            N_z = n(),
            mean_resid = mean(resid),
            var_resid_total = var(resid),  # Total residual variance
            .groups = "drop"
        )
    
    # -------------------------------------------------------------------------
    # Within-strata statistics
    # -------------------------------------------------------------------------
    strata_stats <- dat %>%
        group_by(Z, strata) %>%
        summarise(
            N_zk = n(),
            mean_resid_k = mean(resid),
            var_resid_k = var(resid),
            .groups = "drop"
        ) %>%
        group_by(Z) %>%
        mutate(
            N_z = sum(N_zk),
            w_k = N_zk / N_z  # Stratum weight
        ) %>%
        ungroup()
    
    # -------------------------------------------------------------------------
    # Between-Strata Variance (BS)
    # BS = sum_k w_k * (mean_k - grand_mean)^2
    # -------------------------------------------------------------------------
    bs_component <- strata_stats %>%
        left_join(overall_stats %>% select(Z, grand_mean = mean_resid), by = "Z") %>%
        group_by(Z) %>%
        summarise(
            BS = sum(w_k * (mean_resid_k - grand_mean)^2),
            .groups = "drop"
        )
    
    # -------------------------------------------------------------------------
    # Within-Strata Variance (WS)
    # WS = sum_k w_k * var_k (weighted average of stratum variances)
    # -------------------------------------------------------------------------
    ws_component <- strata_stats %>%
        group_by(Z) %>%
        summarise(
            WS = sum(w_k * var_resid_k),
            .groups = "drop"
        )
    
    # -------------------------------------------------------------------------
    # Heterogeneity metrics (for optimal allocation benefit)
    # -------------------------------------------------------------------------
    heterogeneity <- strata_stats %>%
        group_by(Z) %>%
        summarise(
            n_strata = n(),
            min_var = min(var_resid_k),
            max_var = max(var_resid_k),
            var_ratio = max(var_resid_k) / min(var_resid_k),
            # Coefficient of variation of stratum variances
            cv_var = sd(var_resid_k) / mean(var_resid_k),
            .groups = "drop"
        )
    
    # -------------------------------------------------------------------------
    # Combine all components
    # -------------------------------------------------------------------------
    decomposition_by_arm <- overall_stats %>%
        left_join(bs_component, by = "Z") %>%
        left_join(ws_component, by = "Z") %>%
        left_join(heterogeneity, by = "Z") %>%
        mutate(
            # Verify decomposition (should be close to var_resid_total)
            BS_plus_WS = BS + WS,
            
            # Key metrics
            BS_minus_WS = BS - WS,
            pct_BS = BS / var_resid_total * 100,
            pct_WS = WS / var_resid_total * 100,
            
            # Stratification benefit ratio
            BS_over_WS = BS / WS,
            
        )
    
    
    # -------------------------------------------------------------------------
    # Sum over treatment arms
    # -------------------------------------------------------------------------
    
    mean_resid_total = mean(overall_stats$mean_resid)
    S2_e_total <- sum(overall_stats$var_resid_total)
    
    # BS and WS summed over arms
    BS <- sum(bs_component$BS)
    WS <- sum(ws_component$WS)
    
    decomposition_overall <- tibble(
        N = sum(overall_stats$N_z),
        n_strata = first(heterogeneity$n_strata),
        
        # Variance components (summed over arms)
        mean_resid = mean_resid_total,
        var_resid = S2_e_total,
        
        BS = BS,
        WS = WS,
        
        # Verify decomposition
        BS_plus_WS = BS + WS,
        
        # Key metrics
        BS_minus_WS = BS - WS,
        pct_BS = BS / S2_e_total * 100,
        pct_WS = WS / S2_e_total * 100,
        BS_over_WS = BS / WS,
        
        # Heterogeneity (average across arms)
        var_ratio = mean(heterogeneity$var_ratio)
    )
    
    if(by_tx){
        decomposition = list(by_arm=decomposition_by_arm, overall=decomposition_overall)
    }
    else {
        decomposition=decomposition_overall
    }
    
    return(decomposition)
}



#' Create visualization of residual structure by strata
#' 
#' @param dat Data frame with Yobs, Yhat, Z, strata columns
#' @param strata Stratification variable
#' @return ggplot object
#' 
#' @export
plot_strata_residuals <- function(dat, strata) {
    
    dat$strata=strata
    n_strat=n_distinct(dat$strata)
    
    dat <- dat %>%
        mutate(
            resid = Yobs - Yhat,
            strata_label = factor(strata,labels=paste0("Stratum ", 1:n_strat))
        )
    
    if (!is.null(dat$Z)){
        dat = dat %>% mutate( Z_label = factor(Z, levels = c(0, 1), labels = c("Control", "Treatment")))
        
        
        p <- ggplot(dat, aes(x = strata_label, y = resid, fill = Z_label)) +
            geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
            facet_wrap(~Z_label) +
            labs(
                x = "Stratum",
                y = "Residual (Human - LLM)",
                title = "Prediction Residuals by Stratum",
                fill = "Group"
            ) +
            theme_minimal() +
            theme(legend.position = "none")
    }
    
    else if (is.null(dat$Z)){
        
        p <- ggplot(dat, aes(x = strata_label, y = resid)) +
            geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
            labs(
                x = "Stratum",
                y = "Residual (Human - LLM)",
                title = "Prediction Residuals by Stratum") +
            theme_minimal() +
            theme(legend.position = "none")
        
    }
    
    return(p)
}




