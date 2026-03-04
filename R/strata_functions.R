
# Labels for different numbers of groups
get_labels <- function(n_groups, var_name) {
    base_labels <- switch(as.character(n_groups),
                          "2" = c("Low", "High"),
                          "3" = c("Low", "Med", "High"),
                          "4" = c("Low", "Med-Low", "Med-High", "High"),
                          "5" = c("Very Low", "Low", "Med", "High", "Very High"),
                          paste0("Q", 1:n_groups)
    )
    paste(base_labels, var_name)
}


collapse_ordered_adjacent <- function(x, n_group, var_name) {
    x = as.factor(x)
    lv <- levels(x)
    n = length(x)
    
    counts <- table(x)
    grp <- tibble(level = lv, n = as.integer(counts)) %>%
        mutate(level_i = map( lv, ~  .x ) )

    while (nrow(grp) > n_group) {
        i_min <- which.min(grp$n)
        
        if (i_min == 1) {
            j <- 2
        } else if (i_min == nrow(grp)) {
            j <- nrow(grp) - 1
        } else {
            j <- if (grp$n[i_min - 1] <= grp$n[i_min + 1]) i_min - 1 else i_min + 1
        }
        
        i <- min(i_min, j)
        grp$n[i] <- grp$n[i] + grp$n[i + 1]
        grp$level_i[i] <- list( c( grp$level_i[[i]], grp$level_i[[i+1]]) )
        grp = grp[ -(i+1), ]
    }
    
    fmt_label <- function(x) {
        if (length(x) == 1) {
            return(x)
        } else if (length(x) == 2) {
            return( paste0( x[[1]], ",", x[[2]] ) )
        } else {
            return( paste0( x[[1]], "-", x[[length(x)]] ) )
        }
    }
    
    mapping <- tibble(
        orig = unlist(grp$level_i),
        new  = map_chr(grp$level_i, fmt_label) |> rep(times = lengths(grp$level_i))
    )
    mapping$new = paste0( var_name, "=", mapping$new )
    
    out <- mapping$new[match(x, mapping$orig)]
    factor(out, levels = unique(mapping$new))
}


# Create quantile-based strata (returns NULL if breaks aren't unique)
make_quantile_strata <- function(x, n_groups, var_name = "X") {
    
    n <- length(x)
    
    if ( n_distinct( x ) < n_groups ) {
        return(NULL)
    }
    
    # treat x as categorical if factor/character/logical, or if n is NULL/NA
    is_cat <- is.factor(x) || is.character(x) || is.logical(x) || 
        ( n_distinct(x) <= pmin( 30, 5*n_groups ))
    
    if (is_cat) {
        f <- collapse_ordered_adjacent(x, n_groups, var_name)
        return(f)
    }
    
    probs <- seq(0, 1, length.out = n_groups + 1)
    breaks <- quantile(x, probs = probs, na.rm = TRUE)
    if (dplyr::n_distinct(breaks) < n_groups + 1) return(NULL)
    labels <- get_labels(n_groups, var_name)
    cut(x, breaks = breaks, labels = labels, include.lowest = TRUE)
}


# Create quantile-based strata separately by treatment arm
# NOTE: This function is not currently used by any exported functions.
# The pooled stratification approach (make_quantile_strata) is used
# instead, since arm-specific cutpoints can produce inconsistent
# strata labels across treatment groups.
make_quantile_strata_by_arm <- function(x, z_vec, n_groups, var_name = "X") {

    warning( "make_quantile_strata_by_arm is experimental and not currently recommended" )
    
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

#' Create Interaction-Based Stratification (Continuous or Categorical)
#'
#' Constructs an interaction stratification from two variables. Each variable
#' may be continuous (binned into quantile-based strata) or categorical
#' (used directly as factor levels). Optionally, quantile strata can be
#' constructed within treatment arms.
#'
#' If a variable is a factor, character, or logical vector—or if `n` is
#' `NULL` or `NA`—it is treated as categorical. Otherwise, it is treated
#' as continuous and partitioned into `n` quantile strata.
#'
#' @param x1 First stratification variable.
#' @param n1 Number of quantile strata for `x1` (ignored if categorical).
#' @param var_name1 Character string used in labeling strata for `x1`.
#' @param x2 Second stratification variable.
#' @param n2 Number of quantile strata for `x2` (ignored if categorical).
#' @param var_name2 Character string used in labeling strata for `x2`.
#' @param z_vec Optional treatment indicator. If supplied and a variable
#'   is continuous, quantile strata are constructed within treatment arms.
#'
#' @return A factor representing the interaction of the two stratifications,
#'   with empty levels dropped. Returns `NULL` if either stratification fails.
#'
#' @export
make_interaction_strata <- function(x1, n1, var_name1 = NULL, 
                                    x2 = NULL, n2 = NULL, var_name2 = NULL,
                                    z_vec = NULL) {
    
    if ( !is.null(x2) && ( length(x1) != length(x2) ) ) {
        stop("x1 and x2 must be the same length")
    }
    
    make_strata_one <- function(x, n, var_name, z_vec = NULL) {
        # treat as categorical if factor/character/logical, or if n is NULL/NA
        is_cat <- is.factor(x) || is.character(x) || is.logical(x) || is.null(n) || is.na(n) || n_distinct(x) <= n
        
        if (is_cat) {
            f <- as.factor(x)
            levs <- levels(f)
            levels(f) <- paste0(var_name, "=", levs)
            return(f)
        }
        
        if (is.null(z_vec)) {
            make_quantile_strata(x, n, var_name)
        } else {
            make_quantile_strata(x, n, var_name)
            #make_quantile_strata_by_arm(x, z_vec, n, var_name)
        }
    }
    
    if ( is.null( var_name1 ) ) {
        var_name1 <- "Var1"
    }
    
    s1 <- make_strata_one(x1, n1, var_name1, z_vec)
    if ( is.null(s1) ) {
        return(NULL)
    }
    
    if ( !is.null( x2 ) ) {
        if ( is.null( var_name2 ) ) {
            var_name2 <- "Var2"
        }
        s2 <- make_strata_one(x2, n2, var_name2, z_vec)
        
        if ( is.null(s2) ) {
            return(NULL)
        }
        s_int <- interaction(s1, s2, drop = TRUE, sep = " / ")
        stopifnot( length( s_int ) == length(x1) )
        s_int
    } else {
        s1   
    }
}




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
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Helper functions
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    
    
    
    
    # Apply jitter to a variable
    add_jitter <- function(x) x + rnorm(length(x), 0, sd_jitter)
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Validate and extract treatment indicator if provided
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
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
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Prepare variables (with jitter for specified vars)
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    if (is.null(jitter_vars)) jitter_vars <- character(0)
    
    var_data <- lapply(vars, function(v) {
        x <- dat[[v]]
        if (is.null(x)) stop(paste("Variable", v, "not found in dat"))
        if (v %in% jitter_vars) add_jitter(x) else x
    })
    names(var_data) <- vars
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Generate single-variable stratifications (tertiles, quartiles, quintiles)
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    single_strats <- list()
    
    for (v in vars) {
        # Pooled stratification (quantile breaks computed on full sample)
        single_strats[[paste0(v, "_tertiles")]]  <- make_quantile_strata(var_data[[v]], 3, v)
        single_strats[[paste0(v, "_quartiles")]] <- make_quantile_strata(var_data[[v]], 4, v)
        single_strats[[paste0(v, "_quintiles")]] <- make_quantile_strata(var_data[[v]], 5, v)
    }
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Generate interaction stratifications for all pairs of variables
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
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
                if ( length( interaction_strats[[strat_name]] ) != n ) {
                    interaction_strats[[strat_name]] <- NULL
                }
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
    if (nrow(data) == 0) return(NULL)
    
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
            mutate(N_z=sum(N_zk),
                   w_k = N_zk / N_z,
                   mean_resid_z=weighted.mean(mean_resid, w_k)) %>%
            summarise(
                N_z = sum(N_zk),
                n_strata = n(),
                min_size = min(N_zk),
                max_size = max(N_zk),
                imbalance = max(N_zk) / min(N_zk),
                min_var = min(var_resid),
                max_var = max(var_resid),
                ratio_var = max(var_resid) / min(var_resid),
                range_var = max(var_resid) - min(var_resid),
                avg_within_var = weighted.mean(var_resid, w_k),
                sd_within_var = sqrt(sum(w_k * (var_resid - avg_within_var)^2)),
                var_mean_resid = sum(w_k * (mean_resid - mean_resid_z)^2),
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
            mutate(N=sum(N_k),
                   w_k = N_k / N,
                   avg_mean_resid=weighted.mean(mean_resid, w_k)) %>%
            summarise(
                n_strata = n(),
                min_size = min(N_k),
                max_size = max(N_k),
                imbalance = max(N_k) / min(N_k),
                min_var = min(var_resid),
                max_var = max(var_resid),
                ratio_var = max(var_resid) / min(var_resid),
                range_var = max(var_resid) - min(var_resid),
                avg_within_var = weighted.mean(var_resid, w_k),
                sd_within_var = sqrt(sum(w_k * (var_resid - avg_within_var)^2)),
                var_mean_resid = sum(w_k * (mean_resid - avg_mean_resid)^2),
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
#'   a treatment indicator column, and any columns referenced by the candidate
#'   stratifications.
#' @param candidates Named list of candidate stratification vectors. Each element
#'   should be a vector of the same length as \code{nrow(data)} containing stratum
#'   assignments. Names are used to identify stratifications in the output.
#' @param Z Character string. Name of the treatment indicator column in
#'   \code{data}. Default is \code{"Z"}.
#' @param ... Additional arguments passed to \code{\link{evaluate_stratification}}.
#'
#' @return A data frame of performance metrics for each stratification method,
#'   sorted by promise.
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
            het <- het %>% summarise(n_strata=mean(n_strata),
                                     min_size = min(min_size),
                                     max_size = max(max_size),
                                     imbalance=mean(imbalance),
                                     min_var=mean(min_var),
                                     max_var=mean(max_var),
                                     ratio_var=mean(ratio_var),
                                     range_var=mean(range_var),
                                     avg_within_var=mean(avg_within_var),
                                     sd_within_var=mean(sd_within_var),
                                     var_mean_resid=mean(var_mean_resid))
        }
        
        data.frame(
            stratification = res$stratification,
            n_strata = het$n_strata,
            min_stratum_size = het$min_size,
            max_stratum_size = het$max_size,
            avg_size_ratio = het$imbalance,
            mean_spread = het$var_mean_resid,
            avg_within_var_resid = het$avg_within_var,
            sd_within_var_resid = het$sd_within_var,
            min_var_resid = het$min_var,
            max_var_resid = het$max_var,
            ratio_var_resid = het$ratio_var
        )
    }) %>%
        arrange(desc(mean_spread))
    
    list(
        summary = tibble(summary_df),
        details = results
    )
}


#' Evaluate a candidate stratification strategy using pre-coding
#' proxies
#'
#' Since we don't have residuals, we use proxies:
#' - Variance of Yhat within strata (lower = more homogeneous predictions)
#' - Spread of stratum means (higher = better potential for bias correction)
#' - Balance of stratum sizes
#'
#' @param data Data frame with Yhat column
#' @param strata Stratification vector
#' @param Z Character string specifying the name of the treatment
#'   indicator  or \code{NULL} to ignore treatment groups.
#' @param Yhat Character string specifying the name of the column
#'   containing predicted values (Yhat). Default is \code{"Yhat"}.
#' @return Data frame with proxy quality metrics
#' @export
evaluate_stratification_precoding <- function( data, strata=NULL, Z = "Z",
                                               Yhat = "Yhat" ) {
    
    if (is.null(Z)){
        data$Z=1
    } else {
        data$Z = data[[Z]]
    }
    
    data$Yhat = data[[Yhat]]
    
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
    
    if (nrow(strata_stats) == 0) return(NULL)
    
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
            sd_within_var=sqrt(sum(w_k * (var_Yhat - avg_within_var)^2))
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
#' @param Z Character string specifying the name of the treatment
#'   indicator  or \code{NULL} to ignore treatment groups.
#' @param Yhat Character string specifying the name of the column
#'   containing predicted values (Yhat). Default is \code{"Yhat"}.
#' @return Data frame comparing stratifications
#' @export
compare_stratifications_precoding <- function(data, candidates, Z=NULL,
                                              Yhat = "Yhat" ) {
    
    stopifnot( is.data.frame(data), 
               Yhat %in% colnames(data),
               all( sapply(candidates, length) == nrow(data)),
               (is.null(Z) || Z %in% names(data)) )
    
    
    results <- purrr::map_dfr(names(candidates), function(cname) {
        
        strata <- candidates[[cname]]
        
        metrics <- evaluate_stratification_precoding(data, strata, Z=Z, Yhat = Yhat)
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


#' Compute variance decomposition (BS and WS components)
#'
#' Computes the between-strata (BS) and within-strata (WS) variance
#' components of the residuals. These components are related to the
#' variance reduction from stratified sampling as described in Theorem
#' 3.2 of the manuscript.
#'
#' @param dat Data frame with Yobs, Yhat columns
#' @param Z Character string specifying the treatment indicator column
#'   name, or \code{NULL} to ignore treatment groups. Default is
#'   \code{"Z"}.
#' @param strata Stratification variable. Can be a vector of stratum
#'   assignments (same length as \code{nrow(dat)}) or \code{NULL} if
#'   \code{dat} already contains a column named \code{"strata"}.
#' @param n_code Optional. Total number of units to be sampled for
#'   coding. If provided, the function computes the exact BS and WS
#'   terms from Theorem 3.2 which depend on the sampling fraction. If
#'   \code{NULL} (default), computes simplified population-level
#'   variance decomposition components that serve as proxies.
#'
#'
#' @return If \code{Z} is specified (treatment/control), returns a
#'   list with:
#'   \describe{
#'     \item{\code{by_arm}}{Data frame with decomposition for each treatment arm.}
#'     \item{\code{overall}}{Data frame with overall decomposition (summed over arms).}
#'   }
#'   If \code{Z} is \code{NULL}, returns a single data frame with the
#'   decomposition.  See details for the specific components included in each case.
#'
#'

#'
#'
#' @details 
#' The variance reduction from stratified sampling (Theorem
#' 3.2) is: 
#' 
#' \deqn{Var(\hat{\tau}_{SRS}) - Var(\hat{\tau}_{strat}) = BS - WS}
#'
#' \strong{Simplified (population-level) decomposition} (when
#' \code{n_code = NULL}):
#'
#' This computes the standard ANOVA-style decomposition of residual
#' variance:
#' \itemize{
#'   \item \code{BS = sum_k w_k * (mean_resid_k - grand_mean)^2}
#'   \item \code{WS = sum_k w_k * var_resid_k}
#' }
#' where \code{w_k = N_zk / N_z} is the stratum weight. This satisfies:
#' \deqn{Var(resid) = BS + WS}
#'
#'
#' These simplified components serve as useful diagnostics:
#' \itemize{
#'   \item Large BS suggests systematic differences in prediction bias across
#'     strata, which stratification can exploit.
#'   \item Heterogeneous within-stratum variances (captured by \code{var_ratio})
#'     suggest potential benefit from optimal allocation.
#' }
#'
#'
#' \strong{Exact Theorem 3.2 decomposition} (when \code{n_code} is
#' provided):
#'
#' The exact formulas incorporate finite population corrections and
#' sampling fractions. The BS term from Theorem 3.2 is: 
#' \deqn{BS = sum_k [ (N_z - n_z)/N_z * N_zk / (n_z(N_z-1)) * (mean_resid_k - grand_mean)^2 ]}
#'
#' The WS term is more complex and depends on the sampling design.
#' Under proportional allocation, WS is approximately zero. Under
#' optimal allocation (oversampling high-variance strata), WS can
#' become negative, providing additional variance reduction.
#'
#' The returned decomposition includes:
#' \describe{
#'     \item{\code{BS}}{Between-strata variance component. Measures how much
#'       stratum mean residuals differ from the grand mean. Larger values indicate
#'       greater potential benefit from stratification.}
#'     \item{\code{WS}}{Within-strata variance component. Weighted average of
#'       within-stratum residual variances. Under proportional allocation, this
#'       term is typically small and positive.  NOTE: The exact WS term from Theorem 3.2 can be negative under optimal
#'       allocation (oversampling high-variance strata), but the simplified
#'       population version computed here is always positive.}
#'     \item{\code{BS_minus_WS}}{The difference BS - WS, which approximates
#'       the variance reduction from stratified sampling vs. SRS. Positive values
#'       indicate stratification reduces variance.}
#' }
#'
#'
#' @seealso \code{\link{evaluate_stratification}} for computing
#'   stratification quality metrics,
#'   \code{\link{plot_strata_residuals}} for visualization.
#'
#' @export
compute_variance_decomposition <- function(dat, Z = "Z", strata = NULL, n_code = NULL) {
    
    by_tx <- !is.null(Z)
    
    if (by_tx) {
        stopifnot(Z %in% colnames(dat))
        dat$Z <- dat[[Z]]
    } else {
        dat$Z <- 1
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
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Overall (unstratified) statistics by treatment arm
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    overall_stats <- dat %>%
        group_by(Z) %>%
        summarise(
            N_z = n(),
            mean_resid = mean(resid),
            var_resid_total = var(resid),  # Total residual variance (sample variance)
            .groups = "drop"
        )
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Within-strata statistics
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    strata_stats <- dat %>%
        group_by(Z, strata) %>%
        summarise(
            N_zk = n(),
            mean_resid_k = mean(resid),
            var_resid_k = var(resid),  # Sample variance within stratum
            .groups = "drop"
        ) %>%
        group_by(Z) %>%
        mutate(
            N_z = sum(N_zk),
            w_k = N_zk / N_z  # Stratum weight (proportion of arm in stratum)
        ) %>%
        ungroup()
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Compute BS and WS components ----
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    if (is.null(n_code)) {
        # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # Simplified (population-level) decomposition
        # Standard ANOVA-style variance decomposition
        # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        
        # Between-Strata Variance (BS)
        # BS = sum_k w_k * (mean_k - grand_mean)^2
        # This is the variance of stratum means around the grand mean
        bs_component <- strata_stats %>%
            left_join(overall_stats %>% select(Z, grand_mean = mean_resid), by = "Z") %>%
            group_by(Z) %>%
            summarise(
                BS = sum(w_k * (mean_resid_k - grand_mean)^2),
                .groups = "drop"
            )
        
        # Within-Strata Variance (WS)
        # WS = sum_k w_k * var_k (weighted average of stratum variances)
        ws_component <- strata_stats %>%
            group_by(Z) %>%
            summarise(
                WS = sum(w_k * var_resid_k),
                .groups = "drop"
            )
        
    } else {
        # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # Exact Theorem 3.2 decomposition (with sampling fractions)
        # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        
        # Determine n_z (samples per arm) assuming proportional split
        n_arms <- length(unique(dat$Z))
        n_per_arm <- n_code / n_arms
        
        # Between-Strata Variance (BS) - Theorem 3.2 formula
        # BS = sum_k [ (N_z - n_z)/N_z * N_zk / (n_z(N_z-1)) * (mean_k - grand_mean)^2 ]
        bs_component <- strata_stats %>%
            left_join(overall_stats %>% select(Z, N_z_total = N_z, grand_mean = mean_resid), by = "Z") %>%
            mutate(
                n_z = n_per_arm,
                # Finite population correction factor
                fpc = (N_z_total - n_z) / N_z_total,
                # Weight for BS term
                bs_weight = fpc * N_zk / (n_z * (N_z_total - 1))
            ) %>%
            group_by(Z) %>%
            summarise(
                BS = sum(bs_weight * (mean_resid_k - grand_mean)^2),
                .groups = "drop"
            )
        
        # Within-Strata Variance (WS) - Theorem 3.2 formula
        # Under proportional allocation: n_zk / N_zk approx = n_z / N_z for all k
        # WS = sum_k [ (N_zk/N_z^2 * (N_zk - n_zk)/n_zk - (N_zk-1)/(N_z(N_z-1)) * (N_z - n_z)/n_z) * var_k ]
        # Under proportional allocation, this simplifies to approximately 0
        ws_component <- strata_stats %>%
            left_join(overall_stats %>% select(Z, N_z_total = N_z), by = "Z") %>%
            mutate(
                n_z = n_per_arm,
                # Proportional allocation: n_zk = n_z * (N_zk / N_z)
                n_zk = pmax(1, round(n_z * w_k)),  # At least 1 per stratum
                # First term: (N_zk/N_z^2) * (N_zk - n_zk) / n_zk
                term1 = (N_zk / N_z_total^2) * (N_zk - n_zk) / n_zk,
                # Second term: (N_zk - 1) / (N_z * (N_z - 1)) * (N_z - n_z) / n_z
                term2 = (N_zk - 1) / (N_z_total * (N_z_total - 1)) * (N_z_total - n_z) / n_z,
                # WS weight (can be positive or negative)
                ws_weight = term1 - term2
            ) %>%
            group_by(Z) %>%
            summarise(
                WS = sum(ws_weight * var_resid_k),
                .groups = "drop"
            )
    }
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Heterogeneity metrics (for optimal allocation benefit)
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
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
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Combine all components by treatment arm
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    decomposition_by_arm <- overall_stats %>%
        left_join(bs_component, by = "Z") %>%
        left_join(ws_component, by = "Z") %>%
        left_join(heterogeneity, by = "Z") %>%
        mutate(
            # Verify decomposition (for simplified version, BS + WS should equal var_resid_total)
            BS_plus_WS = BS + WS,
            
            # Key metric: variance reduction from stratification
            # Positive values indicate stratification reduces variance
            BS_minus_WS = BS - WS,
            
            # Percentage of total variance explained by between-strata differences
            pct_BS = BS / var_resid_total * 100,
            pct_WS = WS / var_resid_total * 100,
            
            # Stratification benefit ratio (BS/WS > 1 suggests benefit)
            BS_over_WS = BS / WS
        )
    
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Sum over treatment arms for overall decomposition
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    mean_resid_total <- mean(overall_stats$mean_resid)
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
        
        # Key metric: variance reduction from stratification
        BS_minus_WS = BS - WS,
        pct_BS = BS / S2_e_total * 100,
        pct_WS = WS / S2_e_total * 100,
        BS_over_WS = BS / WS,
        
        # Heterogeneity (average across arms)
        var_ratio = mean(heterogeneity$var_ratio)
    )
    
    if (by_tx) {
        decomposition <- list(by_arm = decomposition_by_arm, overall = decomposition_overall)
    } else {
        decomposition <- decomposition_overall
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


#' Compare Pre-coding Proxy Metrics with Oracle Metrics
#'
#' Creates a combined table and visualization comparing the pre-coding proxy
#' metrics (based only on Yhat) with oracle metrics (based on actual residuals)
#' for multiple candidate stratifications. This validates whether the pre-coding
#' metrics successfully identify stratifications that perform well in practice.
#'
#' @param oracle_comparison Output from \code{\link{compare_stratifications}},
#'   which uses actual residuals (Yobs - Yhat).
#' @param precoding_comparison Output from \code{\link{compare_stratifications_precoding}},
#'   which uses only Yhat as a proxy.
#' @param top_n Integer specifying how many top stratifications to include
#'   in the output. Default is 10.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{combined_table}}{Data frame merging oracle and pre-coding metrics
#'       for each stratification, with columns for both metric types and rank
#'       comparisons.}
#'     \item{\code{correlation}}{Spearman correlation between oracle and pre-coding

#'       rankings, measuring how well pre-coding metrics predict oracle performance.}
#'     \item{\code{plot}}{A ggplot object showing the relationship between
#'       pre-coding and oracle metrics.}
#'   }
#'
#' @details
#' This function is useful for validating that the pre-coding proxy metrics
#' (which can be computed before any human coding is done) successfully identify
#' stratifications that would perform well according to the oracle metrics
#' (which require knowing the true residuals).
#'
#' The key comparison is between:
#' \itemize{
#'   \item \code{mean_spread} (oracle): Variance of mean residuals across strata (BS component)
#'   \item \code{mean_spread} (pre-coding): Variance of mean Yhat across strata (proxy for BS)
#' }
#'
#' High correlation between rankings suggests the pre-coding metrics are reliable
#' for selecting stratifications in practice.
#'
#' @seealso \code{\link{compare_stratifications}}, \code{\link{compare_stratifications_precoding}}
#'
#' @export
compare_oracle_vs_precoding <- function(oracle_comparison, precoding_comparison, top_n = 10) {
    
    # Extract summary tables
    oracle_df <- oracle_comparison$summary
    precoding_df <- precoding_comparison
    
    # Add rank columns
    oracle_df <- oracle_df %>%
        dplyr::mutate(oracle_rank = dplyr::row_number())
    
    precoding_df <- precoding_df %>%
        dplyr::mutate(precoding_rank = dplyr::row_number())
    
    # Merge on stratification name
    combined <- dplyr::inner_join(
        oracle_df %>%
            dplyr::select(stratification, n_strata,
                          oracle_mean_spread = mean_spread,
                          oracle_avg_within_var = avg_within_var_resid,
                          oracle_ratio_var = ratio_var_resid,
                          oracle_rank),
        precoding_df %>%
            dplyr::select(stratification,
                          precoding_mean_spread = mean_spread,
                          precoding_avg_within_var = avg_within_var_Yhat,
                          precoding_rank),
        by = "stratification"
    )
    
    # Calculate rank difference
    combined <- combined %>%
        dplyr::mutate(
            rank_diff = precoding_rank - oracle_rank,
            abs_rank_diff = abs(rank_diff)
        ) %>%
        dplyr::arrange(oracle_rank)
    
    # Calculate Spearman correlation between rankings
    spearman_cor <- cor(combined$oracle_rank, combined$precoding_rank, method = "spearman")
    
    # Calculate correlation between the actual metric values
    metric_cor <- cor(combined$oracle_mean_spread, combined$precoding_mean_spread,
                      method = "spearman", use = "complete.obs")
    
    # Create comparison plot
    plot_data <- combined %>%
        dplyr::slice_head(n = top_n)
    
    # Scatter plot of oracle vs precoding mean_spread
    p1 <- ggplot2::ggplot(combined,
                          ggplot2::aes(x = precoding_mean_spread, y = oracle_mean_spread)) +
        ggplot2::geom_point(ggplot2::aes(color = n_strata), size = 3, alpha = 0.7) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "gray40", linetype = "dashed") +
        #ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "red") +
        ggplot2::labs(
            x = expression("Pre-coding Proxy: " * Var(bar(hat(Y))[k])),
            y = expression("Oracle: " * Var(bar(e)[k])),
            title = "Validation of Pre-coding Proxy vs. Oracle Stratification Metrics",
            subtitle = sprintf("Spearman correlation = %.3f", metric_cor),
            color = "# Strata"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "right") +
        ggplot2::scale_color_viridis_c()
    
    
    list(
        combined_table = combined,
        rank_correlation = spearman_cor,
        metric_correlation = metric_cor,
        scatter_plot = p1    )
}


#' Compare Candidate Stratifications Using Oracle Variance Decomposition
#'
#' For each candidate stratification, computes the BS (between-strata) and WS
#' (within-strata) variance components from Theorem 3.2, which together determine
#' the variance reduction from stratified sampling. Returns a summary table and
#' visualization comparing candidates based on oracle information (true residuals).
#'
#' @param data Data frame containing \code{Yobs}, \code{Yhat}, and optionally
#'   a treatment indicator column.
#' @param candidates Named list of candidate stratification vectors. Each element
#'   should be a factor or character vector of the same length as \code{nrow(data)}
#'   containing stratum assignments.
#' @param n_code Total number of units to sample for human coding. Required for
#'   computing the exact Theorem 3.2 formulas.
#' @param Z Character string specifying the name of the treatment indicator column,
#'   or \code{NULL} to ignore treatment groups. Default is \code{"Z"}.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{summary}}{Data frame with one row per stratification containing:
#'       \itemize{
#'         \item \code{stratification}: Name of the stratification
#'         \item \code{n_strata}: Number of strata
#'         \item \code{BS}: Between-strata variance component
#'         \item \code{WS}: Within-strata variance component
#'         \item \code{BS_minus_WS}: Variance reduction (BS - WS)
#'         \item \code{var_ratio}: Ratio of max to min within-stratum variance
#'         \item \code{rank}: Rank by \code{BS_minus_WS} (1 = best)
#'       }}
#'     \item{\code{decomposition_details}}{List of full decomposition results
#'       for each stratification (output from \code{compute_variance_decomposition}).}
#'   }
#'
#' @details
#' This function uses oracle information (the true residuals Yobs - Yhat) to
#' compute the exact variance decomposition for each candidate stratification.
#' The key metric is BS - WS, which represents the variance reduction achieved
#' by stratified sampling compared to simple random sampling (Theorem 3.2):
#'
#' \deqn{Var(\hat{\tau}_{SRS}) - Var(\hat{\tau}_{strat}) = BS - WS}
#'
#' Stratifications with larger BS - WS values provide greater variance reduction.
#' The BS component measures systematic differences in prediction bias across
#' strata, while the WS component (typically small under proportional allocation)
#' reflects within-stratum residual heterogeneity.

#' @seealso \code{\link{compute_variance_decomposition}} for the underlying
#'   decomposition.
#'
#' @examples
#' \dontrun{
#' # Create candidate stratifications
#' candidates <- create_stratification_candidates(data, vars = c("Yhat", "word_count"))
#'
#' # Compare using oracle decomposition
#' results <- compare_decomposition_oracle(data, candidates, n_code = 200)
#'
#' # View summary table
#' results$summary
#' }
#' @export
compare_decomposition_oracle <- function(data, candidates, n_code, Z = "Z") {
    
    # Validate inputs
    stopifnot("Yobs" %in% names(data))
    stopifnot("Yhat" %in% names(data))
    stopifnot(length(candidates) > 0)
    
    # Compute residuals
    data$resid <- data$Yobs - data$Yhat
    
    # Storage for results
    summary_list <- list()
    decomposition_list <- list()
    
    for (cname in names(candidates)) {
        
        # Compute decomposition for this stratification
        decomp <- compute_variance_decomposition(
            dat = data,
            Z = Z,
            strata = candidates[[cname]],
            n_code = n_code
        )
        
        # Extract overall summary (or single summary if no treatment)
        if (is.list(decomp) && "overall" %in% names(decomp)) {
            overall <- decomp$overall
        } else {
            overall <- decomp
        }
        
        # Store summary row
        summary_list[[cname]] <- tibble::tibble(
            stratification = cname,
            n_strata = overall$n_strata,
            BS = overall$BS,
            WS = overall$WS,
            BS_minus_WS = overall$BS_minus_WS,
            BS_plus_WS = overall$BS_plus_WS,
            pct_BS = overall$pct_BS,
            var_ratio = overall$var_ratio,
            var_resid = overall$var_resid
        )
        
        # Store full decomposition
        decomposition_list[[cname]] <- decomp
    }
    
    # Combine into summary table
    summary_df <- dplyr::bind_rows(summary_list) %>%
        dplyr::arrange(dplyr::desc(BS_minus_WS)) %>%
        dplyr::mutate(rank = dplyr::row_number())
    
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Return results
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    list(
        summary = summary_df,
        decomposition_details = decomposition_list
    )
}


#' Plot Oracle Decomposition Comparison
#'
#' Creates a horizontal bar chart comparing the variance reduction (BS - WS)
#' across candidate stratifications, with bars colored by the relative
#' contribution of BS vs WS.
#'
#' @param decomposition_result Output from \code{\link{compare_decomposition_oracle}}.
#' @param top_n Integer specifying how many top stratifications to display.
#'   Default is \code{NULL} (show all).
#' @param show_components Logical. If \code{TRUE}, shows separate bars for BS
#'   and WS. If \code{FALSE} (default), shows only BS - WS as a single bar.
#'
#' @return A ggplot object.
#'
#' @export
plot_decomposition_comparison <- function(decomposition_result, top_n = NULL,
                                          show_components = FALSE) {
    
    summary_df <- decomposition_result$summary
    
    # Filter to top_n if specified
    if (!is.null(top_n) && top_n < nrow(summary_df)) {
        summary_df <- summary_df %>%
            dplyr::slice_head(n = top_n)
    }
    
    if (show_components) {
        # Grouped bar chart with BS and WS separate
        plot_data <- summary_df %>%
            tidyr::pivot_longer(
                cols = c(BS, WS),
                names_to = "component",
                values_to = "value"
            ) %>%
            dplyr::mutate(
                component = factor(component, levels = c("BS", "WS")),
                stratification = factor(stratification,
                                        levels = rev(summary_df$stratification))
            )
        
        p <- ggplot2::ggplot(plot_data,
                             ggplot2::aes(x = stratification, y = sqrt(value), fill = component)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.7) +
            ggplot2::coord_flip() +
            ggplot2::labs(
                title = "Variance Decomposition Components",
                x = NULL,
                y = "Variance Component Value",
                fill = "Component"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "bottom") +
            ggplot2::scale_fill_manual(values = c("BS" = "#2166AC", "WS" = "#B2182B"))
        
    } else {
        # Single bar chart showing BS - WS
        plot_data <- summary_df %>%
            dplyr::mutate(
                stratification = factor(stratification,
                                        levels = rev(stratification)),
                fill_color = ifelse(BS_minus_WS > 0, "positive", "negative")
            )
        
        p <- ggplot2::ggplot(plot_data,
                             ggplot2::aes(x = stratification, y = sqrt(BS_minus_WS),
                                          fill = fill_color)) +
            ggplot2::geom_bar(stat = "identity", width = 0.7) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
            ggplot2::coord_flip() +
            ggplot2::labs(
                title = "Variance Reduction from Stratification (BS - WS)",
                subtitle = "Positive values indicate stratification reduces variance vs. SRS",
                x = NULL,
                y = "BS - WS (Variance Reduction)"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none") +
            ggplot2::scale_fill_manual(values = c("positive" = "#2166AC", "negative" = "#B2182B"))
        
        # Add value labels
        #p <- p +
        #  ggplot2::geom_text(
        #    ggplot2::aes(label = sprintf("%.4f", BS_minus_WS),
        #                   hjust = ifelse(BS_minus_WS > 0, -0.1, 1.1)),
        #     size = 3
        # )
    }
    
    return(p)
}


