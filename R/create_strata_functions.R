
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
        new  = map_chr(grp$level_i, fmt_label) %>% rep(times = lengths(grp$level_i))
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


# Quantile strata with arm-specific breaks and common rank labels.
# Continuous x: breaks are the within-arm quantiles, so every arm is
# split into n_groups equally sized strata; the labels ("Low X", ...)
# denote within-arm quantile ranks. Categorical x (same rule as
# make_quantile_strata) is collapsed on the pooled sample so its labels
# stay comparable across arms. Returns NULL if the breaks in any arm are
# not unique. Falls back to pooled cuts when z_vec is NULL.
make_quantile_strata_by_arm <- function(x, z_vec, n_groups, var_name = "X") {
    
    if ( is.null(z_vec) ) {
        return( make_quantile_strata(x, n_groups, var_name) )
    }
    stopifnot( length(x) == length(z_vec) )
    
    if ( n_distinct(x) < n_groups ) {
        return(NULL)
    }
    
    is_cat <- is.factor(x) || is.character(x) || is.logical(x) ||
        ( n_distinct(x) <= pmin( 30, 5 * n_groups ) )
    if ( is_cat ) {
        return( collapse_ordered_adjacent(x, n_groups, var_name) )
    }
    
    probs  <- seq(0, 1, length.out = n_groups + 1)
    labels <- get_labels(n_groups, var_name)
    result <- factor(rep(NA_character_, length(x)), levels = labels)
    
    for ( z_val in unique(z_vec) ) {
        idx    <- which(z_vec == z_val)
        breaks <- quantile(x[idx], probs = probs, na.rm = TRUE)
        if ( dplyr::n_distinct(breaks) < n_groups + 1 ) return(NULL)
        result[idx] <- cut(x[idx], breaks = breaks, labels = labels,
                           include.lowest = TRUE)
    }
    
    result
}



#' Create Interaction-Based Stratification (Continuous or Categorical)
#'
#' Constructs an interaction stratification from two variables. Each variable
#' may be continuous (binned into quantile-based strata) or categorical
#' (used directly as factor levels). When \code{z_vec} is supplied and
#' \code{cut_by_arm = TRUE}, quantile strata for continuous variables are
#' constructed within treatment arms.
#'
#' If a variable is a factor, character, or logical vector, or if `n` is
#' `NULL` or `NA`, it is treated as categorical. Otherwise, it is treated
#' as continuous and partitioned into `n` quantile strata.
#'
#' @param x1 First stratification variable.
#' @param n1 Number of quantile strata for `x1` (ignored if categorical).
#' @param var_name1 Character string used in labeling strata for `x1`.
#' @param x2 Second stratification variable.
#' @param n2 Number of quantile strata for `x2` (ignored if categorical).
#' @param var_name2 Character string used in labeling strata for `x2`.
#' @param z_vec Optional treatment indicator. If supplied, a variable is
#'   continuous, and \code{cut_by_arm = TRUE}, its quantile breaks are
#'   computed within each treatment arm and the labels denote within-arm
#'   quantile ranks.
#' @param cut_by_arm Logical; cut continuous variables at within-arm
#'   quantiles when \code{z_vec} is supplied. Default \code{TRUE}.
#'
#' @return A factor representing the interaction of the two stratifications,
#'   with empty levels dropped. Returns `NULL` if either stratification fails.
#'
#' @export
make_interaction_strata <- function(x1, n1, var_name1 = NULL,
                                    x2 = NULL, n2 = NULL, var_name2 = NULL,
                                    z_vec = NULL, cut_by_arm = TRUE) {        # CHANGED
    
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
        
        if ( is.null(z_vec) || !cut_by_arm ) {                                # CHANGED
            make_quantile_strata(x, n, var_name)
        } else {
            make_quantile_strata_by_arm(x, z_vec, n, var_name)                # CHANGED
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



# Cut a numeric variable at user-supplied cutpoints. These are fixed
# values, so the cut is the same in every arm. NULL if fewer than two
# non-empty groups result.
make_break_strata <- function(x, breaks, var_name = "X") {
    if ( !is.numeric(x) ) {
        return(NULL)
    }
    breaks <- sort( unique( breaks ) )
    if ( length(breaks) < 3 ) {
        return(NULL)
    }
    f <- droplevels( cut( x, breaks = breaks, include.lowest = TRUE, dig.lab = 4 ) )
    if ( nlevels(f) < 2 ) {
        return(NULL)
    }
    levels(f) <- paste0( var_name, "=", levels(f) )
    f
}


# Use a variable (or an externally built stratification) at its own
# levels, with "var=" prefixed labels. NULL if fewer than two levels.
make_level_strata <- function(x, var_name = "X") {
    f <- droplevels( as.factor( x ) )
    if ( nlevels(f) < 2 ) {
        return(NULL)
    }
    levels(f) <- paste0( var_name, "=", levels(f) )
    f
}


# Cross two stratifications (same construction make_interaction_strata
# uses); NULL if either is NULL.
cross_strata <- function(s1, s2) {
    if ( is.null(s1) || is.null(s2) ) {
        return(NULL)
    }
    stopifnot( length(s1) == length(s2) )
    interaction( s1, s2, drop = TRUE, sep = " / " )
}


# Canonical key of a partition: two stratifications get the same key
# exactly when they group the units identically, whatever the labels or
# level order.
partition_key <- function(s) {
    s <- as.character(s)
    paste( match( s, unique(s) ), collapse = "," )
}


# Screen a named list of candidate stratifications, in order, for
# (i) failed construction (NULL), (ii) wrong length, (iii) missing
# values, (iv) fewer than two strata, (v) a cell smaller than
# min_cell_size (cell = stratum-by-arm when z_vec is supplied, stratum
# otherwise; structurally empty cells count as 0), and (vi) duplicate
# partitions (the first candidate with a given partition keeps its
# name). Returns list(candidates = <kept, as factors>, screening = <tibble>).
screen_candidates <- function(candidates, n, z_vec = NULL,
                              min_cell_size = 2, drop_duplicates = TRUE) {
    
    empty_log <- tibble::tibble( candidate = character(0), n_strata = integer(0),
                                 min_cell = integer(0), kept = logical(0),
                                 reason = character(0) )
    if ( length(candidates) == 0 ) {
        return( list( candidates = list(), screening = empty_log ) )
    }
    
    if ( is.null(z_vec) ) {
        z_vec <- rep( 1L, n )
    }
    
    seen <- character(0)   # partition key -> name of the first candidate with it
    rows <- vector( "list", length(candidates) )
    
    for ( i in seq_along(candidates) ) {
        nm <- names(candidates)[i]
        s  <- candidates[[i]]
        n_strata <- NA_integer_
        min_cell <- NA_integer_
        kept <- FALSE
        reason <- ""
        
        if ( is.null(s) ) {
            reason <- "construction failed (non-unique breaks or too few distinct values)"
        } else if ( length(s) != n ) {
            reason <- "length does not match nrow(dat)"
        } else if ( anyNA(s) ) {
            reason <- "contains missing values"
        } else {
            s <- droplevels( as.factor(s) )
            n_strata <- nlevels(s)
            min_cell <- as.integer( min( table( s, z_vec ) ) )
            if ( n_strata < 2 ) {
                reason <- "fewer than 2 strata"
            } else if ( min_cell < min_cell_size ) {
                reason <- sprintf( "smallest cell has %d unit(s) (min_cell_size = %d)",
                                   min_cell, min_cell_size )
            } else {
                key <- partition_key(s)
                if ( drop_duplicates && key %in% names(seen) ) {
                    reason <- sprintf( "duplicate of '%s'", seen[[key]] )
                } else {
                    seen[[key]] <- nm
                    kept <- TRUE
                    candidates[[i]] <- s
                }
            }
        }
        
        rows[[i]] <- tibble::tibble( candidate = nm, n_strata = n_strata,
                                     min_cell = min_cell, kept = kept,
                                     reason = reason )
    }
    
    screening <- dplyr::bind_rows( rows )
    list( candidates = candidates[ screening$kept ], screening = screening )
}



#' Create Candidate Stratification Variables
#'
#' Generates a named list of candidate stratifications for use with
#' \code{\link{compare_stratifications}},
#' \code{\link{compare_stratifications_precoding}}, and the sampling
#' functions. By default this reproduces the earlier library (tertiles,
#' quartiles, and quintiles of each variable plus 2x2, 2x3, 3x2, and 3x3
#' crossings of each pair of variables) with three changes: quantile cuts
#' are made within each treatment arm when \code{Z} is supplied,
#' categorical variables are also offered at their native levels, and
#' every candidate is screened so that each stratum-by-arm cell holds at
#' least \code{min_cell_size} units and no two candidates induce the same
#' partition of the units. The library can be extended with other group
#' counts and crossings, user-supplied cutpoints, and externally built
#' stratifications.
#'
#' @param dat Data frame containing the variables to stratify on.
#' @param vars Character vector of column names in \code{dat} to build
#'   stratifications from.
#' @param Z Optional name of the treatment indicator column in \code{dat}.
#'   When supplied, quantile cuts are made within each arm if
#'   \code{cut_by_arm = TRUE}, and the \code{min_cell_size} screen is
#'   applied to stratum-by-arm cells rather than to strata. Default
#'   \code{NULL}.
#' @param cut_by_arm Logical. If \code{TRUE} (default) and \code{Z} is
#'   supplied, quantile breaks for continuous variables are computed
#'   separately within each treatment arm, so every arm is divided into
#'   equally sized strata and the labels denote within-arm quantile ranks
#'   ("Low Yhat" is the bottom tertile of predicted scores within the
#'   arm). If \code{FALSE}, breaks are computed on the pooled sample and
#'   define the same subgroups in both arms (needed, for example, for
#'   stratum-specific treatment effects). Categorical variables and
#'   \code{breaks} cuts are never cut by arm.
#' @param n_groups Integer vector of group counts for single-variable
#'   quantile cuts. Default \code{c(3, 4, 5)}.
#' @param interaction_configs Pairwise quantile crossings, as strings such
#'   as \code{"2x3"} (first variable of the pair in 2 groups, second in 3)
#'   or as length-2 integer vectors. Default
#'   \code{c("2x2", "2x3", "3x2", "3x3")}. \code{NULL} suppresses
#'   quantile-by-quantile crossings.
#' @param breaks Optional named list of numeric cutpoint vectors, one
#'   element per variable (names must be in \code{vars}), e.g.
#'   \code{list(Yhat_raw = c(0, 1, 2, 3, 4, 5, 9))}. Give a named list of
#'   vectors to get several cuts of one variable, e.g.
#'   \code{list(Yhat_raw = list(binned = c(0, 1, 2, 3, 4, 5, 9),
#'   binned2 = c(0, 2, 3, 4, 9)))}. Cutpoints are fixed values applied
#'   identically in both arms.
#' @param custom Optional named list of externally built stratifications
#'   (for example, embedding-cluster labels). Each element is a vector of
#'   length \code{nrow(dat)} or a function of \code{dat} returning one.
#' @param cross_levels Logical; let native levels and \code{breaks} cuts
#'   enter pairwise crossings, with the quantile cuts of the other
#'   variable (at the group counts appearing in
#'   \code{interaction_configs}) and with each other. Default \code{TRUE}.
#' @param cross_custom Logical; cross each \code{custom} stratification
#'   with every single-variable stratification of the variables in
#'   \code{vars}. Default \code{TRUE}.
#' @param max_levels A variable that is a factor, character, or logical,
#'   or numeric with at most this many distinct values, is also offered at
#'   its native levels (\code{{var}_levels}). Default \code{10}.
#' @param min_cell_size Minimum number of units required in every cell of
#'   a candidate, where a cell is a stratum-by-arm combination when
#'   \code{Z} is supplied and a stratum otherwise. Candidates with any
#'   smaller cell, including a stratum that is empty in one arm, are
#'   dropped. Default \code{2}, the smallest cell with an estimable
#'   within-cell variance; \code{0} disables the screen.
#' @param drop_duplicates Logical; drop candidates that induce exactly the
#'   same partition of the units as an earlier candidate (for instance,
#'   the tertiles of a three-level variable duplicate its native levels).
#'   The first construction keeps its name; dropped ones are listed in the
#'   screening log. Default \code{TRUE}.
#' @param jitter_vars Character vector (subset of \code{vars}) to which
#'   small random noise is added before computing quantile breaks, to
#'   break ties. Useful for integer-valued variables such as word counts,
#'   especially with \code{cut_by_arm = TRUE}. Default \code{NULL}.
#' @param sd_jitter Standard deviation of the jitter noise. Default
#'   \code{1e-3}.
#' @param seed Optional integer seed for reproducibility when jitter is
#'   applied. Default \code{NULL}.
#' @param verbose Logical; print a one-line summary of how many candidates
#'   were constructed, kept, and dropped. Default \code{TRUE}.
#'
#' @return A named list of factors, each a candidate stratification of
#'   length \code{nrow(dat)}, with a \code{"screening"} attribute: a tibble
#'   with one row per constructed candidate giving \code{n_strata}, the
#'   smallest cell (\code{min_cell}), whether it was \code{kept}, and the
#'   \code{reason} if not. Retrieve it with
#'   \code{attr(candidates, "screening")}.
#'
#' @details
#' Naming. Single-variable candidates are \code{{var}_levels} (native
#' levels), \code{{var}_cut} or \code{{var}_{name}} (user cutpoints), and
#' \code{{var}_tertiles}, \code{{var}_quartiles}, \code{{var}_quintiles},
#' or \code{{var}_q{k}} (quantile cuts). Pairwise crossings are
#' \code{{var1}_{var2}_{tag1}x{tag2}}, where a tag is a group count,
#' \code{levels}, or the name of a cut (for example
#' \code{Yhat_word_count_2x3} or \code{Yhat_all_agree_3xlevels}). Custom
#' stratifications appear under their own names and, when crossed, as
#' \code{{custom}_x_{var}_{tag}}.
#'
#' Screening runs in construction order (native levels and user cuts,
#' then quantile cuts, then custom stratifications, then crossings that
#' involve levels or cuts, then quantile crossings), so when two
#' candidates coincide the more descriptive name survives. Candidates
#' whose construction fails (non-unique quantile breaks in the pooled
#' sample or in either arm, too few distinct values), that contain
#' missing values, or that have fewer than two strata are dropped and
#' logged as well.
#'
#' The \code{min_cell_size} screen is a necessary condition only: it
#' guarantees that within-cell variances exist, not that they are well
#' estimated. For pilot-based allocation, filter further on
#' \code{min_stratum_size} in the comparison tables.
#'
#' @examples
#' \dontrun{
#' # Default library, cut within arm and screened
#' cands <- create_stratification_candidates(dat, vars = c("Yhat", "word_count"), Z = "Z")
#' attr(cands, "screening")
#'
#' # READS-style library: rubric bins of the raw score, crossed with LLM
#' # agreement, plus an external clustering
#' cands <- create_stratification_candidates(
#'     dat, vars = c("Yhat", "Yhat_raw", "n_distinct_Yhat", "all_agree"), Z = "Z",
#'     breaks = list(Yhat_raw = list(binned = c(0, 1, 2, 3, 4, 5, 9),
#'                                   binned2 = c(0, 2, 3, 4, 9))),
#'     custom = list(embed_k6 = dat$embedding_cluster) )
#' }
#'
#' @export
create_stratification_candidates <- function(dat, vars, Z = NULL,
                                             cut_by_arm = TRUE,
                                             n_groups = c(3, 4, 5),
                                             interaction_configs = c("2x2", "2x3", "3x2", "3x3"),
                                             breaks = NULL,
                                             custom = NULL,
                                             cross_levels = TRUE,
                                             cross_custom = TRUE,
                                             max_levels = 10,
                                             min_cell_size = 2,
                                             drop_duplicates = TRUE,
                                             jitter_vars = NULL, sd_jitter = 1e-3, seed = NULL,
                                             verbose = TRUE) {
    
    if (!is.null(seed)) set.seed(seed)
    
    n <- nrow(dat)
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Argument checks (NEW)
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    stopifnot( is.numeric(min_cell_size), length(min_cell_size) == 1,
               !is.na(min_cell_size), min_cell_size >= 0 )
    
    missing_vars <- setdiff(vars, names(dat))
    if (length(missing_vars) > 0) {
        stop("Variable(s) not found in dat: ", paste(missing_vars, collapse = ", "))
    }
    
    n_groups <- sort(unique(as.integer(n_groups)))
    n_groups <- n_groups[!is.na(n_groups) & n_groups >= 2]
    
    # interaction_configs -> named list of c(n1, n2); "2x3" -> c(2, 3)
    configs <- interaction_configs
    if (is.character(configs)) {
        configs <- lapply(strsplit(configs, "x", fixed = TRUE), as.integer)
    }
    if (length(configs) > 0) {
        ok <- vapply(configs,
                     function(cfg) length(cfg) == 2 && !anyNA(cfg) && all(cfg >= 2),
                     logical(1))
        if (!all(ok)) {
            stop("interaction_configs must be strings like '2x3' or length-2 integer vectors with entries >= 2")
        }
        names(configs) <- vapply(configs, function(cfg) paste0(cfg[1], "x", cfg[2]), character(1))
    }
    cross_ns <- if (length(configs) > 0) sort(unique(unlist(configs))) else integer(0)
    
    if (!is.null(breaks)) {
        if (is.null(names(breaks)) || any(names(breaks) == "")) {
            stop("breaks must be a named list with names in vars")
        }
        bad <- setdiff(names(breaks), vars)
        if (length(bad) > 0) {
            stop("breaks supplied for variables not in vars: ", paste(bad, collapse = ", "))
        }
    }
    
    if (!is.null(custom)) {
        if (is.null(names(custom)) || any(names(custom) == "")) {
            stop("custom must be a named list")
        }
        custom <- lapply(custom, function(cs) {
            if (is.function(cs)) cs <- cs(dat)
            if (length(cs) != n) {
                stop("each custom stratification must have length nrow(dat)")
            }
            cs
        })
    }
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Helper functions
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    # Apply jitter to a variable (ORIGINAL)
    add_jitter <- function(x) x + rnorm(length(x), 0, sd_jitter)
    
    # Name tag for a k-group quantile cut (NEW)
    quant_tag <- function(k) {
        switch(as.character(k), "3" = "tertiles", "4" = "quartiles",
               "5" = "quintiles", paste0("q", k))
    }
    
    # Variables offered at their native levels (NEW)
    is_categorical <- function(x) {
        is.factor(x) || is.character(x) || is.logical(x) || n_distinct(x) <= max_levels
    }
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Validate and extract treatment indicator if provided (ORIGINAL)
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
    
    # Pooled or within-arm quantile cut (NEW)
    cut_one <- function(x, k, v) {
        if (is.null(z_vec) || !cut_by_arm) {
            make_quantile_strata(x, k, v)
        } else {
            make_quantile_strata_by_arm(x, z_vec, k, v)
        }
    }
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Prepare variables (with jitter for specified vars) (ORIGINAL)
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    if (is.null(jitter_vars)) jitter_vars <- character(0)
    
    var_data <- lapply(vars, function(v) {
        x <- dat[[v]]
        if (is.null(x)) stop(paste("Variable", v, "not found in dat"))
        if (v %in% jitter_vars) add_jitter(x) else x
    })
    names(var_data) <- vars
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Generate single-variable stratifications
    # single_parts[[v]] holds every single-variable partition of v by tag;
    # NULL entries are kept (via list(NULL)) so failures reach the log.
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    single_parts <- list()
    
    for (v in vars) {
        x <- var_data[[v]]
        parts <- list()
        
        # (a) NEW: native levels of categorical variables
        if (is_categorical(x)) {
            parts["levels"] <- list(make_level_strata(x, v))
        }
        
        # (b) NEW: user-supplied cutpoints
        if (!is.null(breaks[[v]])) {
            b <- breaks[[v]]
            if (!is.list(b)) b <- list(b)
            tags <- names(b)
            if (is.null(tags)) tags <- rep("", length(b))
            for (j in seq_along(b)) {
                tag <- if (tags[j] != "") tags[j] else if (length(b) == 1) "cut" else paste0("cut", j)
                parts[tag] <- list(make_break_strata(x, b[[j]], v))
            }
        }
        
        # (c) quantile cuts: the ORIGINAL tertiles/quartiles/quintiles,
        #     now for every k in n_groups and within arm when cut_by_arm
        for (k in n_groups) {
            parts[quant_tag(k)] <- list(cut_one(x, k, v))
        }
        
        single_parts[[v]] <- parts
    }
    
    single_strats <- list()
    for (v in vars) {
        for (tag in names(single_parts[[v]])) {
            single_strats[paste0(v, "_", tag)] <- list(single_parts[[v]][[tag]])
        }
    }
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Custom (externally built) stratifications, alone and crossed with
    # every single-variable stratification (NEW)
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    custom_strats <- list()
    
    if (!is.null(custom)) {
        for (cn in names(custom)) {
            cs <- make_level_strata(custom[[cn]], cn)
            custom_strats[cn] <- list(cs)
            if (cross_custom && !is.null(cs)) {
                for (v in vars) {
                    for (tag in names(single_parts[[v]])) {
                        custom_strats[paste0(cn, "_x_", v, "_", tag)] <-
                            list(cross_strata(cs, single_parts[[v]][[tag]]))
                    }
                }
            }
        }
    }
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Generate interaction stratifications for all pairs of variables
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    interaction_strats <- list()   # quantile x quantile (ORIGINAL loop)
    level_cross_strats <- list()   # crossings involving levels or cuts (NEW)
    
    if (length(vars) >= 2) {
        # All unique pairs
        var_pairs <- combn(vars, 2, simplify = FALSE)
        
        # NEW: crossings involving native levels or user cuts. Each such
        # partition of one variable is crossed with the quantile cuts of
        # the other variable at the group counts used in
        # interaction_configs, and with the other variable's levels/cuts.
        if (cross_levels) {
            quantile_tags <- vapply(n_groups, quant_tag, character(1))
            extra_parts <- lapply(vars, function(v) {
                p <- single_parts[[v]]
                p[setdiff(names(p), quantile_tags)]
            })
            names(extra_parts) <- vars
            cross_parts <- lapply(vars, function(v) {
                p <- list()
                for (k in cross_ns) p[as.character(k)] <- list(cut_one(var_data[[v]], k, v))
                p
            })
            names(cross_parts) <- vars
            
            for (pair in var_pairs) {
                v1 <- pair[1]
                v2 <- pair[2]
                e1 <- extra_parts[[v1]]
                e2 <- extra_parts[[v2]]
                q1 <- cross_parts[[v1]]
                q2 <- cross_parts[[v2]]
                others2 <- c(e2, q2)   # levels/cuts first so their names survive duplicates
                for (t1 in names(e1)) {
                    for (t2 in names(others2)) {
                        level_cross_strats[paste0(v1, "_", v2, "_", t1, "x", t2)] <-
                            list(cross_strata(e1[[t1]], others2[[t2]]))
                    }
                }
                for (t1 in names(q1)) {
                    for (t2 in names(e2)) {
                        level_cross_strats[paste0(v1, "_", v2, "_", t1, "x", t2)] <-
                            list(cross_strata(q1[[t1]], e2[[t2]]))
                    }
                }
            }
        }
        
        # ORIGINAL: quantile x quantile crossings, one per interaction config
        for (pair in var_pairs) {
            v1 <- pair[1]
            v2 <- pair[2]
            
            for (config_name in names(configs)) {
                n1 <- configs[[config_name]][1]
                n2 <- configs[[config_name]][2]
                
                strat_name <- paste0(v1, "_", v2, "_", config_name)
                s_int <- make_interaction_strata(
                    var_data[[v1]], n1, v1,
                    var_data[[v2]], n2, v2,
                    z_vec = z_vec, cut_by_arm = cut_by_arm
                )
                if ( length( s_int ) != n ) {
                    s_int <- NULL
                }
                interaction_strats[strat_name] <- list(s_int)
            }
        }
    }
    
    candidates <- c(single_strats, custom_strats, level_cross_strats, interaction_strats)
    
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    # Screen: failed construction, NA, < 2 strata, min_cell_size, duplicates (NEW)
    # -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
    
    scr <- screen_candidates(candidates, n = n, z_vec = z_vec,
                             min_cell_size = min_cell_size,
                             drop_duplicates = drop_duplicates)
    candidates <- scr$candidates
    
    if (verbose) {
        rs <- scr$screening$reason
        n_small <- sum(grepl("^smallest cell", rs))
        n_dup   <- sum(grepl("^duplicate", rs))
        n_other <- nrow(scr$screening) - length(candidates) - n_small - n_dup
        message(sprintf(
            "create_stratification_candidates: %d constructed, %d kept (%d dropped for a cell below %d units, %d duplicate, %d failed or invalid).",
            nrow(scr$screening), length(candidates), n_small, min_cell_size, n_dup, n_other))
    }
    
    attr(candidates, "screening") <- scr$screening
    
    return(candidates)
}



