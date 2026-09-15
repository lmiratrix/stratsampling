# =====================================================================
# 01_run_full_sim_fast.R
#
# FAST runner for the main factorial simulation: identical design and
# identical results to 01_run_full_sim.R (same RNG stream per scenario,
# verified with verify_all_fast() in sim_functions_fast.R), ~50x faster
# per iteration.  Runs the three strata configurations and the
# full-coding oracle, writing the same four results files.
#
# Arms per iteration (one_iter_all order): Subset, SRS, Stratified,
# Optimal, Stratified-Random, Optimal-Permuted, Optimal-Reversed,
# Optimal-Random, and, when PILOT_FRAC is set, Optimal-Pilot-Pooled,
# Optimal-Pilot-Ind, Stratified-PilotSelect, Optimal-PilotSelect,
# Optimal-PilotSelect-Shrink (13 rows; 8 with PILOT_FRAC <- NULL).
#
# Run from the stratsampling/ root.
# =====================================================================
library(foreach)
library(doParallel)
library(doRNG)

FAST_SRC <- normalizePath("scripts/sim_functions_fast.R")
source(FAST_SRC)

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------
n_iter       <- 1000
PILOT_FRAC   <- 1/3      # share of the coding budget spent on the pilot arms; NULL = no pilot arms
INCLUDE_RAND <- TRUE
# Pilot draw for the Ind arm: FALSE = fresh pilot (one_iter_all() as
# written, which calls one_iter_optimal_pilot() twice); TRUE = one pilot
# shared with the Pooled arm.  Must match sim_functions.R.
PAIRED_PILOT <- TRUE

param_grid <- expand.grid(
    N = 1000,
    coding_budget = seq(0.1, 0.9, by = 0.1),
    n_strat = c(4),
    strata_config = "balanced-exact",
    tau = 0,
    bias = c("none", "small", "moderate", "large", "extreme-contrast"),
    residual_var = c("homogeneous", "heterogeneous", "extreme-contrast"),
    sd_y = c(3),
    R2_target = c(0.4, 0.85),
    stringsAsFactors = FALSE
)
cli::cli_inform("Total simulation scenarios: {nrow(param_grid)}\n")

method_lab <- all_method_labels(INCLUDE_RAND, PILOT_FRAC)
nr <- length(method_lab)
nc <- length(FAST_COLS)

cores <- min(parallel::detectCores() - 1, nrow(param_grid))
cl <- makeCluster(cores)
registerDoParallel(cl)
invisible(parallel::clusterCall(cl, function(p) source(p), FAST_SRC))   # workers need the helpers

# ---------------------------------------------------------------------
# One configuration: one task per scenario returning a plain numeric
# matrix (n_iter * nr) x nc; assembled into the long results table once
# at the end, with the same columns as 01_run_full_sim.R produces.
# ---------------------------------------------------------------------
run_config <- function(pg) {
    set.seed(1234)                                   # same seed as the reference runner
    tm <- system.time(
        res <- foreach(i = 1:nrow(pg),
                       .export = c("n_iter", "nr", "nc", "PILOT_FRAC",
                                   "INCLUDE_RAND", "PAIRED_PILOT"),.options.RNG=1234) %dorng% {
                                       p <- pg[i, ]
                                       M <- matrix(NA_real_, n_iter * nr, nc)
                                       for (it in seq_len(n_iter)) {
                                           M[((it - 1L) * nr + 1L):(it * nr), ] <-
                                               suppressWarnings(one_iter_all_fast(p, pilot_frac = PILOT_FRAC,
                                                                                  include_rand = INCLUDE_RAND,
                                                                                  paired_pilot = PAIRED_PILOT))
                                       }
                                       M
                                   }
    )
    print(tm)
    big <- do.call(rbind, res)
    n_s <- nrow(pg); rep_s <- n_iter * nr
    out <- data.frame(
        .n            = rep(rep(seq_len(n_iter), each = nr), times = n_s),
        method        = rep(method_lab, times = n_iter * n_s),
        est           = big[, 1],
        SE            = big[, 2],
        SEhat         = big[, 3],
        SEhat.full    = big[, 4],
        SEhat.samp    = big[, 5],
        S2.0          = big[, 6],
        S2.1          = big[, 7],
        N             = rep(pg$N, each = rep_s),
        n_code        = rep(round(pg$N * pg$coding_budget), each = rep_s),
        n_strat       = rep(pg$n_strat, each = rep_s),
        strata_config = rep(pg$strata_config, each = rep_s),
        tau           = rep(pg$tau, each = rep_s),
        bias          = rep(pg$bias, each = rep_s),
        residual_var  = rep(pg$residual_var, each = rep_s),
        pop.sd        = rep(pg$sd_y, each = rep_s),
        R2_target     = rep(pg$R2_target, each = rep_s),
        stringsAsFactors = FALSE, row.names = NULL
    )
    if (!is.null(PILOT_FRAC)) {
        out$pilot_frac   <- PILOT_FRAC
        out$strat_choice <- CAND_NAMES[big[, 8]]       # NA except pilot-select rows
        out$n_pilot      <- big[, 9]                   # NA for the non-pilot arms
    }
    out
}

cli::cli_h1("balanced-exact")
results <- run_config(param_grid)
saveRDS(results, file = "results/sim_results_balanced_exact_all_fast.rds")
rm(results)

cli::cli_h1("balanced-approx")
param_grid$strata_config <- "balanced-approx"
results2 <- run_config(param_grid)
saveRDS(results2, file = "results/sim_results_balanced_approx_all_fast.rds")
rm(results2)

cli::cli_h1("unbalanced")
param_grid$strata_config <- "unbalanced"
results3 <- run_config(param_grid)
saveRDS(results3, file = "results/sim_results_unbalanced_all_fast.rds")
rm(results3)

# ---------------------------------------------------------------------
# Full-coding oracle (difference in arm means on the whole population)
# ---------------------------------------------------------------------
cli::cli_h1("oracle (full coding)")
params_oracle <- expand.grid(
    N = 1000, n_strat = 4,
    strata_config = c("balanced-exact", "balanced-approx", "unbalanced"),
    tau = 0,
    bias = c("none", "small", "moderate", "large", "extreme-contrast"),
    residual_var = c("homogeneous", "heterogeneous", "extreme-contrast"),
    sd_y = 3,
    R2_target = unique(param_grid$R2_target),
    stringsAsFactors = FALSE
)

set.seed(1234)
res_o <- foreach(i = 1:nrow(params_oracle), .options.RNG=1234) %dorng% {
    p <- params_oracle[i, ]
    M <- matrix(NA_real_, n_iter, 2L)
    for (it in seq_len(n_iter)) M[it, ] <- oracle_fast(p)
    M
}
stopCluster(cl)

big_o <- do.call(rbind, res_o)
n_o <- nrow(params_oracle)
results.oracle <- data.frame(
    .n            = rep(seq_len(n_iter), times = n_o),
    method        = "Oracle",
    est           = big_o[, 1],
    SE            = big_o[, 2],
    N             = rep(params_oracle$N, each = n_iter),
    n_strat       = rep(params_oracle$n_strat, each = n_iter),
    strata_config = rep(params_oracle$strata_config, each = n_iter),
    tau           = rep(params_oracle$tau, each = n_iter),
    bias          = rep(params_oracle$bias, each = n_iter),
    residual_var  = rep(params_oracle$residual_var, each = n_iter),
    pop.sd        = rep(params_oracle$sd_y, each = n_iter),
    R2_target     = rep(params_oracle$R2_target, each = n_iter),
    stringsAsFactors = FALSE, row.names = NULL
)
saveRDS(results.oracle, file = "results/sim_results_oracle_fast.rds")
cli::cli_h1("done")
