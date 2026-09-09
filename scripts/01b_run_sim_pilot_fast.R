# =====================================================================
# Pilot-based optimal allocation simulation, balanced-exact
# configuration only.  FAST runner: identical design and identical
# results to 01b_run_sim_pilot.R (same RNG stream), ~100x faster.
#
# Each iteration computes, on the SAME population:
#   Stratified (proportional), Optimal (oracle Neyman),
#   Optimal-Random (misspecified variances), and for each pilot budget
#   fraction:
#     Optimal-Pilot-Pooled (pilot kept in estimation sample)
#     Optimal-Pilot-Ind    (pilot discarded; robustness arm)
#
# n_code is the TOTAL budget everywhere (pilot included), so all
# methods at a given coding fraction are directly comparable.
# =====================================================================
library(foreach)
library(doParallel)
library(doRNG)

FAST_SRC <- normalizePath("scripts/sim_functions_fast.R")
source(FAST_SRC)

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------
PILOT_FRACS <- sort(c(seq(0.1, 0.7, by = 0.1), 1/4, 1/3, 2/3, 3/4))
n_iter      <- 1000

param_grid <- expand.grid(
    N = 1000,
    coding_budget = seq(0.1, 0.9, by = 0.1),
    n_strat = 4,
    strata_config = "balanced-exact",
    tau = 0,
    bias = c("none", "small", "moderate", "large", "extreme-contrast"),
    residual_var = c("homogeneous", "heterogeneous", "extreme-contrast"),
    sd_y = 3,
    R2_target = c(0.4, 0.85),
    stringsAsFactors = FALSE
)
cli::cli_inform("Total simulation scenarios: {nrow(param_grid)}\n")

# Row labels: constant within every iteration, so they are built once
# here instead of being pasted onto 25 one-row data frames per iteration.
nr <- 3L + 2L * length(PILOT_FRACS)
method_lab <- c("Stratified", "Optimal", "Optimal-Random",
                rep(c("Optimal-Pilot-Pooled", "Optimal-Pilot-Ind"),
                    times = length(PILOT_FRACS)))
pilot_lab  <- c(NA_real_, NA_real_, NA_real_, rep(PILOT_FRACS, each = 2))

# ---------------------------------------------------------------------
# Run: one task per scenario, each returning a plain numeric matrix
# (n_iter * nr) x 7.  No data frames, no bind_rows, inside the workers.
# ---------------------------------------------------------------------
cores <- min(parallel::detectCores() - 1, nrow(param_grid))
cl <- makeCluster(cores)
registerDoParallel(cl)
parallel::clusterCall(cl, function(p) source(p), FAST_SRC)   # workers need the helpers

set.seed(1234)
system.time(
    res <- foreach(i = 1:nrow(param_grid)) %dorng% {
        pg <- param_grid[i, ]
        M <- matrix(NA_real_, n_iter * nr, 7L)
        for (it in seq_len(n_iter)) {
            M[((it - 1L) * nr + 1L):(it * nr), ] <- one_iter_pilot_fast(pg, PILOT_FRACS)
        }
        M
    }
)
stopCluster(cl)

# ---------------------------------------------------------------------
# Assemble the long results table once, at the end
# ---------------------------------------------------------------------
big <- do.call(rbind, res)
rm(res)
n_s  <- nrow(param_grid)
rep_s <- n_iter * nr          # rows contributed by each scenario

results_pilot <- data.frame(
    .n            = rep(rep(seq_len(n_iter), each = nr), times = n_s),
    method        = rep(method_lab, times = n_iter * n_s),
    pilot_frac    = rep(pilot_lab,  times = n_iter * n_s),
    est           = big[, 1],
    SE            = big[, 2],
    SEhat         = big[, 3],
    SEhat.full    = big[, 4],
    SEhat.samp    = big[, 5],
    S2.0          = big[, 6],
    S2.1          = big[, 7],
    N             = rep(param_grid$N, each = rep_s),
    n_code        = rep(round(param_grid$N * param_grid$coding_budget), each = rep_s),
    n_strat       = rep(param_grid$n_strat, each = rep_s),
    strata_config = rep(param_grid$strata_config, each = rep_s),
    tau           = rep(param_grid$tau, each = rep_s),
    bias          = rep(param_grid$bias, each = rep_s),
    residual_var  = rep(param_grid$residual_var, each = rep_s),
    pop.sd        = rep(param_grid$sd_y, each = rep_s),
    R2_target     = rep(param_grid$R2_target, each = rep_s),
    stringsAsFactors = FALSE, row.names = NULL
)

saveRDS(results_pilot, file = "results/sim_results_pilot_balanced_exact_fast.rds")
