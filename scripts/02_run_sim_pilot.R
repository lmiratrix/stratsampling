
# =====================================================================
# Pilot-based optimal allocation simulation, balanced-exact
# configuration only. '
#
# Each iteration computes, on the SAME population:
#   Stratified (proportional), Optimal (oracle Neyman),
#   and for each pilot budget fraction:
#     Optimal-Pilot-Pooled (pilot kept in estimation sample)
#     Optimal-Pilot-Ind  (pilot discarded; robustness arm)
#
# n_code is the TOTAL budget everywhere (pilot included), so all
# methods at a given coding fraction are directly comparable.

# Compared to 01_run_full_sim.R, this script focuses on performance of the pilot-based estimators
# across different pilot fractions (fraction of the total coding budget used for the pilot)
# =====================================================================

library(dplyr)
library(purrr)
library(stratsampling)
source("scripts/sim_functions.R")

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

#PILOT_FRACS <- c(1/5,1/4,1/3, 1/2,3/5,3/4)   # share of the coding budget spent on the pilot

PILOT_FRACS <- sort(c(seq(0.2,0.7,by=0.1),1/4,1/3,2/3,3/4))  # share of the coding budget spent on the pilot
PILOT_FRACS

n_iter      <- 1000             # set to e.g. 100 for a test run

param_grid <- expand.grid(
    N = 1000,
    coding_budget = seq(0.1, 0.9, by = 0.1),
    n_strat = 4,
    strata_config = "balanced-exact",
    tau = 0,
    bias = c("none", "small","moderate", "large","extreme-contrast"),
    residual_var = c("homogeneous", "heterogeneous","extreme-contrast"),
    sd_y = 3,
    R2_target = c(0.4, 0.85),
    stringsAsFactors = FALSE
)

cli::cli_inform("Total simulation scenarios: {nrow(param_grid)}\n")



# ---------------------------------------------------------------------
# One iteration: 3 baselines + 2 pilot allocation arms + 2 pilot selection arms per pilot fraction,
# all on the same population draw
# ---------------------------------------------------------------------


one_iter_pilot_study = function(params, pilot_fracs = PILOT_FRACS) {
    
    dat = generate_population(N = params$N, n_strat = params$n_strat,
                              strata_config = params$strata_config,
                              tau = params$tau,
                              bias = params$bias,
                              residual_var = params$residual_var,
                              sd_y = params$sd_y,
                              R2_target = params$R2_target)
    
    n_code = round(params$N * params$coding_budget)
    
    out.prop   = one_iter_strat(dat, n_code, Z = "Z");    out.prop$method   = "Stratified"
    out.opt    = one_iter_optimal_strat(dat, n_code, Z = "Z"); out.opt$method = "Optimal"
    out.opt_rand = one_iter_optimal_wrong(dat, n_code, Z= "Z", distort="random"); out.opt_rand$method = "Optimal-Random"
    out.opt_perm = one_iter_optimal_wrong(dat, n_code, Z= "Z", distort="permute"); out.opt_perm$method = "Optimal-Permuted"
    
    
    base = bind_rows( out.prop, out.opt, out.opt_rand,out.opt_perm) %>%
        mutate(pilot_frac = NA_real_)
    
    # Pilot arms: fresh pilot draw per (fraction x mode) on the same
    # population. suppressWarnings() silences the expected thin-cell
    # fallback messages at small pilots (verified once interactively).
    # Pilot arms: ONE proportional stratified pilot per fraction, shared
    # by the Pooled and Ind arms (paired comparison); fresh SRS pilot for
    # the selection arms. suppressWarnings() silences the expected
    # thin-cell fallback / below-floor messages at small pilots.
    pilots = purrr::map_dfr(pilot_fracs, function(pf) {
        n_pilot = round(pf * n_code)
        dat$S_pilot = generate_sample_index(dat, n_pilot, Z = "Z", strata = "strata")
        vt = suppressWarnings(
            estimate_pilot_variances(dat, S = "S_pilot", Z = "Z", strata = "strata"))
        
        # Pooled: wave 2 tops the pilot up to the total budget
        d_pool = suppressWarnings(
            calc_optimal_allocation(dat, n_code, Z = "Z", strata = "strata",
                                    var_alloc = vt, S_pilot = "S_pilot"))
        pool = est_ma_strat(d_pool, Yhat = "Yhat", S = "S", Z = "Z", strata = "strata")
        pool$method = "Optimal-Pilot-Pooled"
        
        # Ind: same pilot informs the variances only; fresh draw of the
        # remaining budget from the full population
        d_ind = suppressWarnings(
            calc_optimal_allocation(dat, n_code - sum(dat$S_pilot),
                                    Z = "Z", strata = "strata", var_alloc = vt))
        ind = est_ma_strat(d_ind, Yhat = "Yhat", S = "S", Z = "Z", strata = "strata")
        ind$method = "Optimal-Pilot-Ind"
        
        select = suppressWarnings(
            one_iter_pilot_select(dat, n_code, candidates = NULL,
                                  pilot_select_frac = pf, Z = "Z"))
        
        bind_rows(pool, ind, select) %>% mutate(pilot_frac = pf)
    })
    
    bind_rows(base, pilots) %>%
        relocate(method) %>%
        mutate(N = params$N,
               n_code = n_code,
               n_strat = params$n_strat,
               strata_config = params$strata_config,
               tau = params$tau,
               bias = params$bias,
               residual_var = params$residual_var,
               pop.sd = params$sd_y,
               R2_target = params$R2_target)%>%
        select(method, pilot_frac, everything())
}


if (FALSE) {
    # Test one iteration: expect 4 + 4*length(PILOT_FRACS)  rows
    one_iter_pilot_study(param_grid[1, ])
}

# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------

library(doParallel)
library(foreach)
library(doRNG)

cores <- min(parallel::detectCores() - 1, nrow(param_grid))
cl <- makeCluster(cores)
registerDoParallel(cl)

cli::cli_h1("Running the full pilot simulation" )

set.seed(1234)

system.time(
    results_pilot <- foreach(i = 1:nrow(param_grid),
                             .packages = c("dplyr", "purrr", "stratsampling")) %dorng% {
                                 plyr::rdply(n_iter, one_iter_pilot_study(param_grid[i, ]))
                             }
)

stopCluster(cl)

cli::cli_h1("Simulation complete.  Now saving" )

results_pilot <- dplyr::bind_rows(results_pilot)

saveRDS(results_pilot, file = "results/sim_results_pilot_balanced_exact.rds")
