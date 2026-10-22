# =====================================================================
# test_calibration_single_scenario_v2.R
#
# v2 of test_calibration_single_scenario.R (kept as-is; this is a
# separate, narrower script, not a replacement). The v1 table mixes
# several different "what are we stratifying on" stories together (the
# true DGP strata for every oracle row, Yhat-quantile candidates for the
# pilot-select rows) plus an allocation dimension (proportional vs
# Neyman) on top of that -- which is exactly what made the earlier
# results confusing to read. This script drops allocation ENTIRELY
# (every arm below is PROPORTIONAL; Neyman/Optimal stays in v1) and
# isolates just the "what variable are we stratifying on, and how much
# budget do we spend" story:
#
#   SRS                                    -- no stratification at all
#   True strata (oracle, full budget)      -- upper bound: real DGP strata, known
#   Strat on Yhat, 3/4/5 bins (full budget)-- Yhat-quantile proxy, varying bin count
#   Strat on noise X1, 4 bins (full budget)-- negative control: stratify on a
#                                              variable with NO relationship to
#                                              Yobs/Yhat/strata
#   Strat on Yhat, 4 bins (80% budget)     -- pure "less budget" effect on the
#                                              SAME stratifying variable as above
#   Pilot selects among 3 Yhat bins, 2nd stage only (pilot discarded)
#   Pilot selects + pilot/2nd-stage precision-combined
#
# Uses the REAL (non-fast) package functions throughout --
# create_stratification_candidates() needs a genuine extra covariate
# (X1) that the "fast" backend (sim_functions_fast.R) has no concept
# of, so there is no fast-path shortcut available here; n_iter is kept
# moderate accordingly.
#
# Same two-part save/load structure as v1: Part A simulates and saves;
# Part B loads and prints the calibration table, re-runnable on its own.
# =====================================================================

library(tidyverse)

source("scripts/sim_functions.R")   # one_iter_srs/one_iter_strat/one_iter_pilot_select/...

cli_h1_wrap <- function(text) {
    w <- max(20, cli::console_width() - 4L)
    lines <- strwrap(text, width = w)
    cli::cli_h1(lines[1])
    if (length(lines) > 1) cli::cli_text(paste(lines[-1], collapse = " "))
}

RESULTS_FILE <- "results/test_calibration_single_scenario_v2.rds"
FORCE_RERUN  <- FALSE

# Generic version of one_iter_strat() (sim_functions.R hardcodes
# strata = "strata"): proportional allocation on an ARBITRARY strata column.
one_iter_strat_on <- function(dat, n_code, strata_col, Z = "Z", finite_pop = FALSE) {
    dat$S <- generate_sample_index(dat, n_code, Z = Z, strata = strata_col)
    est_ma_strat(dat, Yhat = "Yhat", S = "S", Z = Z, strata = strata_col, finite_pop = finite_pop)
}

one_iter_strat_only <- function(dat, n_code, pilot_frac, Z = "Z", finite_pop = FALSE) {
    require(stratsampling)

    dat$X1 <- rnorm(nrow(dat))    # pure noise covariate: unrelated to Yobs/Yhat/strata

    cands_yhat <- create_stratification_candidates(dat, vars = "Yhat", Z = Z,
                                                    n_groups = c(3, 4, 5), verbose = FALSE)
    cands_x1   <- create_stratification_candidates(dat, vars = "X1", Z = Z,
                                                    n_groups = 4, verbose = FALSE)

    out_srs <- one_iter_srs(dat, n_code, finite_pop = finite_pop)
    out_srs$method <- "SRS"

    out_true <- one_iter_strat(dat, n_code, Z = Z, finite_pop = finite_pop)
    out_true$method <- "True strata (oracle)"

    out_yhat <- purrr::map_dfr(names(cands_yhat), function(nm) {
        dd <- dat
        dd$cand_strata <- cands_yhat[[nm]]
        r <- one_iter_strat_on(dd, n_code, strata_col = "cand_strata", Z = Z, finite_pop = finite_pop)
        r$method <- paste0("Strat-", nm)
        r
    })

    dat_x1 <- dat
    dat_x1$cand_strata <- cands_x1[["X1_quartiles"]]
    out_x1 <- one_iter_strat_on(dat_x1, n_code, strata_col = "cand_strata", Z = Z, finite_pop = finite_pop)
    out_x1$method <- "Strat-X1_quartiles (noise)"

    n_code_80 <- round(0.8 * n_code)
    dat80 <- dat
    dat80$cand_strata <- cands_yhat[["Yhat_quartiles"]]
    out_80 <- one_iter_strat_on(dat80, n_code_80, strata_col = "cand_strata", Z = Z, finite_pop = finite_pop)
    out_80$method <- "Strat-Yhat_quartiles (80% budget)"

    out_pilot <- suppressWarnings(one_iter_pilot_select(
        dat, n_code, candidates = cands_yhat, pilot_select_frac = pilot_frac,
        Z = Z, finite_pop = finite_pop)) %>%
        filter(method %in% c("Stratified-PilotSelect", "Pilot+Stratified-PilotSelect"))

    bind_rows(out_srs, out_true, out_yhat, out_x1, out_80, out_pilot) %>%
        relocate(method)
}

# ---------------------------------------------------------------------
# Part A: run the simulation (skipped if RESULTS_FILE already exists
# and FORCE_RERUN is FALSE).
# ---------------------------------------------------------------------
if (FORCE_RERUN || !file.exists(RESULTS_FILE)) {

    n_iter     <- 600      # ~0.5 sec/iter with the real (non-fast) package functions;
                            # bump up for tighter MC error if you can wait longer
    PILOT_FRAC <- 1 / 5     # matches v1's current setting

    # Same scenario as v1, for continuity.
    params <- list(
        N             = 1000,
        coding_budget = 0.3,
        n_strat       = 4,
        strata_config = "balanced-exact",
        tau           = 0,
        bias          = "extreme-contrast",
        residual_var  = "extreme-contrast",
        sd_y          = 3,
        R2_target     = 0.55
    )
    n_code <- round(params$N * params$coding_budget)

    set.seed(2026)
    res_list <- vector("list", n_iter)
    for (it in seq_len(n_iter)) {
        dat <- generate_population(N = params$N, n_strat = params$n_strat,
                                   strata_config = params$strata_config,
                                   tau = params$tau, bias = params$bias,
                                   residual_var = params$residual_var,
                                   sd_y = params$sd_y, R2_target = params$R2_target)
        res_list[[it]] <- one_iter_strat_only(dat, n_code, pilot_frac = PILOT_FRAC) %>%
            mutate(.n = it)
        if (it %% 250 == 0) cli::cli_inform("iter {it}/{n_iter}")
    }
    results <- bind_rows(res_list)

    dir.create("results", showWarnings = FALSE)
    saveRDS(list(results = results, params = params, pilot_frac = PILOT_FRAC, n_code = n_code),
            RESULTS_FILE)
    cli::cli_inform("Saved {n_iter} reps x {n_distinct(results$method)} arms to {RESULTS_FILE}")
}

# ---------------------------------------------------------------------
# Part B: load and analyze -- can be re-run on its own once Part A has
# produced RESULTS_FILE at least once.
# ---------------------------------------------------------------------
saved   <- readRDS(RESULTS_FILE)
results <- saved$results
tau     <- saved$params$tau

cli_h1_wrap(sprintf(
    "Scenario: %s (pilot_frac = %s, n_code = %s)",
    paste(names(saved$params), unlist(saved$params), sep = "=", collapse = ", "),
    saved$pilot_frac, saved$n_code
))

arm_order <- c("SRS", "True strata (oracle)",
              "Strat-Yhat_tertiles", "Strat-Yhat_quartiles", "Strat-Yhat_quintiles",
              "Strat-X1_quartiles (noise)", "Strat-Yhat_quartiles (80% budget)",
              "Stratified-PilotSelect", "Pilot+Stratified-PilotSelect")

nice_labels <- c(
    "SRS"                               = "SRS (no stratification)",
    "True strata (oracle)"              = "True strata (oracle, full budget)",
    "Strat-Yhat_tertiles"               = "Strat on Yhat, 3 bins (full budget)",
    "Strat-Yhat_quartiles"              = "Strat on Yhat, 4 bins (full budget)",
    "Strat-Yhat_quintiles"              = "Strat on Yhat, 5 bins (full budget)",
    "Strat-X1_quartiles (noise)"        = "Strat on noise X1, 4 bins (full budget)",
    "Strat-Yhat_quartiles (80% budget)" = "Strat on Yhat, 4 bins (80% budget)",
    "Stratified-PilotSelect"            = "Pilot selects among 3 Yhat bins, 2nd stage only",
    "Pilot+Stratified-PilotSelect"      = "Pilot selects + pilot/2nd-stage combined"
)

calib <- results %>%
    filter(!is.na(est), !is.na(SEhat)) %>%
    mutate(
        CI_lower = est - 2 * SEhat,
        CI_upper = est + 2 * SEhat,
        covers   = as.numeric(CI_lower < tau & CI_upper > tau)
    ) %>%
    group_by(method) %>%
    summarise(
        n          = n(),
        bias       = mean(est - tau),
        emp_SE     = sd(est),
        mean_SEhat = sqrt(mean(SEhat^2)),
        `SEhat/SE` = mean_SEhat / emp_SE,
        coverage   = mean(covers),
        .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = arm_order)) %>%
    arrange(method) %>%
    mutate(Design = nice_labels[as.character(method)]) %>%
    relocate(Design, .before = method) %>%
    mutate(rel_eff = emp_SE / emp_SE[method == "SRS"])

calib %>%
    dplyr::select(-n, -bias) %>%
    knitr::kable(digits = 3) %>%
    print()
