# =====================================================================
# test_calibration_single_scenario.R
#
# Targeted diagnostic: run one_iter_all_fast() many times for a SINGLE
# simulation scenario and report Monte Carlo bias / coverage / SEhat
# calibration for every estimator. This mirrors the "validity6" table
# built in 03a_sim_results_figures.R (section 6d), but for one scenario
# at much higher n_iter, so a calibration problem is visible without
# waiting on the full factorial grid (432 scenarios x 9 budgets).
#
# Motivating observation (validity6 on the full grid): the two pooled
# "+ pilot estimate" arms -- Pilot+Stratified-PilotSelect and
# Pilot+Optimal-PilotSelect -- under-cover (~0.90 vs nominal 0.95) with
# SEhat/SE ~ 0.83, while every other arm (including the un-pooled
# Stratified-PilotSelect / Optimal-PilotSelect rows built from the same
# wave-2 sample) is essentially exactly calibrated.
#
# Likely cause, already flagged as a caveat in sim_functions.R: inside
# one_iter_pilot_select(), the pilot (wave 1) and the wave-2 sample are
# both drawn from the FULL population, independently of each other, so
# they can and do share coded units. combine_estimators() weights them
# as if they were independent draws; any positive covariance induced by
# that overlap makes its pooled SEhat anti-conservative (too small),
# which is exactly the under-coverage / SEhat-too-small pattern seen
# only in the two "+pilot estimate" rows.
#
# Usage: run the whole file top to bottom once (Part A simulates and
# saves; Part B loads and prints calibration). After that, Part B can be
# re-run by itself -- source only from the "Part B" marker down, or just
# re-run the whole file with FORCE_RERUN left FALSE.
# =====================================================================

library(tidyverse)

FAST_SRC <- normalizePath("scripts/sim_functions_fast.R")
source(FAST_SRC)
source("scripts/plotting_functions.R")   # method_levels / method_labels

RESULTS_FILE <- "results/test_calibration_single_scenario.rds"
FORCE_RERUN  <- TRUE

# cli::cli_h1() does not wrap -- a title longer than the console width is
# silently truncated (the rest is dropped, not just visually cut off), which
# lost several scenario parameters below. Wrap it ourselves: the first
# wrapped line becomes the ruled h1, any remaining lines print underneath.
cli_h1_wrap <- function(text) {
    w <- max(20, cli::console_width() - 4L)
    lines <- strwrap(text, width = w)
    cli::cli_h1(lines[1])
    if (length(lines) > 1) cli::cli_text(paste(lines[-1], collapse = " "))
}



# ---------------------------------------------------------------------
# Part A: run the simulation (skipped if RESULTS_FILE already exists ----
# and FORCE_RERUN is FALSE, so Part B can be iterated on without
# rerunning n_iter reps every time).
# ---------------------------------------------------------------------
if (FORCE_RERUN || !file.exists(RESULTS_FILE)) {
    
    n_iter       <- 3000
    PILOT_FRAC   <- 1 / 5     # matches 01_run_full_sim_fast.R
    INCLUDE_RAND <- TRUE
    PAIRED_PILOT <- TRUE
    
    # One representative scenario from the middle of the main grid.
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
    
    method_lab <- all_method_labels(INCLUDE_RAND, PILOT_FRAC)
    nr <- length(method_lab)
    nc <- length(FAST_COLS)
    
    cli::cli_inform("Number methods = {nr} and number FAST_COLS = {nc}" )
    
    if ( FALSE ) {
        one_iter_all_fast(params, pilot_frac = PILOT_FRAC,
                          include_rand = INCLUDE_RAND,
                          paired_pilot = PAIRED_PILOT)
    }
    
    set.seed(1234)
    M <- matrix(NA_real_, n_iter * nr, nc)
    for (it in seq_len(n_iter)) {
        M[((it - 1L) * nr + 1L):(it * nr), ] <-
            suppressWarnings(one_iter_all_fast(params, pilot_frac = PILOT_FRAC,
                                               include_rand = INCLUDE_RAND,
                                               paired_pilot = PAIRED_PILOT))
        if (it %% 500 == 0) cli::cli_inform("iter {it}/{n_iter}")
    }
    
    results <- data.frame(
        .n           = rep(seq_len(n_iter), each = nr),
        method       = rep(method_lab, times = n_iter),
        est          = M[, 1],
        SE           = M[, 2],
        SEhat        = M[, 3],
        SEhat.full   = M[, 4],
        SEhat.samp   = M[, 5],
        S2.0         = M[, 6],
        S2.1         = M[, 7],
        strat_choice = CAND_NAMES[M[, 8]],
        n_pilot      = M[, 9],
        stringsAsFactors = FALSE
    )
    
    dir.create("results", showWarnings = FALSE)
    saveRDS(list(results = results, params = params, pilot_frac = PILOT_FRAC), RESULTS_FILE)
    cli::cli_inform("Saved {n_iter} reps x {nr} methods to {RESULTS_FILE}")
}

# ---------------------------------------------------------------------
# Part B: load and analyze -- can be re-run on its own once Part A has ----
# produced RESULTS_FILE at least once.
# ---------------------------------------------------------------------
saved      <- readRDS(RESULTS_FILE)
results    <- saved$results
tau        <- saved$params$tau
PILOT_FRAC <- saved$pilot_frac   # Parts C/D reuse this so they always match Part A

cli_h1_wrap(sprintf(
    "Scenario: %s",
    paste(names(saved$params), unlist(saved$params), sep = "=", collapse = ", ")
))

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
        mean_SEhat = sqrt( mean(SEhat^2) ),          # matches analysis_results' mean_sehat in 03a
        `SEhat/SE` = mean_SEhat / emp_SE,
        coverage   = mean(covers),
        .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = intersect(method_levels, unique(method)))) %>%
    arrange(method) %>%
    mutate(Design = method_labels[as.character(method)]) %>%
    relocate(Design, .before = method) %>%
    ungroup()  %>%
    mutate( var = emp_SE^2,
            rel_eff    = emp_SE^2 / emp_SE[ method =="Stratified" ]^2 )


calib %>%
    dplyr::select( -n, -bias ) %>%
    knitr::kable( digits = 3 ) %>%
    print()


cli::cli_inform(c(
    "i" = "Watch the two rows 'Proportional/Optimal (pilot-selected strata) + pilot estimate'.",
    "i" = "If they alone under-cover with SEhat/SE < 1 while the un-pooled",
    "i" = "Stratified-PilotSelect/Optimal-PilotSelect rows (same wave-2 sample) are calibrated,",
    "i" = "that isolates the problem to combine_estimators()'s independence assumption in",
    "i" = "one_iter_pilot_select() (pilot and wave-2 samples are both drawn from the full",
    "i" = "population and can overlap)."
))







# =====================================================================
# Part C: pilot-selection diagnostics pt 1 ----
#
# Two follow-up questions on why the pilot-select/pooled rows' TRUE SE
# (not just SEhat) isn't much lower than SRS despite spending the same
# total budget:
#
#  (1) Which candidate stratification of Yhat (tertiles/quartiles/
#      quintiles) does the n_pilot = round(pilot_frac*n_code/2)*2 = 100
#      unit SRS pilot actually pick, iteration to iteration? A noisy
#      choice among the three candidates dilutes the gain relative to
#      always using the single best one.
#  (2) one_iter_pilot_select() draws the pilot (wave 1) and each wave-2
#      sample INDEPENDENTLY from the full population (see the
#      "Counterargument" note above one_iter_pilot_select() in
#      sim_functions.R): any unit landing in both waves is charged
#      against the budget twice even though it's coded only once. That
#      shrinks the realized (deduplicated) estimation sample below
#      n_pilot + n_main, on top of the fact that the un-pooled rows
#      already only use n_main = n_code - n_pilot = 200 (not the full
#      300) since the pilot is discarded from them entirely.
#
# strat_choice for (1) is already recorded in RESULTS_FILE (Part A/B)
# for every pilot-select row; this section adds (2), which needs the
# realized wave-1/wave-2 sample INDICATOR VECTORS (not just counts), so
# it re-derives one_iter_pilot_select()'s fast-path building blocks
# directly (srs_fast/gsi_fast/alloc_fast all already return 0/1 vectors
# over the same N rows in the same order -- see sim_functions_fast.R
# section 9 -- so overlap is just a vector AND/OR).
# ---------------------------------------------------------------------

cli_h1_wrap("Which stratification does the pilot select? (from Part A/B results)")
results %>%
    filter(method == "Optimal-PilotSelect", !is.na(strat_choice)) %>%
    count(strat_choice) %>%
    mutate(share = n / sum(n)) %>%
    print()



# Part D: Rerun sim with extra diagnostics on reuse of units ----

PILOT_DIAG_FILE <- "results/test_calibration_single_scenario_pilot_diag.rds"


if (FORCE_RERUN || !file.exists(PILOT_DIAG_FILE)) {

    n_iter <- 3000
    params <- saved$params            # same scenario as Part A
    n_code <- round(params$N * params$coding_budget)

    # Mirrors pilot_select_fast() (sim_functions_fast.R) but also
    # returns the realized wave-1/wave-2 overlap and union sizes.
    pilot_select_diag <- function(pc, cands, n_code, pilot_select_frac, pc1) {
        n_pilot <- round(pilot_select_frac * n_code / 2) * 2
        n_main  <- n_code - n_pilot
        
        Sp   <- srs_fast(pc, n_pilot)                 # wave 1 (pilot)
        best <- select_fast(cands, Sp)
        pcs  <- cands$cells[[best]]
        
        S_prop <- gsi_fast(pcs, n_main)                # wave 2, proportional
        vt     <- tryCatch(pilotvar_shrink_fast(pcs, Sp), error = function(e) NULL)
        S_opt  <- if (is.null(vt)) NULL else alloc_fast(pcs, n_main, vt)  # wave 2, Neyman
        
        overlap_prop <- sum(Sp == 1L & S_prop == 1L)
        overlap_opt  <- if (is.null(S_opt)) NA_integer_ else sum(Sp == 1L & S_opt == 1L)
        
        data.frame(
            strat_choice = CAND_NAMES[best],
            n_pilot      = sum(Sp),
            n_main_prop  = sum(S_prop),
            overlap_prop = overlap_prop,
            unique_prop  = sum(Sp) + sum(S_prop) - overlap_prop,
            n_main_opt   = if (is.null(S_opt)) NA_integer_ else sum(S_opt),
            overlap_opt  = overlap_opt,
            unique_opt   = if (is.null(S_opt)) NA_integer_ else sum(Sp) + sum(S_opt) - overlap_opt
        )
    }
    
    set.seed(5678)
    diag_rows <- vector("list", n_iter)
    for (it in seq_len(n_iter)) {
        pc    <- pop_fast(params)
        cands <- candidates_fast(pc)
        pc1   <- one_stratum_fast(pc)
        diag_rows[[it]] <- pilot_select_diag(pc, cands, n_code, PILOT_FRAC, pc1)
        if (it %% 500 == 0) cli::cli_inform("pilot-diag iter {it}/{n_iter}")
    }
    pilot_diag <- bind_rows(diag_rows)
    
    dir.create("results", showWarnings = FALSE)
    saveRDS(list(pilot_diag = pilot_diag, params = params, n_code = n_code,
                 pilot_frac = PILOT_FRAC),
            PILOT_DIAG_FILE)
    cli::cli_inform("Saved pilot-selection diagnostics to {PILOT_DIAG_FILE}")
}

pd_saved   <- readRDS(PILOT_DIAG_FILE)
pilot_diag <- pd_saved$pilot_diag
n_code     <- pd_saved$n_code

cli_h1_wrap("Pilot-selection choice, rederived directly (cross-check on the count above)")
print(pilot_diag %>% count(strat_choice) %>% mutate(share = n / sum(n)))

cli_h1_wrap("Effective (deduplicated) coded sample vs nominal budget")
print(pilot_diag %>%
          summarise(
              nominal_budget    = n_code,
              mean_n_pilot      = mean(n_pilot),
              mean_n_main_prop  = mean(n_main_prop),
              mean_overlap_prop = mean(overlap_prop),
              mean_unique_prop  = mean(unique_prop),
              pct_budget_wasted_prop = 100 * mean(overlap_prop) / n_code,
              mean_n_main_opt   = mean(n_main_opt, na.rm = TRUE),
              mean_overlap_opt  = mean(overlap_opt, na.rm = TRUE),
              mean_unique_opt   = mean(unique_opt, na.rm = TRUE),
              pct_budget_wasted_opt = 100 * mean(overlap_opt, na.rm = TRUE) / n_code
          ))

cli::cli_inform(c(
    "i" = "'pct_budget_wasted' = the share of the total coding budget spent re-coding a unit",
    "i" = "that the pilot had already coded (charged twice in the sim, coded once in practice).",
    "i" = "'mean_unique_*' vs 'nominal_budget' shows how much SMALLER the truly-informative",
    "i" = "deduplicated sample is than the nominal n_code, on top of the un-pooled rows already",
    "i" = "only using n_main = n_code - n_pilot for estimation (the pilot itself is discarded there)."
))

# =====================================================================
# Part E: is Stratified-PilotSelect's gap vs Stratified just "less ----
# budget," or is stratifying on QUANTILES OF Yhat inherently weaker
# than knowing the population's TRUE (DGP) strata?
#
# Every oracle row in Part A/B ("Stratified", "Optimal", ...) stratifies
# on the population's actual `strata` label -- the same label
# generate_surrogate_noise() used to assign each unit's bias and
# residual-variance multiplier. Stratified-PilotSelect/Optimal-PilotSelect
# never see that label at all: they stratify on quantile bins of Yhat,
# chosen by the pilot from {tertiles, quartiles, quintiles}. That is a
# different (and, especially under bias/residual_var = "extreme-contrast",
# considerably noisier) proxy for the true strata, so rescaling the
# oracle Stratified row's variance by n_code/n_main is not the right
# benchmark for what Stratified-PilotSelect "should" achieve. This
# isolates the two effects, with NO pilot and NO selection noise in
# either one:
#   (D1) True strata at n_main = n_code - n_pilot (same budget as the
#        un-pooled pilot-select rows, but the TRUE strata instead of a
#        Yhat-quantile proxy) -- tests "is it just less budget?"
#   (D2) Each of the 3 Yhat-quantile candidates on their own, at the
#        FULL n_code budget AND at n_main -- tests "how much does
#        proxying strata with Yhat quantiles cost, by itself?" The
#        n_main versions are the closest apples-to-apples comparison to
#        Stratified-PilotSelect (same budget, same candidate variable),
#        minus the pilot's noisy choice among the three.
# ---------------------------------------------------------------------
STRAT_COMPARE_FILE <- "results/test_calibration_single_scenario_strat_compare.rds"

if (FORCE_RERUN || !file.exists(STRAT_COMPARE_FILE)) {

    n_iter     <- 3000
    params     <- saved$params        # PILOT_FRAC comes from Part B (== Part A's value)
    n_code     <- round(params$N * params$coding_budget)
    n_pilot    <- round(PILOT_FRAC * n_code / 2) * 2
    n_main     <- n_code - n_pilot

    arm_labels <- c("True strata (full budget)",
                    "True strata (n_main budget)",
                    paste0(CAND_NAMES, " (full budget)"),
                    paste0(CAND_NAMES, " (n_main budget)"))
    n_arms <- length(arm_labels)

    set.seed(4242)
    M <- matrix(NA_real_, n_iter * n_arms, 2)   # est, SEhat
    for (it in seq_len(n_iter)) {
        pc    <- pop_fast(params)
        cands <- candidates_fast(pc)

        row0 <- (it - 1L) * n_arms
        M[row0 + 1L, ] <- est_fast(pc, gsi_fast(pc, n_code))[c(1, 3)]
        M[row0 + 2L, ] <- est_fast(pc, gsi_fast(pc, n_main))[c(1, 3)]
        for (j in 1:3) {
            pcs <- cands$cells[[j]]
            M[row0 + 2L + j, ]      <- est_fast(pcs, gsi_fast(pcs, n_code))[c(1, 3)]
            M[row0 + 5L + j, ]      <- est_fast(pcs, gsi_fast(pcs, n_main))[c(1, 3)]
        }
        if (it %% 500 == 0) cli::cli_inform("strat-compare iter {it}/{n_iter}")
    }

    strat_compare <- data.frame(
        arm   = rep(arm_labels, times = n_iter),
        est   = M[, 1],
        SEhat = M[, 2]
    )

    dir.create("results", showWarnings = FALSE)
    saveRDS(list(strat_compare = strat_compare, params = params,
                n_code = n_code, n_pilot = n_pilot, n_main = n_main),
            STRAT_COMPARE_FILE)
    cli::cli_inform("Saved {n_iter} reps x {n_arms} arms to {STRAT_COMPARE_FILE}")
}

sc_saved      <- readRDS(STRAT_COMPARE_FILE)
strat_compare <- sc_saved$strat_compare

arm_order <- c("True strata (full budget)", "True strata (n_main budget)",
              paste0(CAND_NAMES, " (full budget)"),
              paste0(CAND_NAMES, " (n_main budget)"))

cli_h1_wrap("True strata vs Yhat-quantile proxies -- full budget vs n_main, no pilot/selection noise")
sc_tab <- strat_compare %>%
    group_by(arm) %>%
    summarise(emp_SE = sd(est), mean_SEhat = mean(SEhat), .groups = "drop") %>%
    mutate(arm = factor(arm, levels = arm_order)) %>%
    arrange(arm) %>%
    mutate(rel_to_true_full = emp_SE / emp_SE[arm == "True strata (full budget)"])

print(sc_tab, n = Inf)

cli::cli_inform(c(
    "i" = "n_code = {sc_saved$n_code}, n_pilot = {sc_saved$n_pilot}, n_main = {sc_saved$n_main}.",
    "i" = "Row 2 vs row 1 isolates the pure 'less budget, same (true) strata' effect.",
    "i" = "Rows 3-5 (full budget) isolate the pure 'Yhat-quantile proxy instead of true strata'",
    "i" = "effect, with no budget loss at all.",
    "i" = "Rows 6-8 (n_main budget) are the closest apples-to-apples comparison to the",
    "i" = "Stratified-PilotSelect row in the Part A/B calibration table above -- same budget,",
    "i" = "same candidate stratification variable, just without the pilot's noisy choice among",
    "i" = "the three candidates."
))
