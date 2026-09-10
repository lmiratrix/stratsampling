# `scripts/`

This directory holds the research pipeline behind the paper (Mozer, Pashley, & Miratrix, 2026), as opposed to the installable `stratsampling` package itself (`R/`). Nothing here is exported, tested, or loaded by `library(stratsampling)` — these scripts *use* the package to run the paper's simulations and applied case studies. Numeric prefixes indicate execution order within each pipeline stage.

Most scripts read from or write to gitignored local directories (`raw data/`, `generated data/`, `results/`, `figures/`) that are not part of the repo. You will need those directories populated locally (or the raw source files, e.g. from Dropbox) before a script will run.

## Simulation pipeline (Monte Carlo comparisons of sampling designs)

- **`sim_functions.R`** — shared helpers sourced by the runner scripts below: builds one simulated population/sample/estimate per scenario given a parameter row.
- **`sim_functions_fast.R`** — drop-in base-R replacements for the same inner loop, used only by the `01b`/pilot runners. Deliberately consumes the RNG in the same order as the package/`sim_functions.R` versions so results agree to \~1e-12 for a given seed, at roughly 100x the speed.
- **`01_run_full_sim.R`** — main factorial simulation: builds a `param_grid` (coding budget, number of strata, bias pattern, residual-variance heterogeneity, target R², etc.) and runs stratified vs. SRS designs across it.
- **`02_run_sim_mean.R`** — companion simulation focused on mean-outcome (rather than treatment-effect) estimation.
- **`01b_run_sim_pilot.R`** / **`01b_run_sim_pilot_fast.R`** — pilot-based optimal (Neyman) allocation simulation: for each pilot budget fraction, compares proportional stratified allocation, oracle-optimal allocation, and pilot-estimated optimal allocation (pilot pooled into the estimation sample vs. discarded). `_fast` is the `sim_functions_fast.R`-powered version and should be preferred for actually running the grid; `01b_run_sim_pilot.R` is the reference/slow implementation used to verify `_fast` agrees with it.
- **`03a_sim_results_figures.R`**, **`03b_sim_performance_by_method.R`** — turn `01`/`02` simulation output into the paper's main comparison figures (MDES/frontier plots, performance-by-method breakdowns).
- **`03c_sim_pilot_results.R`** — turns `01b`/pilot simulation output (`results/sim_results_pilot_*.rds`) into tables/figures comparing pilot-based allocation methods.
- **`test_run_sim_mean_single_scenario.R`**, **`test_simple_calibration_simulation.R`** — ad hoc, interactive debugging scripts (not part of `testthat`): running one scenario by hand to debug optimal allocation with a single stratum, and a calibration check that the package's variance estimators are correct. Not run as part of any pipeline or CI.
- **`plotting_functions.R`** — shared plotting helpers (canonical method-level ordering/colors, frontier/coverage/bias/delta plot builders) used by the `03*` figure scripts.

## Applied data pipeline (three real case studies: READS, PERSUADE, Barbera)

Each case study follows the same three-stage shape: prep raw data → get LLM/surrogate scores → run the applied analysis.

- **`00_prep_barbera_data.R`** — cleans the raw Barbera (2020) NYT-article coder data (strips HTML breaks, normalizes punctuation) into the finite population used by the Barbera case study. (READS and PERSUADE do not have separate prep scripts; their prep happens inline at the top of the corresponding `04*`/example script.)
- **`04a_get_llm_scores_READS.R`**, **`04b_get_llm_scores_persuade.R`**, **`04c_get_llm_scores_barbera.R`** — call an LLM (via `ellmer`) with a scoring rubric to generate the surrogate outcome (`Yhat`) for each dataset, and save the scored data to `generated data/`.
- **`05a_READS_application.R`**, **`05b_persuade_application.R`**, **`05c_barbera_application.R`** — the applied analyses proper: run the precoding stratification/sampling/estimation workflow (`create_stratification_candidates()` → `compare_stratifications_precoding()` → `generate_sample_index()` → `est_ma_strat()`, etc.) on each real coded dataset and produce the paper's applied-analysis tables/figures.
- **`application_functions.R`** — shared helpers for the `05*` scripts: oracle variance calculation for SRS/stratified-proportional/stratified-optimal designs, single- and multi-sample estimate computation, and forest-plot/comparison-table builders.
- **`luke_persuade_example.R`** — standalone example/scratch script for loading and feature-generating the PERSUADE corpus (via `rcttext::generate_features()`); not part of the numbered pipeline.

## Conventions for new scripts

- Match the numeric-prefix (`NN_description.R`) or `NNx_description.R` naming convention for a new pipeline stage, and the `04*_get_llm_scores_*`/`05*_*_application` convention for a new applied dataset.
- Put shared helpers in the matching `*_functions.R` file (`sim_functions.R`, `application_functions.R`, `plotting_functions.R`) rather than inlining them in a numbered script.
- Never hardcode API keys (OpenAI, etc.) in these scripts — read them from the environment (e.g. `Sys.getenv("OPENAI_API_KEY")`, set via `.Renviron`), not `Sys.setenv()` with a literal key string.
