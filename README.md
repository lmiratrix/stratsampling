
# stratsampling

<!-- badges: start -->
<!-- badges: end -->

Tools for **stratified sampling designs that leverage surrogate
scores** (e.g., machine learning predictions) to efficiently estimate
treatment effects when human coding is expensive.

In text-as-data experiments, researchers often have access to a
machine-coded surrogate outcome for the full sample but can only
afford to have human coders score a subset of observations. The
`stratsampling` package helps researchers decide which observations to
code and how to combine surrogate predictions with the coded subsample
to produce efficient, unbiased estimates.

The methods are described in:

> Mozer, R., Pashley, N., & Miratrix, L. (2026). Stratified sampling
> for machine-assisted text coding in experiments. *Working paper*.


## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("lmiratrix/stratsampling")
```


## Key features

- **Stratification tools**: Generate and evaluate candidate
  stratifications based on surrogate scores and covariates, using
  either pre-coding proxy metrics or oracle (full-data) evaluation.

- **Sampling functions**: Draw stratified samples with proportional
  or optimal (Neyman) allocation across strata and treatment arms.

- **Model-assisted estimators**: Compute treatment effect estimates
  that combine surrogate predictions with bias-corrected adjustments
  from the coded subsample, under both SRS and stratified designs.

- **Power analysis**: Compare the statistical power (MDES curves)
  of stratified vs. simple random sampling designs as a function of
  the number of coded units.

- **Simulation framework**: Generate populations with configurable
  surrogate accuracy, bias patterns, and residual variance
  heterogeneity for Monte Carlo simulation studies.


## Quick example

```r
library(stratsampling)

# Load built-in data
data(toy_data)
data(toy_data_coded)

# Create a 3x3 interaction stratification
toy_data$strata <- make_interaction_strata(
  x1 = toy_data$Yhat, n1 = 3, "Yhat",
  x2 = toy_data$X3, n2 = 3, "X3",
  z_vec = toy_data$Tx
)

# Draw a stratified sample of 100 units
toy_data$S <- generate_sample_index(toy_data, 100,
                                     Z = "Tx", strata = "strata",
                                     seed = 42)

# Merge coded outcomes and estimate treatment effect
full_data <- dplyr::left_join(toy_data, toy_data_coded, by = "ID")
est_ma_strat(full_data, Yobs = "Yobs", Z = "Tx",
             S = "S", strata = "strata", Yhat = "Yhat")
```


## Vignettes

- `vignette("package_overview")` -- Getting started with a
  complete workflow.
- `vignette("simulation_framework")` -- Running simulation studies
  comparing sampling designs.
- `vignette("applied_analysis")` -- Evaluating and selecting
  stratifications in practice.


## Authors

- Reagan Mozer (maintainer)
- Nicole Pashley
- Luke Miratrix
