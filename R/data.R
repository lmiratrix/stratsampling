#' Toy dataset for stratified sampling examples
#'
#' A simulated dataset containing surrogate predictions and baseline covariates
#' for 1000 observations from a randomized controlled trial. This dataset is
#' used in the package vignettes to demonstrate the stratification, sampling,
#' and estimation workflow.
#'
#' @format A data frame with 1000 rows and 6 columns:
#' \describe{
#'   \item{X1}{Baseline covariate 1 (continuous).}
#'   \item{X2}{Baseline covariate 2 (continuous).}
#'   \item{X3}{Baseline covariate 3 (continuous). In this toy example, X3
#'     moderates the accuracy of the surrogate predictions.}
#'   \item{Yhat}{Surrogate prediction of the outcome, e.g., from a machine
#'     learning model applied to the full sample.}
#'   \item{Tx}{Treatment assignment indicator (0 = control, 1 = treatment).}
#'   \item{ID}{Unique observation identifier.}
#' }
#'
#' @examples
#' data(toy_data)
#' head(toy_data)
"toy_data"


#' Coded outcomes for the toy dataset subsample
#'
#' A subset of \code{\link{toy_data}} containing human-coded outcomes for the
#' observations selected for coding. This dataset is used in the package
#' vignettes to demonstrate merging coded outcomes back into the full dataset
#' for model-assisted estimation.
#'
#' @format A data frame with approximately 100 rows and 2 columns:
#' \describe{
#'   \item{ID}{Unique observation identifier, matching \code{toy_data$ID}.}
#'   \item{Yobs}{Human-coded outcome for the selected observation.}
#' }
#'
#' @examples
#' data(toy_data_coded)
#' head(toy_data_coded)
"toy_data_coded"
