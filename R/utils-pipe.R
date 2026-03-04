#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom dplyr %>% mutate filter select group_by summarise summarize
#'   ungroup rename left_join bind_rows arrange n pull slice_head across
#'   everything n_distinct sample_n if_else row_number desc first
#' @importFrom tidyr pivot_longer pivot_wider separate crossing
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map map_dfr map_dbl map2_dfr map_chr
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_col geom_boxplot
#'   geom_hline geom_vline geom_abline geom_tile geom_bar geom_text
#'   facet_wrap facet_grid labs theme theme_minimal theme_classic
#'   scale_fill_viridis_d scale_color_manual scale_y_continuous
#'   scale_x_continuous element_text expand_limits coord_flip
#' @importFrom scales percent comma
#' @importFrom glue glue
#' @importFrom stats var sd rnorm quantile lm residuals cor pnorm qnorm
#'   runif median setNames weighted.mean
#' @importFrom utils combn
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


# Global variable declarations to avoid R CMD check NOTEs for
# non-standard evaluation (NSE) column references used by dplyr/tidyr/ggplot2
utils::globalVariables(c(
  # Column names used across the package
  ".", "Z", "S", "strata", "Yobs", "Yhat", "Y0", "Y1", "ID",
  "Tx", "Yobs_full",

  # Variables from estimation functions
  "est", "SEhat", "SE", "method", "Method", "MDES",
  "h", "SE_full", "n_code",
  "Ybar", "adj", "Ybar.adj", "frac", "N", "N.z", "N.zk",
  "var.resid", "var.resid.coded", "Ybar2.coded", "Ybar.obs.zk",
  "s2.obs.zk", "Y2bar.obs", "Ybar.obs", "pt.3", "S2.hat", "n.z",

  # Variables from strata / evaluate functions
  "resid", "n", "n_k", "N_k", "var_resid", "mean_resid",
  "sd_resid", "bias", "err", "val", "stat",
  "var_e", "var_Y", "mean_Y", "weight", "weighted_var",
  "BS_var", "WS_var", "total_var", "BS_share", "WS_share",
  "value", "component", "stratification", "rank_oracle", "rank_precoding",
  "metric", "source", "name",
  "level", "level_i", "i_min",
  "strata_eval", "N_zk", "N_z", "w_k", "avg_within_var",
  "mean_resid_z", "avg_mean_resid",
  "Yhat_bar_z", "var_Yhat",

  # Variables from compute_variance_decomposition
  "mean_resid_k", "grand_mean", "var_resid_k", "N_z_total",
  "n_z", "fpc", "bs_weight", "n_zk", "term1", "term2",
  "ws_weight", "var_resid_total",

  # Variables from compare functions
  "n_strata", "min_size", "max_size", "imbalance",
  "min_var", "max_var", "ratio_var", "range_var",
  "sd_within_var", "var_mean_resid", "mean_spread",
  "BS_minus_WS", "fill_color",
  "avg_within_var_resid", "ratio_var_resid",
  "oracle_rank", "avg_within_var_Yhat", "precoding_rank",
  "rank_diff", "precoding_mean_spread", "oracle_mean_spread",

  # Variables from describe_population
  "w", "b", "sd_e", "Var_b", "E_sig2",

  # Variables from sampling functions
  "pi", ".strata", "n.zk", "Nk", "sigma_k", "alloc",
  ".micro_strata", ".Z", "p", "prop", "BS", "WS",

  # Variables from plot functions
  "strata_label", "Z_label",

  # Variables from dat_gen functions
  "noise_sd", "bias_shift",

  # Variables used in vignettes / simulation
  "iter", "R2", "R2_label", "bias_level", "var_pattern",
  "Yhat_strat", "X3_strat", "mean_Yhat",
  "ney.hat", "varY", "n_un", "n_samp",
  "resid_sd", "RMSE", "mean_est", "empirical_SE", "mean_SEhat"
))
