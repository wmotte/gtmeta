# metrics.R -- evaluation metrics for the Phase-4 simulation studies
#
# Small, dependency-free building blocks: each takes vectors over
# simulation replications (or over studies) and returns a scalar or a
# short named list. metric_table() is the per-estimator convenience
# aggregator used by sim/*.R.

#' Bias of an estimator over replications
#' @param est Estimates.
#' @param truth True value (scalar or vector like `est`).
#' @return Mean of `est - truth`.
#' @export
metric_bias <- function(est, truth) mean(est - truth)

#' Mean squared error over replications
#' @inheritParams metric_bias
#' @return A single numeric value: the mean of the squared errors
#'   `(est - truth)^2`.
#' @export
metric_mse <- function(est, truth) mean((est - truth)^2)

#' Empirical coverage of interval estimates
#' @param lo,hi Interval bounds per replication.
#' @inheritParams metric_bias
#' @return Fraction of intervals containing the truth.
#' @export
metric_coverage <- function(lo, hi, truth) mean(lo <= truth & truth <= hi)

#' Mean interval width
#' @param lo,hi Interval bounds per replication.
#' @return A single numeric value: the mean of the interval widths
#'   `hi - lo`.
#' @export
metric_ci_width <- function(lo, hi) mean(hi - lo)

#' Calibration of probability integral transform (PIT) values
#'
#' Under a well-calibrated predictive distribution the PIT values are
#' U(0,1). Returns the Kolmogorov-Smirnov distance to U(0,1) (with the
#' `ks.test` p-value), plus the mean and standard deviation of the
#' corresponding z-scores `qnorm(pit)` (targets 0 and 1) -- the two
#' summaries used for the LOO calibration checks ([loo_predict()]).
#'
#' @param pit PIT values in (0, 1).
#' @return Named list: `ks`, `ks_p`, `z_mean`, `z_sd`.
#' @export
metric_calibration <- function(pit) {
  stopifnot(all(pit > 0 & pit < 1))
  ks <- suppressWarnings(stats::ks.test(pit, "punif"))
  z <- stats::qnorm(pit)
  list(ks = unname(ks$statistic), ks_p = ks$p.value,
       z_mean = mean(z), z_sd = stats::sd(z))
}

#' Detection power (or type-I error) of a test over replications
#' @param pval P-values over replications.
#' @param alpha Nominal level.
#' @return Rejection rate. Under the null this is the empirical type-I
#'   error; under an alternative it is power.
#' @export
detection_power <- function(pval, alpha = 0.05) mean(pval <= alpha)

#' Localization accuracy of per-candidate anomaly scores
#'
#' For localization tasks (which loop is inconsistent, which study is
#' contaminated): given a score matrix over replications and the index of
#' the true anomalous candidate, returns top-1 accuracy and the mean rank
#' of the true candidate (1 = always ranked first; higher scores must
#' mean more anomalous).
#'
#' @param scores Matrix (replications x candidates) of anomaly scores,
#'   or a single numeric vector for one replication.
#' @param true_index Column index (or vector per replication) of the
#'   truly anomalous candidate.
#' @return Named list: `top1`, `mean_rank`.
#' @export
localization_accuracy <- function(scores, true_index) {
  if (is.vector(scores)) scores <- matrix(scores, nrow = 1)
  n <- nrow(scores)
  true_index <- rep_len(true_index, n)
  picked <- max.col(scores, ties.method = "first")
  rank_true <- vapply(seq_len(n), function(i) {
    sum(scores[i, ] >= scores[i, true_index[i]])
  }, numeric(1))
  list(top1 = mean(picked == true_index),
       mean_rank = mean(rank_true))
}

#' Per-estimator summary table over simulation replications
#'
#' @param est,se Point estimates and their model standard errors.
#' @param lo,hi Interval bounds.
#' @param truth True value.
#' @return One-row data frame: bias, empirical SE, model-SE/empirical-SE
#'   ratio (target 1), RMSE, coverage, mean CI width.
#' @export
metric_table <- function(est, se, lo, hi, truth) {
  emp_se <- stats::sd(est)
  data.frame(
    bias = metric_bias(est, truth),
    emp_se = emp_se,
    se_ratio = mean(se) / emp_se,
    rmse = sqrt(metric_mse(est, truth)),
    coverage = metric_coverage(lo, hi, truth),
    ci_width = metric_ci_width(lo, hi)
  )
}
