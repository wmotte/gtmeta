# igmi_prediction.R -- clinical prediction-model performance as IGMI
# study objects.
#
# igmi_notes Section 8.2 (Definition 8.4, Corollary 8.5). An external
# validation study reports discrimination (concordance index C) and
# calibration (observed:expected ratio O:E, or calibration slope). On the
# variance-stabilizing scales the study object is
# N(theta_i, Sigma_i(rho)) with theta_i = (logit C_i, log(O:E)_i) and
# Sigma_i(rho) = D_i R(rho) D_i built from the equicorrelation R(rho),
# exactly as igmi_study_list() builds Sigma_i(rho) for the generic
# multivariate corpus. The within-study correlation rho between
# discrimination and calibration is not identified from marginal reports
# and enters as a plug-in; by Corollary 8.5 the trace-precision-weighted
# barycenter mean (hence the IGMI pooled operating point) is exactly
# invariant to it.

#' Prediction-model-performance study objects from C-statistic and O:E
#'
#' Turns per-study discrimination and calibration summaries into IGMI study
#' objects on the `(logit C, log O:E)` scale (igmi_notes Definition 8.4).
#' Inputs are given on their natural scales and transformed by the delta
#' method: `se(logit C) = se(C) / (C(1-C))` and
#' `se(log O:E) = se(O:E) / (O:E)`. Pass already-transformed inputs with
#' `transform = FALSE` (then `c_stat`/`oe` are taken as `logit C`/`log O:E`
#' and `c_se`/`oe_se` as their SEs on those scales).
#'
#' The within-study correlation `rho` between discrimination and
#' calibration is a plug-in (unidentified from marginal reports); the
#' trace-precision-weighted BW point estimate is invariant to it
#' (Corollary 8.5), which the accompanying tests assert exactly.
#'
#' @param c_stat Concordance statistics (in `(0,1)` when `transform`).
#' @param c_se Standard errors of `c_stat` (on the same scale).
#' @param oe Observed-to-expected ratios (`> 0` when `transform`).
#' @param oe_se Standard errors of `oe` (on the same scale).
#' @param study Optional study labels (defaults to `seq_along(c_stat)`).
#' @param rho Plug-in within-study correlation between `logit C` and
#'   `log O:E` (default 0); must satisfy `-1 < rho < 1`.
#' @param transform If `TRUE` (default), delta-method transform natural-
#'   scale inputs to `(logit C, log O:E)`; if `FALSE`, inputs are already
#'   on those scales.
#' @return List of study objects `list(study =, theta =, Sigma =)`, the
#'   same contract as [igmi_study_list()].
#' @seealso [igmi_study_list()], [igmi_dta_studies()]
#' @export
igmi_pred_studies <- function(c_stat, c_se, oe, oe_se, study = NULL,
                              rho = 0, transform = TRUE) {
  c_stat <- as.numeric(c_stat); c_se <- as.numeric(c_se)
  oe <- as.numeric(oe); oe_se <- as.numeric(oe_se)
  K <- length(c_stat)
  stopifnot(length(c_se) == K, length(oe) == K, length(oe_se) == K,
            K >= 1L)
  stopifnot(is.numeric(rho), length(rho) == 1L, rho > -1, rho < 1)
  if (is.null(study)) study <- seq_len(K)

  if (transform) {
    if (any(c_stat <= 0 | c_stat >= 1, na.rm = TRUE)) {
      stop("c_stat must lie in (0, 1) when transform = TRUE")
    }
    if (any(oe <= 0, na.rm = TRUE)) {
      stop("oe must be positive when transform = TRUE")
    }
    theta1 <- stats::qlogis(c_stat)
    s1 <- c_se / (c_stat * (1 - c_stat))     # delta method for logit
    theta2 <- log(oe)
    s2 <- oe_se / oe                          # delta method for log
  } else {
    theta1 <- c_stat; s1 <- c_se
    theta2 <- oe;     s2 <- oe_se
  }

  R <- matrix(c(1, rho, rho, 1), 2, 2)
  out <- list()
  for (i in seq_len(K)) {
    if (!is.finite(theta1[i]) || !is.finite(theta2[i]) ||
        !is.finite(s1[i]) || !is.finite(s2[i]) ||
        s1[i] <= 0 || s2[i] <= 0) {
      warning("study ", study[i], " dropped (non-finite / nonpositive SE)")
      next
    }
    D <- diag(c(s1[i], s2[i]))
    out[[length(out) + 1L]] <- list(
      study = study[i],
      theta = stats::setNames(c(theta1[i], theta2[i]),
                              c("logitC", "logOE")),
      Sigma = D %*% R %*% D
    )
  }
  out
}
