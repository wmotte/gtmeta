# igmi_dta.R -- diagnostic test accuracy (sensitivity/specificity) as
# IGMI study objects, plus the IGMI summary ROC (SROC) curve.
#
# igmi_notes Section 8.1 (Definition 8.1, Lemma 8.2, Remark 8.3) and
# Section 8.3 (Definition 8.7, Proposition 8.8). A 2x2 count table
# (TP, FN | FP, TN) becomes a bivariate Gaussian on
# (logit-sensitivity, logit-specificity) whose within-study covariance is
# *diagonal* -- exactly, not by assumption, because sensitivity and
# specificity are estimated from disjoint patient groups (Lemma 8.2). The
# effect/variance recomputation reuses the escalc() pattern of
# read_rm5_estimates() in data_adapter.R (here measure = "PLO", the logit
# proportion). The summary ROC curve is the conditional-mean line of the
# pooled Gaussian, read off the barycenter mean and the between-study
# covariance T_hat (igmi_mm_cov, Rem 5.12); by Proposition 8.8 this is the
# Reitsma bivariate-model SROC evaluated at (theta_bar, T_hat).

#' Diagnostic-accuracy study objects from 2x2 counts
#'
#' Turns per-study diagnostic 2x2 tables into IGMI study objects on the
#' logit scale (igmi_notes Definition 8.1). Each study becomes
#' `N(theta_i, Sigma_i)` with
#' `theta_i = (logit Se_i, logit Sp_i)` and **diagonal**
#' `Sigma_i = diag(1/TP + 1/FN, 1/TN + 1/FP)`; the off-diagonal is exactly
#' zero because sensitivity and specificity come from disjoint patient
#' groups (Lemma 8.2), so unlike the generic multivariate case there is no
#' plug-in within-study correlation here.
#'
#' Cells are corrected before transformation only in studies that contain a
#' zero cell: `cc` (default 1/2) is added to all four cells of such a study
#' (Remark 8.3). Studies with a zero group margin after correction are
#' dropped with a warning.
#'
#' @param tp,fn,fp,tn Equal-length numeric vectors of true positives, false
#'   negatives, false positives, true negatives (diseased: `tp`, `fn`;
#'   non-diseased: `tn`, `fp`).
#' @param study Optional study labels (defaults to `seq_along(tp)`).
#' @param cc Continuity correction added to all four cells of any study
#'   with a zero cell (default `0.5`); set `0` to disable.
#' @return List of study objects `list(study =, theta =, Sigma =)`, the
#'   same contract as [igmi_study_list()]; convert to `igmi_gaussian`
#'   objects with `lapply(x, function(s) igmi_gaussian(s$theta, s$Sigma))`.
#' @seealso [igmi_sroc()], [igmi_study_list()]
#' @export
igmi_dta_studies <- function(tp, fn, fp, tn, study = NULL, cc = 0.5) {
  tp <- as.numeric(tp); fn <- as.numeric(fn)
  fp <- as.numeric(fp); tn <- as.numeric(tn)
  K <- length(tp)
  stopifnot(length(fn) == K, length(fp) == K, length(tn) == K, K >= 1L)
  stopifnot(is.numeric(cc), length(cc) == 1L, cc >= 0)
  if (is.null(study)) study <- seq_len(K)
  cells <- cbind(tp, fn, fp, tn)
  if (any(cells < 0, na.rm = TRUE)) stop("counts must be nonnegative")

  # continuity correction: only studies carrying a zero cell (Rem 8.3)
  has_zero <- apply(cells, 1L, function(r) any(r == 0, na.rm = TRUE)) & cc > 0
  if (any(has_zero, na.rm = TRUE)) {
    message(sum(has_zero, na.rm = TRUE),
            " study/studies with a zero cell corrected by +", cc, ".")
    cells[has_zero, ] <- cells[has_zero, ] + cc
  }
  tp <- cells[, 1]; fn <- cells[, 2]; fp <- cells[, 3]; tn <- cells[, 4]

  # logit sensitivity and specificity via escalc (measure = "PLO"),
  # correction already applied, so add = 0 (cf. read_rm5_estimates()).
  se_es <- metafor::escalc(measure = "PLO", xi = tp, mi = fn, add = 0)
  sp_es <- metafor::escalc(measure = "PLO", xi = tn, mi = fp, add = 0)
  logit_se <- as.numeric(se_es$yi); v_se <- as.numeric(se_es$vi)
  logit_sp <- as.numeric(sp_es$yi); v_sp <- as.numeric(sp_es$vi)

  out <- list()
  for (i in seq_len(K)) {
    if (!is.finite(logit_se[i]) || !is.finite(logit_sp[i]) ||
        !is.finite(v_se[i]) || !is.finite(v_sp[i]) ||
        v_se[i] <= 0 || v_sp[i] <= 0) {
      warning("study ", study[i], " dropped (degenerate after correction)")
      next
    }
    out[[length(out) + 1L]] <- list(
      study = study[i],
      theta = stats::setNames(c(logit_se[i], logit_sp[i]),
                              c("logitSe", "logitSp")),
      Sigma = diag(c(v_se[i], v_sp[i]))
    )
  }
  out
}

#' IGMI summary ROC (SROC) curve from a pooled DTA Gaussian
#'
#' Reads the summary ROC curve off the pooled diagnostic Gaussian
#' (igmi_notes Definition 8.7): the conditional-mean line of
#' `N(theta_bar, T_hat)`,
#' `logit Se = mu_S + (t_SC / t_CC)(logit Sp - mu_C)`,
#' mapped through `expit` into ROC coordinates. By Proposition 8.8 this is
#' the Reitsma bivariate-model SROC evaluated at `(theta_bar, T_hat)`.
#'
#' The mean order must be `(logit Se, logit Sp)`, as produced by
#' [igmi_dta_studies()] and carried through [bw_barycenter()].
#'
#' @param barycenter Either a [bw_barycenter()] result (its `$theta` is
#'   used) or a length-2 numeric `c(logit Se, logit Sp)`.
#' @param T_hat 2x2 between-study covariance on the logit scale, e.g.
#'   `igmi_mm_cov(studies, weights)$T_hat`.
#' @param k Half-width of the plotting range in SD units of logit Sp
#'   (default 2).
#' @param n Number of curve points (default 200).
#' @return Data frame with columns `logit_sp, logit_se, spec, sens, fpr,
#'   tpr`; the summary operating point is attached as attribute
#'   `"summary_point"` (a named numeric `c(sens, spec)`), and the slope as
#'   attribute `"slope"`.
#' @seealso [igmi_dta_studies()], [igmi_mm_cov()]
#' @export
igmi_sroc <- function(barycenter, T_hat, k = 2, n = 200L) {
  mu <- if (is.list(barycenter) && !is.null(barycenter$theta)) {
    barycenter$theta
  } else {
    as.numeric(barycenter)
  }
  mu <- as.numeric(mu)
  stopifnot(length(mu) == 2L)
  T_hat <- as.matrix(T_hat)
  stopifnot(nrow(T_hat) == 2L, ncol(T_hat) == 2L)
  mu_S <- mu[1]; mu_C <- mu[2]
  t_SC <- T_hat[1, 2]; t_CC <- T_hat[2, 2]

  # slope of E[logit Se | logit Sp]; if specificity is (near) homogeneous
  # between studies the conditional line degenerates to a horizontal one.
  slope <- if (t_CC > .Machine$double.eps) t_SC / t_CC else 0
  if (t_CC <= .Machine$double.eps) {
    warning("between-study logit-Sp variance ~ 0; SROC drawn as a ",
            "horizontal line through the summary point")
  }

  sd_C <- sqrt(max(t_CC, 0))
  half <- if (sd_C > 0) k * sd_C else 1
  logit_sp <- seq(mu_C - half, mu_C + half, length.out = n)
  logit_se <- mu_S + slope * (logit_sp - mu_C)

  spec <- stats::plogis(logit_sp)
  sens <- stats::plogis(logit_se)
  out <- data.frame(
    logit_sp = logit_sp, logit_se = logit_se,
    spec = spec, sens = sens,
    fpr = 1 - spec, tpr = sens
  )
  attr(out, "summary_point") <- stats::setNames(
    c(stats::plogis(mu_S), stats::plogis(mu_C)), c("sens", "spec"))
  attr(out, "slope") <- slope
  out
}
