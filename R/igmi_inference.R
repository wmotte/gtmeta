# igmi_inference.R -- location inference for the IGMI barycenter
#
# Implements, and is tested against:
#   Prop 5.13  exact Gaussian law of theta_bar; chi^2_m ellipsoids
#   Prop 5.16  exact studentized pivot via the Frechet scatter:
#              t_{K-1} (m = 1), Hotelling F_{m,K-m} (m >= 2) under
#              proportional total covariances -- no plug-in T
#   Prop 5.17  estimated precision weights (m = 1): exact unbiasedness,
#              variance inflation 1 + 2 sum w(1-w)/nu, naive-interval
#              deficit factor 1 + 4 sum w(1-w)/nu (min nu >= 9)
#   Rem 5.15(a) nonparametric bootstrap over studies = Phase-3 default

#' Exact location confidence region (known covariances, fixed weights)
#'
#' `theta_bar ~ N_m(theta, sum_i w_i^2 (Sigma_i + T))` exactly at every
#' K (igmi_notes Prop 5.13), giving the exact chi^2_m confidence
#' ellipsoid; for `m = 1` with precision weights and `T = 0` this is the
#' classical fixed-effect interval.
#'
#' @param studies List of `igmi_gaussian` objects.
#' @param weights Positive weights (default equal).
#' @param T_mat Between-study covariance `T` (default 0; matrix or
#'   scalar for `m = 1`).
#' @param level Confidence level.
#' @return List: `theta_bar`, `Sigma_bar_theta` (covariance of the
#'   estimate), `chi2_quantile`, `level`; for `m = 1` also `ci`.
#' @export
igmi_location_ci <- function(studies, weights = NULL, T_mat = NULL,
                             level = 0.95) {
  m <- .check_studies(studies)
  K <- length(studies)
  w <- .check_weights(weights, K)
  if (is.null(T_mat)) T_mat <- matrix(0, m, m)
  if (!is.matrix(T_mat)) T_mat <- diag(as.numeric(T_mat), m)
  theta_bar <- Reduce(`+`, Map(function(s, wi) wi * s$theta, studies, w))
  V <- Reduce(`+`, Map(function(s, wi) wi^2 * (s$Sigma + T_mat),
                       studies, w))
  q <- stats::qchisq(level, df = m)
  out <- list(theta_bar = theta_bar, Sigma_bar_theta = V,
              chi2_quantile = q, level = level)
  if (m == 1L) {
    hw <- sqrt(q * V[1, 1])
    out$ci <- c(theta_bar - hw, theta_bar + hw)
  }
  out
}

#' Frechet scatter of the study effects
#'
#' `Omega_hat = sum_i w_i (theta_i - theta_bar)(theta_i - theta_bar)^T`
#' with `tr Omega_hat = V_loc` (igmi_notes Prop 5.16).
#'
#' @inheritParams igmi_location_ci
#' @return m x m matrix.
#' @export
frechet_scatter <- function(studies, weights = NULL) {
  m <- .check_studies(studies)
  K <- length(studies)
  w <- .check_weights(weights, K)
  theta_bar <- Reduce(`+`, Map(function(s, wi) wi * s$theta, studies, w))
  Reduce(`+`, Map(function(s, wi) {
    d <- s$theta - theta_bar
    wi * tcrossprod(d)
  }, studies, w))
}

.warn_if_pivot_conditions_fail <- function(studies, weights, tol = 1e-8) {
  m <- studies[[1]]$m
  K <- length(studies)
  Sig <- lapply(studies, `[[`, "Sigma")
  if (m == 1L) {
    v <- vapply(Sig, function(S) S[1, 1], numeric(1))
    ref <- weights * v
    ok <- max(abs(ref - mean(ref))) <= tol * max(1, abs(mean(ref)))
  } else {
    S0 <- Sig[[1]]
    denom <- sum(S0^2)
    mult <- vapply(Sig, function(S) sum(S * S0) / denom, numeric(1))
    cov_ok <- all(vapply(seq_len(K), function(i) {
      max(abs(Sig[[i]] - mult[i] * S0)) <=
        tol * max(1, max(abs(Sig[[i]])))
    }, logical(1)))
    ref <- weights * mult
    weight_ok <- max(abs(ref - mean(ref))) <= tol * max(1, abs(mean(ref)))
    ok <- cov_ok && weight_ok
  }
  if (!ok) {
    warning("igmi_pivot_ci() exact Hotelling/t pivot requires proportional ",
            "total covariances and inverse-proportional weights; use ",
            "igmi_location_ci() with known covariances or igmi_bootstrap() ",
            "outside these conditions")
  }
  invisible(ok)
}

#' Exact studentized confidence region via the Frechet scatter
#'
#' Under proportional total covariances (`V_i = c_i V_0`, weights
#' proportional to `1/c_i` -- automatic for `m = 1` with total-variance
#' weights), igmi_notes Prop 5.16 gives, exactly at every `K >= m + 1`,
#' `((K-m)/m) (theta_bar - theta)' Omega_hat^{-1} (theta_bar - theta)
#' ~ F_{m, K-m}`. For `m = 1` this is
#' `theta_bar +/- t_{K-1} sqrt(V_loc / (K-1))` -- the HKSJ interval in
#' IGMI notation (Rem 5.18(a)): no estimate of the between-study
#' covariance `T` enters. Outside proportionality the region is
#' approximate; [igmi_bootstrap()] is the Phase-3 default.
#'
#' @inheritParams igmi_location_ci
#' @return List: `theta_bar`, `Omega` (Frechet scatter),
#'   `f_quantile` (`F_{m,K-m}` at `level`), `level`, `df`; for `m = 1`
#'   also `ci` and `se` (`sqrt(V_loc/(K-1))`).
#' @export
igmi_pivot_ci <- function(studies, weights = NULL, level = 0.95) {
  m <- .check_studies(studies)
  K <- length(studies)
  stopifnot(K >= m + 1)
  w <- .check_weights(weights, K)
  .warn_if_pivot_conditions_fail(studies, w)
  theta_bar <- Reduce(`+`, Map(function(s, wi) wi * s$theta, studies, w))
  Omega <- frechet_scatter(studies, w)
  fq <- stats::qf(level, df1 = m, df2 = K - m)
  out <- list(theta_bar = theta_bar, Omega = Omega, f_quantile = fq,
              level = level, df = c(m, K - m))
  if (m == 1L) {
    se <- sqrt(Omega[1, 1] / (K - 1))
    tq <- stats::qt(1 - (1 - level) / 2, df = K - 1)
    out$se <- se
    out$ci <- c(theta_bar - tq * se, theta_bar + tq * se)
  }
  out
}

#' Coverage correction factors for estimated precision weights (m = 1)
#'
#' igmi_notes Prop 5.17(3): with `hat s_i^2 ~ s_i^2 chi^2_{nu_i}/nu_i`
#' and plug-in precision weights, the true variance of the pooled
#' estimate is inflated by `1 + 2 sum_i w_i(1-w_i)/nu_i + O(nu^-2)`,
#' the naive plug-in variance `1/hat S` is biased down by the same
#' leading amount, and the naive interval's variance deficit is the
#' factor `1 + 4 sum_i w_i(1-w_i)/nu_i + O(nu^-2)`. The point estimate
#' stays exactly unbiased (Prop 5.17(1)). Requires `min(nu) >= 9` for
#' the remainder control.
#'
#' @param weights True (or plug-in) normalised precision weights.
#' @param nu Within-study degrees of freedom (recycled to length of
#'   `weights`).
#' @return List: `variance_inflation`, `naive_bias_factor`,
#'   `deficit_factor` (multiply the naive squared standard error by this
#'   to correct the interval at leading order).
#' @export
igmi_weight_correction <- function(weights, nu) {
  w <- weights / sum(weights)
  nu <- rep_len(nu, length(w))
  if (min(nu) < 9) {
    warning("Prop 5.17(3) requires min(nu) >= 9 for the O(nu^-2) ",
            "remainder control")
  }
  s <- sum(w * (1 - w) / nu)
  list(variance_inflation = 1 + 2 * s,
       naive_bias_factor = 1 - 2 * s,
       deficit_factor = 1 + 4 * s)
}

#' Nonparametric bootstrap over studies for the IGMI barycenter
#'
#' The Phase-3 default of igmi_notes Rem 5.15(a): resample the pairs
#' `(theta_i, Sigma_i)` with replacement, recompute the
#' Bures-Wasserstein barycenter by the fixed-point iteration, and report
#' percentile intervals. Weights are recomputed on each resample when
#' `weight_fun` is supplied (default: precision weights), matching the
#' plug-in practice whose leading-order cost is quantified by
#' [igmi_weight_correction()].
#'
#' @inheritParams igmi_location_ci
#' @param B Number of bootstrap resamples.
#' @param weight_fun Function `studies -> weights` applied to each
#'   resample; default [igmi_precision_weights()]. Pass `NULL` to keep
#'   the (resampled) original weights.
#' @param level Confidence level for the percentile intervals.
#' @return List: `theta_draws` (B x m), `scale_draws` (B, `tr Sigma_bar`),
#'   `theta_ci` (percentile, per component), `scale_ci`, `B`.
#' @export
igmi_bootstrap <- function(studies, weights = NULL, B = 2000L,
                           weight_fun = igmi_precision_weights,
                           level = 0.95) {
  m <- .check_studies(studies)
  K <- length(studies)
  w <- .check_weights(weights, K)
  alpha <- 1 - level
  theta_draws <- matrix(NA_real_, B, m)
  scale_draws <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(K, K, replace = TRUE)
    st <- studies[idx]
    wb <- if (is.null(weight_fun)) w[idx] / sum(w[idx]) else weight_fun(st)
    bar <- bw_barycenter(st, wb)
    theta_draws[b, ] <- bar$theta
    scale_draws[b] <- sum(diag(bar$Sigma))
  }
  theta_ci <- t(apply(theta_draws, 2, stats::quantile,
                      probs = c(alpha / 2, 1 - alpha / 2)))
  scale_ci <- stats::quantile(scale_draws,
                              probs = c(alpha / 2, 1 - alpha / 2))
  list(theta_draws = theta_draws, scale_draws = scale_draws,
       theta_ci = theta_ci, scale_ci = scale_ci, B = B)
}
