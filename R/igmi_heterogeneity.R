# igmi_heterogeneity.R -- Frechet variance, decomposition, calibration
#
# Implements, and is tested against:
#   Thm 5.2   V_F = V_loc + V_scale (location/scale decomposition, BW)
#   Cor 5.3   scalar case: V_loc = Q/S, V_scale = sum w_i (s_i - s_bar)^2
#   Prop 5.9  null law S V_loc = Q ~ chi^2_{K-1} (m = 1, precision w)
#   Cor 5.10  I_F^2 = Higgins-Thompson I^2 exactly; V_F-moment estimator
#             = DerSimonian-Laird tau^2 exactly
#   Prop 5.11 multivariate null law: V_loc =d= weighted chi^2 mixture with
#             eigenvalues of M(0); E_0[V_loc] = sum w_i(1-w_i) tr Sigma_i;
#             tr-T moment estimator specializing to DL
#   Rem 5.12  (Phase-5.B update, 4 July 2026) matrix moment estimator:
#             E[Omega_hat] = sum w_i(1-w_i) Sigma_i + (1 - sum w_i^2) T,
#             T_hat = PSD-projection of the unbiased raw estimator;
#             tr T_hat_raw == tr_T_hat_raw of igmi_heterogeneity() exactly

#' Frechet variance and its location/scale decomposition
#'
#' `V_F = sum_i w_i W_2^2(mu_bar, N_i)` evaluated at the
#' Bures-Wasserstein barycenter, split as
#' `V_loc = sum_i w_i ||theta_i - theta_bar||^2` (disagreement in
#' effects) plus `V_scale = sum_i w_i B^2(Sigma_i, Sigma_bar)`
#' (disagreement in uncertainty structure) -- igmi_notes Thm 5.2. In the
#' standing regime (known covariances, fixed weights) all sampling
#' randomness sits in `V_loc`; `V_scale` is a deterministic design
#' functional invisible to every classical index (Cor 5.3).
#'
#' @param studies List of `igmi_gaussian` objects.
#' @param weights Positive weights (default equal).
#' @param barycenter Optional precomputed [bw_barycenter()] result.
#' @return List: `V_F`, `V_loc`, `V_scale`, `barycenter`, `weights`.
#' @export
frechet_variance <- function(studies, weights = NULL, barycenter = NULL) {
  K <- length(studies)
  .check_studies(studies)
  w <- .check_weights(weights, K)
  if (is.null(barycenter)) barycenter <- bw_barycenter(studies, w)
  v_loc <- sum(vapply(seq_len(K), function(i) {
    w[i] * sum((studies[[i]]$theta - barycenter$theta)^2)
  }, numeric(1)))
  v_scale <- sum(vapply(seq_len(K), function(i) {
    w[i] * bures_dist2(studies[[i]]$Sigma, barycenter$Sigma)
  }, numeric(1)))
  list(V_F = v_loc + v_scale, V_loc = v_loc, V_scale = v_scale,
       barycenter = barycenter, weights = w)
}

#' Null mean of the location Frechet variance
#'
#' `E_0[V_loc] = sum_i w_i (1 - w_i) tr Sigma_i` (igmi_notes
#' Prop 5.11(2)), the exact homogeneity-null expectation under fixed
#' weights and known covariances.
#'
#' @inheritParams frechet_variance
#' @return Scalar.
#' @export
vf_null_mean <- function(studies, weights = NULL) {
  K <- length(studies)
  .check_studies(studies)
  w <- .check_weights(weights, K)
  sum(vapply(seq_len(K), function(i) {
    w[i] * (1 - w[i]) * sum(diag(studies[[i]]$Sigma))
  }, numeric(1)))
}

#' Calibrated heterogeneity: I_F^2 and the trace moment estimator
#'
#' `I_F^2 = max(0, 1 - E_0[V_loc] / V_loc)` -- for `m = 1` with
#' precision weights this equals Higgins-Thompson `I^2 = (Q - (K-1))/Q`
#' exactly, and the moment estimator
#' `tr_T_hat = (V_loc - E_0[V_loc]) / (1 - sum w_i^2)` equals
#' DerSimonian-Laird `tau^2_DL` exactly (igmi_notes Cor 5.10,
#' Prop 5.11(2)); both identities are pinned in the tests.
#'
#' @inheritParams frechet_variance
#' @return List: `V_loc`, `E0_V_loc`, `I_F2`, `tr_T_hat` (truncated at
#'   zero), `tr_T_hat_raw`, plus the full [frechet_variance()] split.
#' @export
igmi_heterogeneity <- function(studies, weights = NULL, barycenter = NULL) {
  K <- length(studies)
  .check_studies(studies)
  w <- .check_weights(weights, K)
  fv <- frechet_variance(studies, w, barycenter)
  e0 <- vf_null_mean(studies, w)
  denom <- 1 - sum(w^2)
  if (K < 2L || denom <= .Machine$double.eps) {
    warning("igmi_heterogeneity: between-study heterogeneity is undefined ",
            "for K < 2 or weight fully concentrated on one study ",
            "(1 - sum(w^2) = 0); returning NA for I_F2 and the moment ",
            "estimator.", call. = FALSE)
    return(c(fv, list(E0_V_loc = e0, I_F2 = NA_real_,
                      tr_T_hat = NA_real_, tr_T_hat_raw = NA_real_)))
  }
  i2 <- max(0, 1 - e0 / fv$V_loc)
  raw <- (fv$V_loc - e0) / denom
  c(fv, list(E0_V_loc = e0, I_F2 = i2,
             tr_T_hat = max(0, raw), tr_T_hat_raw = raw))
}

#' Matrix method-of-moments estimator of the between-study covariance
#'
#' The matrix analogue of the trace moment estimator of Prop 5.11(2)
#' (igmi_notes Rem 5.12, Phase-5.B update): with the Frechet scatter
#' `Omega_hat = sum_i w_i (theta_i - theta_bar)(theta_i - theta_bar)'`,
#' exactly `E[Omega_hat] = sum_i w_i(1-w_i) Sigma_i + (1 - sum w_i^2) T`,
#' so the raw estimator
#' `T_raw = (Omega_hat - sum_i w_i(1-w_i) Sigma_i) / (1 - sum w_i^2)`
#' is unbiased for `T` (fixed weights, known `Sigma_i`); `T_hat` is its
#' Frobenius projection onto the PSD cone (eigenvalue truncation at
#' zero). `tr(T_raw)` equals `tr_T_hat_raw` of [igmi_heterogeneity()]
#' exactly, hence for `m = 1` with precision weights the untruncated
#' DerSimonian-Laird `tau^2`; the matrix form is a fixed-weights variant
#' of the multivariate DL family of Jackson, White & Thompson (2010).
#'
#' @inheritParams frechet_variance
#' @return List: `T_hat` (PSD-projected), `T_raw` (unbiased, possibly
#'   indefinite), `Omega` (Frechet scatter), `E0_Omega`
#'   (`sum_i w_i(1-w_i) Sigma_i`), `denom` (`1 - sum w_i^2`).
#' @export
igmi_mm_cov <- function(studies, weights = NULL) {
  m <- .check_studies(studies)
  K <- length(studies)
  w <- .check_weights(weights, K)
  Omega <- frechet_scatter(studies, w)
  E0 <- Reduce(`+`, Map(function(s, wi) wi * (1 - wi) * s$Sigma,
                        studies, w))
  denom <- 1 - sum(w^2)
  if (K < 2L || denom <= .Machine$double.eps) {
    warning("igmi_mm_cov: the between-study covariance is undefined for ",
            "K < 2 or weight fully concentrated on one study ",
            "(1 - sum(w^2) = 0); returning NA matrices.", call. = FALSE)
    NA_mat <- matrix(NA_real_, m, m)
    return(list(T_hat = NA_mat, T_raw = NA_mat, Omega = Omega,
                E0_Omega = E0, denom = denom))
  }
  T_raw <- (Omega - E0) / denom
  eg <- eigen((T_raw + t(T_raw)) / 2, symmetric = TRUE)
  T_hat <- eg$vectors %*% (pmax(eg$values, 0) * t(eg$vectors))
  T_hat <- (T_hat + t(T_hat)) / 2
  list(T_hat = T_hat, T_raw = T_raw, Omega = Omega, E0_Omega = E0,
       denom = denom)
}

#' Null-law matrix M(0) of the location Frechet variance
#'
#' The `mK x mK` matrix of igmi_notes Prop 5.11 (at `T = 0`): under
#' homogeneity `V_loc =d= sum_a lambda_a z_a^2` with `lambda_a` the
#' eigenvalues of `M(0)` and `z_a` iid standard normal.
#'
#' @inheritParams frechet_variance
#' @return List: `M` (the matrix), `lambda` (its eigenvalues, descending).
#' @export
vf_null_matrix <- function(studies, weights = NULL) {
  m <- .check_studies(studies)
  K <- length(studies)
  w <- .check_weights(weights, K)
  Sig <- lapply(studies, `[[`, "Sigma")
  Swk2 <- Reduce(`+`, Map(function(S, wi) wi^2 * S, Sig, w))
  M <- matrix(0, m * K, m * K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      blk <- -w[i] * Sig[[i]] - w[j] * Sig[[j]] + Swk2
      if (i == j) blk <- blk + Sig[[i]]
      M[(i - 1) * m + 1:m, (j - 1) * m + 1:m] <- sqrt(w[i] * w[j]) * blk
    }
  }
  M <- (M + t(M)) / 2
  lambda <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  list(M = M, lambda = pmax(0, lambda))
}

#' Null p-value for the location Frechet variance
#'
#' Upper-tail probability of the weighted chi-square mixture of
#' [vf_null_matrix()] at the observed `V_loc`, by Monte Carlo (a
#' dependency-free stand-in for an Imhof/Davies evaluation, cf.
#' igmi_notes Rem 5.12(b)). For `m = 1` with precision weights this is
#' the classical `Q`-test p-value up to MC error (Prop 5.9).
#'
#' @inheritParams frechet_variance
#' @param V_loc Observed location Frechet variance (default: computed).
#' @param nsim Monte Carlo sample size.
#' @return List: `p`, `V_loc`, `lambda`, `nsim`.
#' @export
vf_null_test <- function(studies, weights = NULL, V_loc = NULL,
                         nsim = 1e5) {
  K <- length(studies)
  w <- .check_weights(weights, K)
  if (is.null(V_loc)) V_loc <- frechet_variance(studies, w)$V_loc
  if (K < 2L) {
    warning("vf_null_test: the location heterogeneity test is undefined ",
            "for K < 2; returning NA.", call. = FALSE)
    return(list(p = NA_real_, V_loc = V_loc, lambda = numeric(0),
                nsim = nsim))
  }
  lambda <- vf_null_matrix(studies, w)$lambda
  lambda <- lambda[lambda > 1e-12 * max(lambda)]
  draws <- as.numeric(
    matrix(stats::rnorm(nsim * length(lambda))^2, nsim) %*% lambda
  )
  # (r + 1) / (nsim + 1): the Monte Carlo p-value never collapses to an
  # exact 0 (which would send -log(p) to Inf downstream).
  list(p = (sum(draws >= V_loc) + 1) / (nsim + 1), V_loc = V_loc,
       lambda = lambda, nsim = nsim)
}
