# simulate.R -- ground-truth generators for the Phase-4 simulation studies
#
# Four generators, one per design axis of TODO.md Phase 2:
#   sim_meta()          univariate random-effects truth (heterogeneity)
#   sim_multivariate()  correlated multi-outcome truth (IGMI setting)
#   sim_contaminated()  outlier/contamination mixtures on top of sim_meta()
#   sim_network()       pairwise contrasts on an evidence graph with
#                       injected loop inconsistency (GSMA setting; the
#                       generator is data-free, only the *empirical* GSMA
#                       arm waits on looped network data)
# All generators draw from the caller's RNG stream: set.seed() outside.

.rmvnorm <- function(n, mu, S) {
  m <- length(mu)
  # Symmetric (eigen) square root, so a singular positive-semidefinite S
  # -- e.g. a rank-deficient Tau with perfectly correlated effects -- is
  # handled per the documented contract, unlike chol() which needs S > 0.
  eS <- eigen(S, symmetric = TRUE)
  ev <- eS$values
  if (any(ev < -1e-8 * max(1, abs(ev[1])))) {
    stop("S is not positive semidefinite")
  }
  R <- eS$vectors %*% (sqrt(pmax(ev, 0)) * t(eS$vectors))
  matrix(stats::rnorm(n * m), n, m) %*% R +
    matrix(mu, n, m, byrow = TRUE)
}

#' Simulate a univariate random-effects meta-analysis
#'
#' Truth: `theta_i ~ N(mu, tau2)`, `yi | theta_i ~ N(theta_i, sei^2)`.
#' `tau2 = 0` gives the fixed-effect null used by the exact calibration
#' results (igmi_notes Prop 5.9: `V_F` null law; Prop 5.16: exact pivot).
#'
#' @param k Number of studies.
#' @param mu True common/mean effect.
#' @param tau2 Between-study variance.
#' @param sei Optional vector of within-study standard errors (length k);
#'   if `NULL`, drawn uniformly from `sei_range`.
#' @param sei_range Range to draw `sei` from when not supplied.
#' @return Data frame `(study, yi, sei, theta)` with attributes
#'   `mu`, `tau2`.
#' @export
sim_meta <- function(k, mu = 0, tau2 = 0, sei = NULL,
                     sei_range = c(0.1, 0.5)) {
  stopifnot(k >= 2L, tau2 >= 0)
  if (is.null(sei)) {
    sei <- stats::runif(k, sei_range[1], sei_range[2])
  }
  stopifnot(length(sei) == k, all(sei > 0))
  theta <- mu + stats::rnorm(k, 0, sqrt(tau2))
  yi <- theta + stats::rnorm(k, 0, sei)
  out <- data.frame(study = seq_len(k), yi = yi, sei = sei, theta = theta)
  attr(out, "mu") <- mu
  attr(out, "tau2") <- tau2
  out
}

#' Simulate a multivariate (multi-outcome) meta-analysis
#'
#' Truth: `theta_i ~ N_m(mu, Tau)`,
#' `y_i | theta_i ~ N_m(theta_i, Sigma_i)` with
#' `Sigma_i = D_i R(rho_within) D_i`, `D_i = diag(sei_i)` drawn per
#' outcome from `sei_range` and `R` an equicorrelation matrix. This is
#' the generating regime of the IGMI study objects `N(theta_i, Sigma_i)`
#' (igmi_notes section 1); the returned `Sigma` list feeds the
#' barycenter code directly.
#'
#' @param k Number of studies.
#' @param mu True mean vector (length m).
#' @param Tau Between-study covariance (m x m, positive semidefinite).
#' @param rho_within Common within-study cross-outcome correlation.
#' @param sei_range Range for per-outcome within-study standard errors.
#' @return List: `y` (k x m matrix of estimates), `Sigma` (list of k
#'   within-study covariances), `theta` (k x m true effects), `mu`, `Tau`.
#' @export
sim_multivariate <- function(k, mu, Tau, rho_within = 0,
                             sei_range = c(0.1, 0.5)) {
  m <- length(mu)
  stopifnot(k >= 2L, m >= 1L, all(dim(Tau) == c(m, m)),
            rho_within > -1 / max(1, m - 1), rho_within < 1)
  R <- matrix(rho_within, m, m)
  diag(R) <- 1
  theta <- if (all(Tau == 0)) {
    matrix(mu, k, m, byrow = TRUE)
  } else {
    .rmvnorm(k, mu, Tau)
  }
  Sigma <- vector("list", k)
  y <- matrix(NA_real_, k, m)
  for (i in seq_len(k)) {
    sei <- stats::runif(m, sei_range[1], sei_range[2])
    Sigma[[i]] <- diag(sei, m) %*% R %*% diag(sei, m)
    y[i, ] <- .rmvnorm(1, theta[i, ], Sigma[[i]])
  }
  list(y = y, Sigma = Sigma, theta = theta, mu = mu, Tau = Tau)
}

#' Simulate a contaminated meta-analysis
#'
#' Two contamination mechanisms on top of [sim_meta()], applied to a
#' random fraction `eps` of studies:
#' * a mean shift of `shift` (absolute, on the effect scale), and/or
#' * an understated standard error: the *reported* `sei` is the true
#'   sampling standard deviation divided by `understate` (so
#'   `understate = 2` means the study is twice as noisy as it claims).
#' This is the breakdown/robustness regime of TODO.md Phase 4.B.
#'
#' @inheritParams sim_meta
#' @param eps Contamination fraction in `[0, 1)`.
#' @param shift Mean shift added to contaminated studies' true effects.
#' @param understate Factor by which contaminated studies understate
#'   their standard error (`1` = none).
#' @return As [sim_meta()], plus logical column `is_outlier`; the
#'   reported `sei` column carries the understated values.
#' @export
sim_contaminated <- function(k, mu = 0, tau2 = 0, sei = NULL,
                             sei_range = c(0.1, 0.5),
                             eps = 0.1, shift = 1, understate = 1) {
  stopifnot(eps >= 0, eps < 1, understate >= 1)
  out <- sim_meta(k, mu = mu, tau2 = tau2, sei = sei,
                  sei_range = sei_range)
  n_out <- if (eps > 0) ceiling(eps * k) else 0L
  idx <- if (n_out > 0) sample.int(k, n_out) else integer(0)
  out$is_outlier <- seq_len(k) %in% idx
  if (n_out > 0) {
    true_sd <- out$sei[idx]
    out$yi[idx] <- out$theta[idx] + shift +
      stats::rnorm(n_out, 0, true_sd)
    out$sei[idx] <- true_sd / understate
  }
  attr(out, "mu") <- mu
  attr(out, "tau2") <- tau2
  out
}

#' Simulate pairwise contrasts on an evidence graph with loop inconsistency
#'
#' Truth: vertex potentials `d` (absolute effects, gauge-fixed at the
#' reference) generate consistent edge contrasts
#' `d[t2] - d[t1]`; an inconsistency vector `omega` (one entry per edge,
#' default all zero) is added on top. With `omega = 0` every cycle sum of
#' true contrasts vanishes (`H^1`-consistency, gsma_notes Thm 3.1); a
#' nonzero `omega` on a cycle edge injects exactly the loop inconsistency
#' that the GSMA harmonic statistic must detect (gsma_notes Thm 4.3).
#'
#' @param edges Two-column character matrix (or data frame) of treatment
#'   pairs `(t1, t2)`, orientation t1 -> t2.
#' @param d Named numeric vector of true absolute effects; names must
#'   cover all treatments in `edges`.
#' @param omega Numeric vector of per-edge inconsistency offsets
#'   (recycled; default 0).
#' @param n_per_edge Number of studies per edge (recycled).
#' @param sei_range Range for study standard errors.
#' @return Data frame `(study, t1, t2, edge, yi, sei)` with attributes
#'   `d`, `omega`.
#' @export
sim_network <- function(edges, d, omega = 0, n_per_edge = 3,
                        sei_range = c(0.1, 0.5)) {
  edges <- as.matrix(edges)
  stopifnot(ncol(edges) == 2L, all(edges %in% names(d)))
  n_e <- nrow(edges)
  omega <- rep_len(omega, n_e)
  n_per_edge <- rep_len(n_per_edge, n_e)
  rows <- list()
  for (e in seq_len(n_e)) {
    t1 <- edges[e, 1]; t2 <- edges[e, 2]
    truth_e <- unname(d[t2] - d[t1]) + omega[e]
    sei <- stats::runif(n_per_edge[e], sei_range[1], sei_range[2])
    rows[[e]] <- data.frame(
      t1 = t1, t2 = t2, edge = paste(t1, t2, sep = ":"),
      yi = truth_e + stats::rnorm(n_per_edge[e], 0, sei),
      sei = sei
    )
  }
  out <- do.call(rbind, rows)
  out <- cbind(study = seq_len(nrow(out)), out)
  rownames(out) <- NULL
  attr(out, "d") <- d
  attr(out, "omega") <- stats::setNames(omega,
                                        paste(edges[, 1], edges[, 2],
                                              sep = ":"))
  out
}
