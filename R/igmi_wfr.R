# igmi_wfr.R -- Wasserstein-Fisher-Rao distances via the pinned POT bridge
#
# Implements the POT recipe of igmi_notes Rem 4.9(b), in the fixed
# convention of Def 4.1:
#   cost      l_delta(x, y) = -log cos^2( (||x-y|| / (2 delta)) ^ pi/2 )
#   marginals unnormalized KL with reg_m = 1
#   value     4 delta^2 x (optimal entropy-transport value)
# Exact (non-entropic) solver: ot.unbalanced.mm_unbalanced with
# div = 'kl' and reg = 0; the objective is recomputed in R from the
# returned plan, so solver inexactness can only overestimate the value.
#
# Oracles used by test_wfr_limit.R:
#   Prop 4.8  two-atom closed form 4 delta^2 (a0 + a1
#             - 2 sqrt(a0 a1) cos((d/(2 delta)) ^ pi/2)); unit-atom
#             birth cost 2 delta; transport cutoff exactly pi delta
#   Cor 4.13  monotone bracket WFR_{delta1} <= WFR_{delta2} <= W_2 for
#             Gaussian pairs, ceiling in closed form (Thm 2.4)
#   Rem 4.9(b)(ii) pairs beyond the cutoff must carry no plan mass
#
# wfr_pool() (tested in test_wfr_pool.R) needs no Python at all: the
# one-atom Frechet objective is closed-form via Prop 4.8.
#
# The environment is the pinned virtualenv of wfr_setup() (wfr_bridge.R).

#' WFR ground cost matrix
#'
#' `l_delta(x, y) = -log cos^2(min(||x - y|| / (2 delta), pi/2))`
#' (igmi_notes Rem 4.9(a)); the `+Inf` beyond the cutoff `pi delta` is
#' encoded as the finite `cap` (Rem 4.9(b)(ii)) -- callers must check
#' that the returned plan puts no mass on capped pairs, which
#' [wfr_dist()] does.
#'
#' @param x,y Support points: numeric vectors (m = 1) or matrices with
#'   one point per row.
#' @param delta WFR length scale.
#' @param cap Finite stand-in for `+Inf` beyond the cutoff.
#' @return List: `M` (cost matrix), `capped` (logical matrix).
#' @export
wfr_cost <- function(x, y, delta, cap = 1e6) {
  if (!is.matrix(x)) x <- matrix(x, ncol = 1)
  if (!is.matrix(y)) y <- matrix(y, ncol = 1)
  stopifnot(ncol(x) == ncol(y), delta > 0)
  d2 <- outer(rowSums(x^2), rowSums(y^2), `+`) - 2 * tcrossprod(x, y)
  d <- sqrt(pmax(d2, 0))   # pmax(d2, .) keeps the matrix dims
  theta <- d / (2 * delta)
  capped <- theta >= pi / 2
  M <- -log(pmax(cos(pmin(theta, pi / 2)), 0)^2)
  M[capped] <- cap
  list(M = M, capped = capped)
}

#' WFR distance between discrete measures (exact solver via POT)
#'
#' Computes `WFR_delta(sum a_i delta_{x_i}, sum b_j delta_{y_j})` by the
#' entropy-transport form (igmi_notes Rem 4.9(a)): minimize
#' `KL(gamma_1 | a) + KL(gamma_2 | b) + <gamma, l_delta>` over
#' nonnegative plans, value times `4 delta^2`. Uses POT's exact
#' majorization-minimization solver (`mm_unbalanced`, `div = 'kl'`,
#' `reg = 0`); the objective is recomputed from the returned plan with
#' the unnormalized KL, so the result is a certified upper bound that
#' matches the optimum to solver tolerance. Errors if the plan puts
#' mass on beyond-cutoff pairs (Rem 4.9(b)(ii)).
#'
#' @param a,b Nonnegative masses.
#' @param x,y Support points (vectors or matrices, one point per row).
#' @param delta WFR length scale.
#' @param max_iter,tol Solver control.
#' @return List: `wfr2` (squared distance), `wfr`, `plan`, `delta`.
#' @export
wfr_dist <- function(a, x, b, y, delta, max_iter = 20000L, tol = 1e-15) {
  stopifnot(wfr_available())
  ot <- reticulate::import("ot")
  np <- reticulate::import("numpy")
  cost <- wfr_cost(x, y, delta)
  G <- ot$unbalanced$mm_unbalanced(
    np$array(as.numeric(a)), np$array(as.numeric(b)),
    np$array(cost$M), reg_m = 1.0, div = "kl",
    numItermax = as.integer(max_iter), stopThr = tol
  )
  G <- as.matrix(G)
  if (any(G[cost$capped] > 1e-10 * max(sum(a), sum(b)))) {
    stop("optimal plan puts mass beyond the transport cutoff pi*delta ",
         "(igmi Rem 4.9(b)(ii)) - increase the cost cap")
  }
  kl <- function(g, mu) {
    # unnormalized KL(g | mu) = sum g log(g/mu) - g + mu, with 0 log 0 = 0
    pos <- g > 0
    sum(g[pos] * log(g[pos] / mu[pos])) - sum(g) + sum(mu)
  }
  g1 <- rowSums(G)
  g2 <- colSums(G)
  val <- sum(G * cost$M) + kl(g1, as.numeric(a)) + kl(g2, as.numeric(b))
  wfr2 <- 4 * delta^2 * val
  list(wfr2 = wfr2, wfr = sqrt(max(0, wfr2)), plan = G, delta = delta)
}

#' Two-atom WFR closed form (oracle)
#'
#' `WFR_delta^2(a0 delta_{x0}, a1 delta_{x1}) = 4 delta^2 (a0 + a1 -
#' 2 sqrt(a0 a1) cos(min(d/(2 delta), pi/2)))` with `d = ||x0 - x1||`
#' (igmi_notes Prop 4.8). Special cases: co-located atoms give the
#' Hellinger cost `2 delta |sqrt(a0) - sqrt(a1)|`; beyond the cutoff
#' `d >= pi delta` the value is the pure kill-and-create cost
#' `4 delta^2 (a0 + a1)`.
#'
#' @param a0,a1 Atom masses.
#' @param d Distance between the atoms.
#' @param delta WFR length scale.
#' @return Squared WFR distance.
#' @export
wfr_two_atom_oracle <- function(a0, a1, d, delta) {
  theta <- min(d / (2 * delta), pi / 2)
  4 * delta^2 * (a0 + a1 - 2 * sqrt(a0 * a1) * cos(theta))
}

#' WFR pooling: the one-atom WFR Frechet mean of the study estimates
#'
#' Robust pooled estimate from the WFR geometry: place each study as a
#' unit-mass atom `delta_{y_i}`, and minimize the Frechet objective
#' `sum_i lambda_i WFR_delta^2(delta_{y_i}, b delta_x)` over a single
#' candidate atom (location `x`, mass `b`). By the two-atom closed form
#' (igmi_notes Prop 4.8) the objective is
#' `4 delta^2 sum_i lambda_i (1 + b - 2 sqrt(b) c_i(x))` with
#' `c_i(x) = cos(min(|y_i - x| / (2 delta), pi/2))`; the first-order
#' condition in the mass gives `sqrt(b_hat) = g(x) := sum_i lambda_i
#' c_i(x)`, profiled value `4 delta^2 (1 - g(x)^2)`, so the estimate is
#' `argmax_x g(x)` (always inside `[min yi, max yi]`: moving toward the
#' hull increases every `c_i`).
#'
#' Two exact structural facts drive its use in Phase 4.B:
#' * `argmax g` is a redescending M-estimator of Andrews sine-wave type
#'   (`psi(u)` proportional to `sin(u / (2 delta))` for `|u| <= pi
#'   delta`, zero beyond): a study farther than the transport cutoff
#'   `pi delta` from the estimate (Prop 4.8(2)) has exactly zero
#'   influence -- the geometry kills its mass instead of transporting it.
#' * As `delta -> Inf`, `WFR_delta` increases to `W_2` (igmi_notes
#'   Thm 4.12 / Cor 4.13) and the objective becomes the weighted sum of
#'   `(y_i - x)^2`, so with the default precision weights the estimate
#'   converges to the fixed-effect estimate -- which is also the scalar
#'   Bures-Wasserstein barycenter mean (Cor 3.6). `delta` therefore
#'   interpolates between FE pooling and robust mode-seeking; the
#'   breakdown behaviour is quantified in `sim/sim_contamination.R`.
#'
#' Following igmi_notes Rem 4.15 (report far-field mass, never
#' renormalize), the retained-mass fraction `mass_ratio = g(x_hat)^2 =
#' b_hat < 1` is returned as a contamination diagnostic, together with
#' the set of studies within the cutoff.
#'
#' The profiled objective can be multimodal (as for any redescending
#' M-estimator); the maximizer is located on an `n_grid` grid over
#' `[min yi, max yi]` and refined with [stats::optimize()] around the
#' best grid point, so `n_grid` must resolve the bump width `2 pi
#' delta`.
#'
#' @param yi Study effect estimates.
#' @param sei Study standard errors; used for the default precision
#'   weights (may be omitted when `weights` is supplied).
#' @param delta WFR length scale (the robustness tuning constant:
#'   rejection point `pi delta` on the `yi` scale). `Inf` gives the
#'   exact fixed-effect / W2 limit.
#' @param weights Optional Frechet weights (default `1 / sei^2`);
#'   normalized internally.
#' The `se` field is the classical M-estimation sandwich standard error
#' (Huber 1964; Godambe form, Hampel et al. 1986) of the stationarity
#' equation `sum_i lambda_i sin((y_i - x)/(2 delta)) = 0` over the
#' active set (igmi_notes Rem 4.9, Phase-5.B update):
#' `se^2 = (2 delta)^2 sum_A lambda_i^2 sin^2(u_i/(2 delta)) /
#' (sum_A lambda_i cos(u_i/(2 delta)))^2`. As `delta -> Inf` it reduces
#' to the HC0 variance of the weighted mean
#' `sqrt(sum lambda_i^2 (y_i - x_hat)^2)`, and for equal weights it is
#' invariant under adding studies beyond the rejection point. First
#' order only: the active set is treated as fixed (psi jumps at the
#' cutoff), valid when little data mass sits near the boundary.
#'
#' @param n_grid Grid size for the global search.
#' @return List: `estimate`, `se` (sandwich standard error),
#'   `mass_ratio` (retained mass `b_hat`,
#'   1 = nothing killed), `value` (Frechet objective at the optimum),
#'   `g_max`, `active` (logical: within the cutoff of the estimate),
#'   `n_active`, `delta`, `weights` (normalized).
#' @export
wfr_pool <- function(yi, sei = NULL, delta, weights = NULL,
                     n_grid = 512L) {
  stopifnot(is.numeric(yi), length(yi) >= 2L, all(is.finite(yi)),
            length(delta) == 1L, delta > 0)
  if (is.null(weights)) {
    stopifnot(!is.null(sei), length(sei) == length(yi), all(sei > 0))
    weights <- 1 / sei^2
  }
  stopifnot(length(weights) == length(yi), all(weights > 0))
  lam <- weights / sum(weights)
  k <- length(yi)
  if (!is.finite(delta)) {
    est <- sum(lam * yi)
    return(list(estimate = est,
                se = sqrt(sum(lam^2 * (yi - est)^2)),
                mass_ratio = 1,
                value = sum(lam * (yi - est)^2), g_max = 1,
                active = rep(TRUE, k), n_active = k,
                delta = Inf, weights = lam))
  }
  g_fun <- function(x) {
    th <- pmin(abs(outer(yi, x, `-`)) / (2 * delta), pi / 2)
    as.numeric(crossprod(lam, cos(th)))
  }
  lo <- min(yi); hi <- max(yi)
  if (hi - lo < .Machine$double.eps * max(1, abs(lo))) {
    est <- lo
  } else {
    grid <- seq(lo, hi, length.out = n_grid)
    j <- which.max(g_fun(grid))
    est <- stats::optimize(g_fun, maximum = TRUE,
                           lower = grid[max(1L, j - 1L)],
                           upper = grid[min(n_grid, j + 1L)],
                           tol = .Machine$double.eps^0.5)$maximum
  }
  g <- g_fun(est)
  active <- abs(yi - est) < pi * delta
  u2 <- (yi[active] - est) / (2 * delta)
  se <- 2 * delta * sqrt(sum(lam[active]^2 * sin(u2)^2)) /
    sum(lam[active] * cos(u2))
  list(estimate = est, se = se, mass_ratio = g^2,
       value = 4 * delta^2 * (1 - g^2), g_max = g,
       active = active, n_active = sum(active),
       delta = delta, weights = lam)
}

#' Grid discretization of a univariate Gaussian for WFR tests
#'
#' Equal-spaced grid over `+/- span` standard deviations with weights
#' proportional to the density, normalised to total mass `mass`. Used
#' by the Phase-3 oracle tests (igmi_notes Cor 4.13); the
#' discretization error is controlled by comparing against the same
#' discretization's exact `W_2` (linear program), not the continuous
#' closed form.
#'
#' @param mu,sd Gaussian parameters.
#' @param n Number of grid points.
#' @param span Half-width of the grid in standard deviations.
#' @param mass Total mass.
#' @return List: `a` (masses), `x` (support).
#' @export
wfr_gauss_grid <- function(mu, sd, n = 200L, span = 5, mass = 1) {
  x <- seq(mu - span * sd, mu + span * sd, length.out = n)
  a <- stats::dnorm(x, mu, sd)
  list(a = mass * a / sum(a), x = x)
}

#' Exact balanced W2 between discrete measures (via POT emd)
#'
#' Reference ceiling for the monotone bracket
#' `WFR_delta <= W_2` (igmi_notes Thm 4.12(2) / Cor 4.13): the exact
#' linear-program optimal transport value with squared Euclidean cost,
#' on the same discrete supports as [wfr_dist()]. Masses must be equal
#' up to roundoff (they are renormalised to the common mass).
#'
#' @param a,b Nonnegative masses with equal totals.
#' @param x,y Support points.
#' @return Squared W2 (unnormalised convention: 1-homogeneous in mass).
#' @export
w2_dist2_discrete <- function(a, x, b, y) {
  stopifnot(wfr_available())
  stopifnot(abs(sum(a) - sum(b)) < 1e-10 * max(sum(a), sum(b)))
  ot <- reticulate::import("ot")
  np <- reticulate::import("numpy")
  if (!is.matrix(x)) x <- matrix(x, ncol = 1)
  if (!is.matrix(y)) y <- matrix(y, ncol = 1)
  d2 <- outer(rowSums(x^2), rowSums(y^2), `+`) - 2 * tcrossprod(x, y)
  mass <- sum(a)
  val <- ot$emd2(np$array(as.numeric(a) / mass),
                 np$array(as.numeric(b) / mass),
                 np$array(pmax(d2, 0)))
  mass * val
}
