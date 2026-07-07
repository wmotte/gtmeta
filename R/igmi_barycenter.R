# igmi_barycenter.R -- Bures-Wasserstein fixed point and Fisher-Rao means
#
# Implements, and is tested against, the following results of
# math/igmi_notes.tex:
#   Thm 2.4 / Lem 3.3   BW distance splitting; mean decoupling
#   Thm 3.4 / Prop 3.5  covariance fixed point; scalar scale rule
#   Cor 3.6             scalar reduction to the fixed-effect estimate
#   Lem 7.3 / Thm 7.5 / Cor 7.8
#                       one-step descent, D(S) stopping criterion,
#                       unconditional convergence; the iteration must show
#                       F_cov nonincreasing AND tr S_k nondecreasing
#                       (AEBCM eq. (21)) -- both monitored here
#   Prop 2.2            Fisher-Rao closed forms (m = 1; centered family)
#   Thm 6.11            exact variational d_FR over Skew(m)
#   Lem 6.6 / Cor 6.21  checkable Karcher-ball criterion; admissible
#                       radius (pi/2) sqrt(7/2) for every m >= 2

# ---- symmetric matrix helpers -------------------------------------------

.eig_fun <- function(S, f) {
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  E <- e$vectors %*% (f(e$values) * t(e$vectors))
  (E + t(E)) / 2
}

.sqrtm <- function(S) .eig_fun(S, sqrt)
.invsqrtm <- function(S) .eig_fun(S, function(v) 1 / sqrt(v))
.logm_spd <- function(S) .eig_fun(S, log)
.expm_sym <- function(S) .eig_fun(S, exp)

# ---- Bures-Wasserstein ---------------------------------------------------

#' Squared Bures distance between covariance matrices
#'
#' `B^2(A, B) = tr A + tr B - 2 tr (A^{1/2} B A^{1/2})^{1/2}`
#' (igmi_notes Thm 2.4). Clamped at zero against roundoff.
#'
#' @param A,B Positive-definite matrices.
#' @return Nonnegative scalar.
#' @export
bures_dist2 <- function(A, B) {
  Ah <- .sqrtm(A)
  cross <- sum(diag(.sqrtm(Ah %*% B %*% Ah)))
  max(0, sum(diag(A)) + sum(diag(B)) - 2 * cross)
}

#' Squared 2-Wasserstein distance between Gaussian study objects
#'
#' `W_2^2 = ||theta_1 - theta_2||^2 + B^2(Sigma_1, Sigma_2)`
#' (splitting, igmi_notes Thm 2.4).
#'
#' @param N1,N2 `igmi_gaussian` objects.
#' @return Nonnegative scalar.
#' @export
bw_dist2 <- function(N1, N2) {
  sum((N1$theta - N2$theta)^2) + bures_dist2(N1$Sigma, N2$Sigma)
}

#' Bures-Wasserstein barycenter of Gaussian studies
#'
#' Mean: `theta_bar = sum_i w_i theta_i` (mean decoupling, igmi_notes
#' Lem 3.3). Covariance: the fixed point of
#' `Sigma = sum_i w_i (Sigma^{1/2} Sigma_i Sigma^{1/2})^{1/2}`
#' (Thm 3.4) computed by the iteration `S <- G(S) = Tbar S Tbar`,
#' `Tbar(S) = sum_i w_i T^{S -> Sigma_i}`, which converges
#' unconditionally from any positive-definite start (Thm 7.5 + Cor 7.8).
#'
#' Stopping criterion: `D(S) = ||(I - Tbar(S)) S^{1/2}||_F^2`, which is
#' the exact one-step decrement of the covariance functional (Lem 7.3)
#' and vanishes exactly at the barycenter (Thm 7.4). Two exact
#' monotonicity invariants are monitored every step and a warning is
#' thrown if either fails beyond roundoff: `F_cov(S_k)` nonincreasing
#' (Lem 7.3) and `tr S_k` nondecreasing (Alvarez-Esteban et al. 2016,
#' eq. (21); close-reading note after Cor 7.8).
#'
#' @param studies List of `igmi_gaussian` objects.
#' @param weights Positive weights (default equal); normalised internally.
#' @param tol Convergence tolerance on `D(S)`.
#' @param max_iter Maximum iterations.
#' @param S0 Starting covariance (default: linear mean of the
#'   `Sigma_i`).
#' @return Object of class `igmi_barycenter`: list with `theta`, `Sigma`,
#'   `weights`, `D` (final decrement), `iterations`, `converged`,
#'   `F_cov` (trace of objective values), `fixed_point_residual`.
#' @export
bw_barycenter <- function(studies, weights = NULL, tol = 1e-12,
                          max_iter = 500L, S0 = NULL) {
  m <- .check_studies(studies)
  K <- length(studies)
  w <- .check_weights(weights, K)

  theta_bar <- Reduce(`+`, Map(function(s, wi) wi * s$theta, studies, w))

  Sig <- lapply(studies, `[[`, "Sigma")
  if (is.null(S0)) S0 <- Reduce(`+`, Map(`*`, Sig, w))
  S <- .check_spd(S0, "S0")

  fcov <- function(S) {
    sum(vapply(seq_len(K), function(i) w[i] * bures_dist2(S, Sig[[i]]),
               numeric(1)))
  }

  f_trace <- fcov(S)
  D <- Inf
  converged <- FALSE
  iter <- 0L
  while (iter < max_iter) {
    iter <- iter + 1L
    Sh <- .sqrtm(S)
    Shi <- .invsqrtm(S)
    M <- Reduce(`+`, lapply(seq_len(K), function(i) {
      w[i] * .sqrtm(Sh %*% Sig[[i]] %*% Sh)
    }))
    Tbar <- Shi %*% M %*% Shi
    Tbar <- (Tbar + t(Tbar)) / 2
    D <- sum(((diag(m) - Tbar) %*% Sh)^2)
    if (D < tol) {
      converged <- TRUE
      break
    }
    S_new <- Tbar %*% S %*% Tbar
    S_new <- (S_new + t(S_new)) / 2
    f_new <- fcov(S_new)
    slack <- 1e-10 * max(1, abs(f_trace[length(f_trace)]))
    if (f_new > f_trace[length(f_trace)] + slack) {
      warning("F_cov increased along the iteration (Lem 7.3 violated ",
              "beyond roundoff) - check inputs/conditioning")
    }
    # the trace invariant tr S_k <= tr S_{k+1} (AEBCM eq. (21)) holds
    # along iterates k >= 1; the first step from an arbitrary S0 may
    # decrease the trace (m = 1: from sum w s^2 to (sum w s)^2, Jensen)
    if (iter > 1L &&
        sum(diag(S_new)) < sum(diag(S)) - 1e-10 * max(1, sum(diag(S)))) {
      warning("tr S decreased along the iteration (AEBCM eq. (21) ",
              "violated beyond roundoff)")
    }
    S <- S_new
    f_trace <- c(f_trace, f_new)
  }
  if (!converged) {
    warning("fixed-point iteration did not reach tol in ", max_iter,
            " iterations (D = ", format(D), ")")
  }

  resid <- {
    Sh <- .sqrtm(S)
    Mfp <- Reduce(`+`, lapply(seq_len(K), function(i) {
      w[i] * .sqrtm(Sh %*% Sig[[i]] %*% Sh)
    }))
    sqrt(sum((S - Mfp)^2))
  }

  structure(list(theta = theta_bar, Sigma = S, weights = w, D = D,
                 iterations = iter, converged = converged,
                 F_cov = f_trace, fixed_point_residual = resid, m = m),
            class = "igmi_barycenter")
}

# ---- Fisher-Rao distances ------------------------------------------------

#' Fisher-Rao distance for m = 1 (closed form)
#'
#' `sqrt(2) arccosh(1 + ((mu1-mu2)^2/2 + (s1-s2)^2) / (2 s1 s2))`
#' (igmi_notes Prop 2.2(1); sqrt(2) x Poincare half-plane).
#'
#' @param mu1,s1,mu2,s2 Means and standard deviations.
#' @return Nonnegative scalar.
#' @export
fr_dist_m1 <- function(mu1, s1, mu2, s2) {
  stopifnot(s1 > 0, s2 > 0)
  sqrt(2) * acosh(1 + ((mu1 - mu2)^2 / 2 + (s1 - s2)^2) / (2 * s1 * s2))
}

#' Fisher-Rao distance for the centered family (closed form)
#'
#' `(1/sqrt(2)) ||log(S1^{-1/2} S2 S1^{-1/2})||_F`
#' (igmi_notes Prop 2.2(2); affine-invariant metric scaled by 1/2).
#'
#' @param S1,S2 Positive-definite covariances.
#' @return Nonnegative scalar.
#' @export
fr_dist_centered <- function(S1, S2) {
  Shi <- .invsqrtm(S1)
  sqrt(sum(.logm_spd(Shi %*% S2 %*% Shi)^2)) / sqrt(2)
}

# lift z(C) = M_{theta,C} diag(Sigma^{-1}, 1, Sigma) M_{theta,C}^T of
# igmi_notes eq. (6.2) (normal form, Lem 6.9); C skew-symmetric m x m
.fr_lift <- function(theta, Sigma, C = NULL) {
  m <- length(theta)
  if (is.null(C)) C <- matrix(0, m, m)
  M <- diag(2 * m + 1)
  M[m + 1, 1:m] <- theta
  M[(m + 2):(2 * m + 1), 1:m] <- -tcrossprod(theta) / 2 + C
  M[(m + 2):(2 * m + 1), m + 1] <- -theta
  D <- matrix(0, 2 * m + 1, 2 * m + 1)
  D[1:m, 1:m] <- solve(Sigma)
  D[m + 1, m + 1] <- 1
  D[(m + 2):(2 * m + 1), (m + 2):(2 * m + 1)] <- Sigma
  M %*% D %*% t(M)
}

.skew_from_par <- function(par, m) {
  C <- matrix(0, m, m)
  C[upper.tri(C)] <- par
  C - t(C)
}

#' Fisher-Rao distance on the full Gaussian manifold
#'
#' Exact variational formula (igmi_notes Thm 6.11):
#' `d_FR(N1, N2) = (1/2) min over C in Skew(m) of
#' ||log(s1^{-1/2} z2(C) s1^{-1/2})||_F`, where `z2(C)` runs over the
#' explicit fiber of the Riemannian submersion (normal form of Lem 6.9)
#' and `s1 = z1(0)`. For `m = 1` the fiber is a point and the formula is
#' closed (it must agree with [fr_dist_m1()] -- tested); for `m = 2` it
#' is a one-variable minimization. The minimization over
#' `m(m-1)/2` fiber parameters is done numerically from `C = 0`
#' (the canonical section, also an upper bound).
#'
#' @param N1,N2 `igmi_gaussian` objects.
#' @return Nonnegative scalar.
#' @export
fr_dist <- function(N1, N2) {
  m <- N1$m
  stopifnot(m == N2$m)
  if (m == 1L) {
    return(fr_dist_m1(N1$theta, sqrt(N1$Sigma[1, 1]),
                      N2$theta, sqrt(N2$Sigma[1, 1])))
  }
  if (max(abs(N1$theta - N2$theta)) == 0) {
    return(fr_dist_centered(N1$Sigma, N2$Sigma))
  }
  s1 <- .fr_lift(N1$theta, N1$Sigma)
  s1i <- .invsqrtm(s1)
  obj <- function(par) {
    z2 <- .fr_lift(N2$theta, N2$Sigma, .skew_from_par(par, m))
    sqrt(sum(.logm_spd(s1i %*% z2 %*% s1i)^2))
  }
  npar <- m * (m - 1) / 2
  if (npar == 1L) {
    # m = 2: one-variable coercive minimization (Thm 6.11)
    opt <- stats::optimize(function(c) obj(c), interval = c(-20, 20))
    val <- opt$objective
    # refine around the optimum
    opt2 <- stats::optim(opt$minimum, obj, method = "BFGS")
    val <- min(val, opt2$value)
  } else {
    opt <- stats::optim(rep(0, npar), obj, method = "BFGS")
    opt2 <- stats::optim(opt$par, obj, method = "Nelder-Mead")
    val <- min(opt$value, opt2$value, obj(rep(0, npar)))
  }
  val / 2
}

# ---- Karcher-ball admissibility (Lem 6.6 + Cor 6.21) ---------------------

#' Checkable Fisher-Rao Karcher-ball criterion
#'
#' Upper-bounds each `d_FR(N*, N_i)` by the two-leg bound of igmi_notes
#' Lem 6.6 (Mahalanobis leg + affine-invariant/sqrt(2) leg) and compares
#' with the fully explicit admissible radius `pi/(2 sqrt(kappa_m)) =
#' (pi/2) sqrt(7/2) ~ 2.938` of Cor 6.21 (`kappa_m = 2/7` for every
#' `m >= 2`, Thm 6.19; no geodesic loops, Lem 6.20). Within the ball the
#' weighted Fisher-Rao Frechet mean exists, is unique, and equals the
#' Karcher mean. For `m = 1` the manifold is Hadamard and the radius is
#' infinite (Cor 6.2).
#'
#' @param studies List of `igmi_gaussian` objects.
#' @param center Optional `igmi_gaussian` ball center `N*`; default is
#'   the Bures-Wasserstein barycenter (Rem 6.7(b)).
#' @return List: `bounds` (per-study distance bounds), `radius_limit`,
#'   `admissible` (all bounds strictly below the limit), `center`.
#' @export
fr_karcher_check <- function(studies, center = NULL) {
  m <- .check_studies(studies)
  if (is.null(center)) {
    b <- bw_barycenter(studies, igmi_precision_weights(studies))
    center <- igmi_gaussian(b$theta, b$Sigma)
  }
  Sci <- .invsqrtm(center$Sigma)
  bounds <- vapply(studies, function(s) {
    maha <- sqrt(sum((Sci %*% (s$theta - center$theta))^2))
    covleg <- sqrt(sum(.logm_spd(Sci %*% s$Sigma %*% Sci)^2)) / sqrt(2)
    maha + covleg
  }, numeric(1))
  radius_limit <- if (m == 1L) Inf else (pi / 2) * sqrt(7 / 2)
  list(bounds = bounds, radius_limit = radius_limit,
       admissible = all(bounds < radius_limit), center = center)
}

# ---- Fisher-Rao means ----------------------------------------------------

#' Fisher-Rao Frechet mean
#'
#' Three regimes, matching the settled mathematics:
#' * `m = 1`: the manifold is hyperbolic (Hadamard), so the weighted
#'   Frechet mean exists and is unique (igmi_notes Cor 6.2(1)); computed
#'   by direct 2-parameter minimization of the Frechet functional with
#'   the closed-form distance [fr_dist_m1()].
#' * centered family (all `theta_i` equal): Cartan-Hadamard again
#'   (Cor 6.2(2)); computed by the classical Karcher fixed-point
#'   iteration for the affine-invariant mean,
#'   `X <- X^{1/2} exp(sum_i w_i log(X^{-1/2} Sigma_i X^{-1/2})) X^{1/2}`
#'   (the 1/sqrt(2) metric scaling does not change the minimizer).
#' * general `m >= 2`: unique within the Karcher ball of Cor 6.21;
#'   [fr_karcher_check()] is run first (warning if not verifiably
#'   admissible), then the Frechet functional is minimized numerically
#'   over `(theta, chol Sigma)` with the variational distance
#'   [fr_dist()] (nested optimization -- small `m`, small `K` only).
#'
#' @param studies List of `igmi_gaussian` objects.
#' @param weights Positive weights (default equal).
#' @param tol Convergence tolerance (gradient norm / step criterion).
#' @param max_iter Maximum iterations (centered fixed point).
#' @return Object of class `igmi_gaussian` (the mean), with attributes
#'   `frechet_value` (the minimized functional) and `regime`.
#' @export
fr_mean <- function(studies, weights = NULL, tol = 1e-10,
                    max_iter = 200L) {
  m <- .check_studies(studies)
  K <- length(studies)
  w <- .check_weights(weights, K)

  if (m == 1L) {
    mu_i <- vapply(studies, function(s) s$theta, numeric(1))
    s_i <- vapply(studies, function(s) sqrt(s$Sigma[1, 1]), numeric(1))
    obj <- function(p) {
      sum(w * fr_dist_m1(p[1], exp(p[2]), mu_i, s_i)^2)
    }
    p0 <- c(sum(w * mu_i), log(sum(w * s_i)))
    opt <- stats::optim(p0, obj, method = "BFGS",
                        control = list(reltol = 1e-14))
    out <- igmi_gaussian(opt$par[1], sei = exp(opt$par[2]))
    attr(out, "frechet_value") <- opt$value
    attr(out, "regime") <- "m1_hyperbolic"
    return(out)
  }

  thetas <- vapply(studies, `[[`, numeric(m), "theta")
  centered <- max(abs(thetas - thetas[, 1])) == 0
  if (centered) {
    X <- Reduce(`+`, Map(function(s, wi) wi * s$Sigma, studies, w))
    for (k in seq_len(max_iter)) {
      Xh <- .sqrtm(X)
      Xhi <- .invsqrtm(X)
      L <- Reduce(`+`, Map(function(s, wi) {
        wi * .logm_spd(Xhi %*% s$Sigma %*% Xhi)
      }, studies, w))
      X <- Xh %*% .expm_sym(L) %*% Xh
      X <- (X + t(X)) / 2
      if (sqrt(sum(L^2)) < tol) break
    }
    out <- igmi_gaussian(studies[[1]]$theta, X)
    attr(out, "frechet_value") <- sum(vapply(seq_len(K), function(i) {
      w[i] * fr_dist_centered(X, studies[[i]]$Sigma)^2
    }, numeric(1)))
    attr(out, "regime") <- "centered_hadamard"
    return(out)
  }

  chk <- fr_karcher_check(studies)
  if (!chk$admissible) {
    warning("configuration not verifiably inside the Karcher ball of ",
            "Cor 6.21 (max bound ", format(max(chk$bounds)),
            " >= ", format(chk$radius_limit),
            "); the computed mean may not be the unique Frechet mean")
  }
  start <- chk$center
  pack <- function(theta, Sigma) c(theta, chol(Sigma)[upper.tri(Sigma, diag = TRUE)])
  unpack <- function(p) {
    theta <- p[1:m]
    R <- matrix(0, m, m)
    R[upper.tri(R, diag = TRUE)] <- p[-(1:m)]
    list(theta = theta, Sigma = crossprod(R))
  }
  obj <- function(p) {
    q <- unpack(p)
    if (min(diag(q$Sigma)) <= 0 || min(eigen(q$Sigma, symmetric = TRUE,
        only.values = TRUE)$values) <= 1e-12) {
      return(1e10)
    }
    Nq <- igmi_gaussian(q$theta, q$Sigma)
    sum(vapply(seq_len(K), function(i) {
      w[i] * fr_dist(Nq, studies[[i]])^2
    }, numeric(1)))
  }
  opt <- stats::optim(pack(start$theta, start$Sigma), obj,
                      method = "Nelder-Mead",
                      control = list(maxit = 2000, reltol = 1e-10))
  q <- unpack(opt$par)
  out <- igmi_gaussian(q$theta, q$Sigma)
  attr(out, "frechet_value") <- opt$value
  attr(out, "regime") <- "karcher_ball"
  out
}
