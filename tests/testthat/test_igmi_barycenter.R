# Anchors: igmi_notes Thm 3.4/7.4/7.5 (fixed point + invariants),
# Prop 2.2 (FR closed forms), Thm 6.11 (variational d_FR),
# Lem 6.6 / Cor 6.21 (Karcher ball), Cor 6.2 (Hadamard means).

rand_spd <- function(m, scale = 1) {
  A <- matrix(rnorm(m * m), m)
  crossprod(A) / m * scale + diag(m) * 0.2
}

test_that("BW fixed point: invariants and residual (Thm 7.4/7.5)", {
  set.seed(21)
  m <- 3; K <- 4
  studies <- lapply(seq_len(K), function(i) {
    igmi_gaussian(rnorm(m), rand_spd(m))
  })
  w <- c(0.4, 0.3, 0.2, 0.1)
  bar <- bw_barycenter(studies, w, tol = 1e-14)
  expect_true(bar$converged)
  expect_lt(bar$D, 1e-13)
  # fixed-point equation residual (eq. 2 of the notes)
  expect_lt(bar$fixed_point_residual, 1e-6)
  # F_cov nonincreasing (Lem 7.3), exact in exact arithmetic
  expect_true(all(diff(bar$F_cov) <= 1e-10))
  # barycenter covariance is PD with the a-priori det bound (Lem 7.7)
  beta_det <- sum(w * vapply(studies, function(s) {
    det(s$Sigma)^(1 / (2 * m))
  }, numeric(1)))
  expect_gte(det(bar$Sigma), beta_det^(2 * m) * (1 - 1e-8))
})

test_that("commuting covariances give the closed form (sum w sqrt)^2", {
  # for simultaneously diagonal Sigma_i the fixed point is
  # S = (sum_i w_i Sigma_i^{1/2})^2 -- direct from eq. (2)
  studies <- list(igmi_gaussian(c(0, 0), diag(c(1, 4))),
                  igmi_gaussian(c(1, 1), diag(c(9, 1))),
                  igmi_gaussian(c(-1, 2), diag(c(4, 25))))
  w <- c(0.5, 0.3, 0.2)
  bar <- bw_barycenter(studies, w, tol = 1e-14)
  Sref <- (0.5 * diag(c(1, 2)) + 0.3 * diag(c(3, 1)) +
             0.2 * diag(c(2, 5)))^2
  expect_equal(bar$Sigma, Sref, tolerance = 1e-8)
})

test_that("variational d_FR reproduces the m = 1 closed form (Thm 6.11)", {
  # Skew(1) = {0}: the lift formula is closed and must equal Prop 2.2(1)
  cases <- list(c(0, 1, 1, 2), c(0.5, 0.8, -1.2, 1.7), c(3, 0.3, 2.9, 0.4))
  for (cc in cases) {
    s1 <- gtmeta:::.fr_lift(cc[1], matrix(cc[2]^2, 1, 1))
    z2 <- gtmeta:::.fr_lift(cc[3], matrix(cc[4]^2, 1, 1))
    s1i <- gtmeta:::.invsqrtm(s1)
    d_lift <- sqrt(sum(gtmeta:::.logm_spd(s1i %*% z2 %*% s1i)^2)) / 2
    expect_equal(d_lift, fr_dist_m1(cc[1], cc[2], cc[3], cc[4]),
                 tolerance = 1e-9)
  }
})

test_that("d_FR: centered pairs and totally geodesic mean slice (m = 2)", {
  # centered pairs: variational formula must return the affine
  # closed form (Prop 2.2(2))
  set.seed(22)
  S1 <- rand_spd(2); S2 <- rand_spd(2)
  N1 <- igmi_gaussian(c(0.3, -0.1), S1)
  N2 <- igmi_gaussian(c(0.3, -0.1), S2)
  expect_equal(fr_dist(N1, N2), fr_dist_centered(S1, S2),
               tolerance = 1e-10)

  # equal identity covariances, means differing along one axis: the
  # univariate slice {N((t,0), diag(a,1))} is totally geodesic (fixed
  # set of an isometric involution + product splitting), so d_FR equals
  # the m = 1 closed form
  d <- 1.3
  Na <- igmi_gaussian(c(0, 0), diag(2))
  Nb <- igmi_gaussian(c(d, 0), diag(2))
  expect_equal(fr_dist(Na, Nb), fr_dist_m1(0, 1, d, 1), tolerance = 1e-6)

  # symmetry of the variational formula (nontrivial: the lift is
  # asymmetric between the two arguments)
  N3 <- igmi_gaussian(c(0.5, 0.2), rand_spd(2))
  N4 <- igmi_gaussian(c(-0.3, 0.8), rand_spd(2))
  expect_equal(fr_dist(N3, N4), fr_dist(N4, N3), tolerance = 1e-6)

  # C = 0 section is an upper bound (Thm 6.11)
  s1 <- gtmeta:::.fr_lift(N3$theta, N3$Sigma)
  z2 <- gtmeta:::.fr_lift(N4$theta, N4$Sigma)
  s1i <- gtmeta:::.invsqrtm(s1)
  ub <- sqrt(sum(gtmeta:::.logm_spd(s1i %*% z2 %*% s1i)^2)) / 2
  expect_lte(fr_dist(N3, N4), ub + 1e-10)
})

test_that("d_FR m = 2 remains finite and below the C = 0 upper bound for large separation", {
  Na <- igmi_gaussian(c(0, 0), diag(2))
  Nb <- igmi_gaussian(c(100, 0), diag(2))
  d <- fr_dist(Na, Nb)
  expect_true(is.finite(d))
  expect_lte(d, fr_dist_m1(0, 1, 100, 1) + 1e-10)
})

test_that("d_FR handles general m >= 3 configurations", {
  set.seed(231)
  N1 <- igmi_gaussian(c(0.2, -0.1, 0.4), rand_spd(3))
  N2 <- igmi_gaussian(c(-0.3, 0.5, 0.1), rand_spd(3))
  d12 <- fr_dist(N1, N2)
  expect_true(is.finite(d12))
  expect_gt(d12, 0)
  expect_equal(d12, fr_dist(N2, N1), tolerance = 1e-5)
})

test_that("d_FR is invariant under the affine action", {
  # N -> (A theta + b, A Sigma A^T) is a Fisher-Rao isometry (Lem 6.4)
  set.seed(23)
  N1 <- igmi_gaussian(c(0.4, -0.2), rand_spd(2))
  N2 <- igmi_gaussian(c(-0.1, 0.6), rand_spd(2))
  A <- matrix(c(1.2, 0.3, -0.4, 0.9), 2)
  b <- c(0.7, -1.1)
  M1 <- igmi_gaussian(A %*% N1$theta + b, A %*% N1$Sigma %*% t(A))
  M2 <- igmi_gaussian(A %*% N2$theta + b, A %*% N2$Sigma %*% t(A))
  expect_equal(fr_dist(N1, N2), fr_dist(M1, M2), tolerance = 1e-5)
})

test_that("Karcher-ball check: bounds, radius, admissibility (Cor 6.21)", {
  # hand-checkable configuration
  Nstar <- igmi_gaussian(c(0, 0), diag(2))
  N1 <- igmi_gaussian(c(1, 0), diag(2))          # Mahalanobis 1, cov leg 0
  N2 <- igmi_gaussian(c(0, 0), diag(c(4, 1)))    # cov leg log(4)/sqrt(2)
  chk <- fr_karcher_check(list(N1, N2), center = Nstar)
  expect_equal(chk$bounds[1], 1, tolerance = 1e-10)
  expect_equal(chk$bounds[2], log(4) / sqrt(2), tolerance = 1e-10)
  expect_equal(chk$radius_limit, (pi / 2) * sqrt(7 / 2))
  expect_true(chk$admissible)

  # m = 1 is Hadamard: radius infinite
  chk1 <- fr_karcher_check(igmi_studies(yi = c(0, 5), sei = c(1, 1)))
  expect_equal(chk1$radius_limit, Inf)

  # a wildly spread configuration must fail the check (m >= 2)
  Nfar <- igmi_gaussian(c(50, 0), diag(2))
  chk2 <- fr_karcher_check(list(N1, Nfar), center = Nstar)
  expect_false(chk2$admissible)
})

test_that("fr_mean m = 1: Hadamard uniqueness anchors (Cor 6.2(1))", {
  # symmetric configuration: mean at mu = 0 by symmetry
  st <- igmi_studies(yi = c(-1, 1), sei = c(1, 1))
  mn <- fr_mean(st)
  expect_equal(mn$theta, 0, tolerance = 1e-6)
  # equal-weight two-point mean is the geodesic midpoint
  d12 <- fr_dist(st[[1]], st[[2]])
  expect_equal(fr_dist(mn, st[[1]]), d12 / 2, tolerance = 1e-5)
  expect_equal(fr_dist(mn, st[[2]]), d12 / 2, tolerance = 1e-5)
  # minimality against perturbations
  F_at <- function(mu, s) {
    sum(0.5 * fr_dist_m1(mu, s, c(-1, 1), c(1, 1))^2)
  }
  f0 <- attr(mn, "frechet_value")
  expect_lt(f0, F_at(0.1, sqrt(mn$Sigma[1, 1])))
  expect_lt(f0, F_at(0, sqrt(mn$Sigma[1, 1]) * 1.1))
})

test_that("fr_mean centered: Karcher fixed point (Cor 6.2(2))", {
  set.seed(24)
  S1 <- rand_spd(2); S2 <- rand_spd(2)
  st <- list(igmi_gaussian(c(0, 0), S1), igmi_gaussian(c(0, 0), S2))
  mn <- fr_mean(st)
  expect_equal(attr(mn, "regime"), "centered_hadamard")
  # two-point equal-weight mean = affine-invariant geodesic midpoint
  S1h <- gtmeta:::.sqrtm(S1)
  S1hi <- gtmeta:::.invsqrtm(S1)
  mid <- S1h %*% gtmeta:::.sqrtm(S1hi %*% S2 %*% S1hi) %*% S1h
  expect_equal(mn$Sigma, (mid + t(mid)) / 2, tolerance = 1e-7)

  # commuting case: geometric mean exp(sum w log Sigma_i)
  stD <- list(igmi_gaussian(c(0, 0), diag(c(1, 4))),
              igmi_gaussian(c(0, 0), diag(c(9, 1))))
  mnD <- fr_mean(stD, weights = c(0.3, 0.7))
  ref <- diag(exp(0.3 * log(c(1, 4)) + 0.7 * log(c(9, 1))))
  expect_equal(mnD$Sigma, ref, tolerance = 1e-7)
})

test_that("fr_mean general m = 2: equidistance and optimality", {
  st <- list(igmi_gaussian(c(0, 0), diag(2)),
             igmi_gaussian(c(1, 0.3), diag(c(1.5, 0.8))))
  mn <- fr_mean(st)
  expect_equal(attr(mn, "regime"), "karcher_ball")
  d1 <- fr_dist(mn, st[[1]])
  d2 <- fr_dist(mn, st[[2]])
  # equal-weight two-point Frechet mean is the midpoint: equidistant,
  # and the Frechet value is 2 * (1/2) * (d/2)^2 = d^2/4
  expect_equal(d1, d2, tolerance = 5e-3)
  d12 <- fr_dist(st[[1]], st[[2]])
  expect_equal(attr(mn, "frechet_value"), d12^2 / 4, tolerance = 1e-3)
})
