# Theorem anchors for wfr_pool(), the one-atom WFR Frechet mean:
#   Prop 4.8      the pooling objective is the sum of two-atom closed
#                 forms; the profiled mass sqrt(b) = g(x) solves the
#                 mass minimization; cross-checked against the exact
#                 POT solver where the pinned env exists
#   Prop 4.8(2)   a study beyond the transport cutoff pi*delta has
#                 exactly zero influence (redescending/breakdown)
#   Thm 4.12 /    delta -> Inf recovers W_2, whose one-atom Frechet
#   Cor 4.13      mean under precision weights is the fixed-effect
#                 estimate = scalar BW barycenter mean (Cor 3.6)
# Only the solver cross-check needs Python; everything else is exact R.

skip_if_no_wfr <- function() {
  skip_if_not_installed("reticulate")
  skip_if(!wfr_available(), "pinned WFR virtualenv not set up")
}

# the raw two-parameter objective J(x, b) via the Prop 4.8 oracle
pool_objective <- function(x, b, yi, lam, delta) {
  sum(lam * vapply(yi, function(y) {
    wfr_two_atom_oracle(1, b, abs(y - x), delta)
  }, numeric(1)))
}

test_that("profiled mass and value match the Prop 4.8 objective", {
  set.seed(41)
  yi <- rnorm(7); sei <- runif(7, 0.1, 0.5)
  lam <- (1 / sei^2) / sum(1 / sei^2)
  delta <- 0.8
  for (x in c(-0.5, 0, 0.4)) {
    th <- pmin(abs(yi - x) / (2 * delta), pi / 2)
    g <- sum(lam * cos(th))
    # closed-form profile: b_hat = g^2 beats a numeric b-minimization
    num <- optimize(function(b) pool_objective(x, b, yi, lam, delta),
                    interval = c(0, 4))
    expect_equal(num$minimum, g^2, tolerance = 1e-5)
    expect_equal(num$objective, 4 * delta^2 * (1 - g^2),
                 tolerance = 1e-8)
  }
  # and the reported optimum is the objective at (estimate, mass_ratio)
  fit <- wfr_pool(yi, sei, delta = delta)
  expect_equal(fit$value,
               pool_objective(fit$estimate, fit$mass_ratio, yi, lam,
                              delta),
               tolerance = 1e-10)
  expect_equal(fit$mass_ratio, fit$g_max^2, tolerance = 1e-12)
})

test_that("the pooling objective matches the exact POT solver", {
  skip_if_no_wfr()
  set.seed(42)
  yi <- rnorm(5); sei <- runif(5, 0.2, 0.4)
  lam <- (1 / sei^2) / sum(1 / sei^2)
  delta <- 0.6
  fit <- wfr_pool(yi, sei, delta = delta)
  solver <- sum(lam * vapply(yi, function(y) {
    wfr_dist(1, y, fit$mass_ratio, fit$estimate, delta = delta)$wfr2
  }, numeric(1)))
  expect_equal(fit$value, solver, tolerance = 1e-6)
})

test_that("delta -> Inf recovers the FE estimate (Cor 4.13 + Cor 3.6)", {
  set.seed(43)
  yi <- rnorm(9, 0.3, 0.4); sei <- runif(9, 0.1, 0.5)
  fe <- fit_fe(yi, sei)$estimate
  # exact at delta = Inf
  fit_inf <- wfr_pool(yi, sei, delta = Inf)
  expect_equal(fit_inf$estimate, fe, tolerance = 1e-12)
  expect_equal(fit_inf$mass_ratio, 1)
  # at large finite delta the profiled objective is numerically flat
  # near its maximum (curvature O(delta^-2)), so pin the estimate to
  # the stationarity oracle g'(x) = 0 instead of to FE directly ...
  delta <- 100
  lam <- (1 / sei^2) / sum(1 / sei^2)
  star <- uniroot(function(x) sum(lam * sin((yi - x) / (2 * delta))),
                  range(yi), tol = 1e-14)$root
  expect_equal(wfr_pool(yi, sei, delta = delta)$estimate, star,
               tolerance = 1e-4)
  # ... where the oracle itself has already converged to FE
  # (leading bias is -m3 / (24 delta^2), third weighted moment)
  expect_equal(star, fe, tolerance = 1e-6)
  d_small <- abs(wfr_pool(yi, sei, delta = 5)$estimate - fe)
  d_large <- abs(wfr_pool(yi, sei, delta = 50)$estimate - fe)
  expect_lt(d_large, d_small + 1e-12)
})

test_that("a study beyond the cutoff has zero influence (Prop 4.8(2))", {
  set.seed(44)
  yi <- rnorm(8, 0, 0.3); sei <- runif(8, 0.1, 0.3)
  delta <- 1
  f20 <- wfr_pool(c(yi, 20), c(sei, 0.1), delta = delta)
  f200 <- wfr_pool(c(yi, 200), c(sei, 0.1), delta = delta)
  fcl <- wfr_pool(yi, sei, delta = delta)
  # moving the far outlier does not move the estimate at all, and the
  # estimate equals the clean-data estimate (argmax g is invariant
  # under the positive rescaling induced by weight renormalization)
  expect_equal(f20$estimate, f200$estimate, tolerance = 1e-5)
  expect_equal(f20$estimate, fcl$estimate, tolerance = 1e-5)
  # the kill shows up in the diagnostics, not the estimate (Rem 4.15)
  expect_identical(f20$n_active, 8L)
  expect_false(f20$active[9])
  expect_lt(f20$mass_ratio, fcl$mass_ratio)
  # FE, by contrast, is dragged by the same outlier
  fe_cl <- fit_fe(yi, sei)$estimate
  fe_20 <- fit_fe(c(yi, 20), c(sei, 0.1))$estimate
  expect_gt(abs(fe_20 - fe_cl), 100 * abs(f20$estimate - fcl$estimate))
})

test_that("two-point configurations match the stationarity oracle", {
  delta <- 1
  # equal precisions: the estimate is the midpoint by symmetry
  fit <- wfr_pool(c(0, 1), weights = c(1, 1), delta = delta)
  expect_equal(fit$estimate, 0.5, tolerance = 1e-6)
  # unequal precisions: g'(x) = 0 <=> lam0 sin(x/(2d)) =
  # lam1 sin((1-x)/(2d)) -- solve independently with uniroot
  lam <- c(2, 1) / 3
  root <- uniroot(function(x) {
    lam[1] * sin(x / (2 * delta)) - lam[2] * sin((1 - x) / (2 * delta))
  }, c(1e-9, 1 - 1e-9), tol = 1e-12)$root
  fit2 <- wfr_pool(c(0, 1), weights = c(2, 1), delta = delta)
  expect_equal(fit2$estimate, root, tolerance = 1e-6)
})

test_that("sandwich se: closed-form limit, invariance, MC calibration", {
  # (igmi_notes Rem 4.9, Phase-5.B update)
  set.seed(42)
  yi <- rnorm(12, 0.3, 0.5)
  sei <- runif(12, 0.1, 0.4)
  lam <- (1 / sei^2) / sum(1 / sei^2)

  # delta = Inf: se is exactly the HC0 variance of the weighted mean
  f_inf <- wfr_pool(yi, sei, delta = Inf)
  expect_equal(f_inf$se,
               sqrt(sum(lam^2 * (yi - f_inf$estimate)^2)),
               tolerance = 1e-12)
  # large delta approaches the same closed form; tolerance is set by the
  # optimizer's location precision on the numerically flat profiled
  # objective (O(1e-4) at delta = 1e4, cf. the stationarity-oracle test),
  # which enters the se at first order -- not by the formula itself
  f_big <- wfr_pool(yi, sei, delta = 1e4)
  expect_equal(f_big$se, f_inf$se, tolerance = 1e-3)

  # equal weights: a dead outlier changes neither estimate nor se
  # (the 1/K factors cancel between numerator and denominator)
  y0 <- yi[1:9]
  cl <- wfr_pool(y0, weights = rep(1, 9), delta = 0.5)
  ct <- wfr_pool(c(y0, 50), weights = rep(1, 10), delta = 0.5)
  expect_identical(ct$n_active, cl$n_active)
  expect_equal(ct$estimate, cl$estimate, tolerance = 1e-6)
  expect_equal(ct$se, cl$se, tolerance = 1e-5)

  # MC calibration on clean data (all points active, boundary mass
  # negligible): mean sandwich se matches the empirical sd of the
  # estimator to a few percent
  set.seed(4242)
  K <- 40
  reps <- 2000
  res <- vapply(seq_len(reps), function(r) {
    y <- rnorm(K)
    f <- wfr_pool(y, weights = rep(1, K), delta = 3)
    c(f$estimate, f$se)
  }, numeric(2))
  expect_equal(mean(res[2, ]), sd(res[1, ]), tolerance = 0.05)
})
