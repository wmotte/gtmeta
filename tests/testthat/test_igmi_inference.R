# Anchors: igmi_notes Prop 5.13 (exact chi^2_m ellipsoid),
# Prop 5.16 (exact t/Hotelling pivot via the Frechet scatter -- the MC
# checks promoted from the scratch verification of Rem 5.18(d)),
# Prop 5.17 (estimated-weights identities and O(1/nu) corrections).

test_that("location ellipsoid is exact in MC (Prop 5.13, m = 2)", {
  set.seed(41)
  m <- 2; K <- 4
  theta <- c(0.3, -0.2)
  Sig <- list(diag(c(1, 2)) / 4, matrix(c(2, 0.6, 0.6, 1), 2) / 4,
              diag(c(0.5, 0.5)), diag(c(1.5, 0.7)))
  w <- c(0.3, 0.3, 0.2, 0.2)
  Sh <- lapply(Sig, gtmeta:::.sqrtm)
  nrep <- 3000L
  cover <- logical(nrep)
  V <- Reduce(`+`, Map(function(S, wi) wi^2 * S, Sig, w))
  Vi <- solve(V)
  q <- qchisq(0.95, df = m)
  for (r in seq_len(nrep)) {
    th <- Map(function(R) theta + as.numeric(R %*% rnorm(m)), Sh)
    tb <- Reduce(`+`, Map(`*`, th, w))
    cover[r] <- as.numeric(t(tb - theta) %*% Vi %*% (tb - theta)) <= q
  }
  expect_equal(mean(cover), 0.95, tolerance = 0.015)

  # the function returns exactly this ellipsoid
  studies <- Map(function(t, S) igmi_gaussian(t, S),
                 list(theta, theta, theta, theta), Sig)
  ci <- igmi_location_ci(studies, w)
  expect_equal(ci$Sigma_bar_theta, V, tolerance = 1e-12)
  expect_equal(ci$chi2_quantile, q)
})

test_that("location ellipsoid includes nonzero T_mat in the exact covariance", {
  theta <- c(0.2, -0.1)
  Sig <- list(diag(c(0.2, 0.3)), matrix(c(0.4, 0.05, 0.05, 0.5), 2))
  T_mat <- matrix(c(0.10, 0.03, 0.03, 0.20), 2)
  w <- c(0.25, 0.75)
  studies <- Map(function(S) igmi_gaussian(theta, S), Sig)
  ci <- igmi_location_ci(studies, weights = w, T_mat = T_mat)
  V_ref <- Reduce(`+`, Map(function(S, wi) wi^2 * (S + T_mat), Sig, w))
  expect_equal(ci$Sigma_bar_theta, V_ref, tolerance = 1e-12)
})

test_that("t pivot is exact at K = 5 with unequal variances (Prop 5.16)", {
  set.seed(42)
  K <- 5; mu <- 0.4
  sei <- c(0.10, 0.15, 0.20, 0.25, 0.30)
  w <- (1 / sei^2) / sum(1 / sei^2)
  nrep <- 4000L
  cover <- logical(nrep)
  for (r in seq_len(nrep)) {
    yi <- mu + rnorm(K, 0, sei)
    piv <- igmi_pivot_ci(igmi_studies(yi = yi, sei = sei), w)
    cover[r] <- piv$ci[1] <= mu && mu <= piv$ci[2]
  }
  expect_equal(mean(cover), 0.95, tolerance = 0.014)
})

test_that("t pivot stays exact with tau2 > 0 under total-variance
           proportionality (Prop 5.16)", {
  set.seed(43)
  K <- 5; mu <- -0.2; tau <- 0.5
  sei <- rep(0.25, K)   # equal sei: totals proportional for any tau
  nrep <- 4000L
  cover <- logical(nrep)
  for (r in seq_len(nrep)) {
    yi <- mu + rnorm(K, 0, sqrt(sei^2 + tau^2))
    piv <- igmi_pivot_ci(igmi_studies(yi = yi, sei = sei))
    cover[r] <- piv$ci[1] <= mu && mu <= piv$ci[2]
  }
  expect_equal(mean(cover), 0.95, tolerance = 0.014)
})

test_that("Hotelling-F pivot is exact at m = 2, K = 6 (Prop 5.16)", {
  set.seed(44)
  m <- 2; K <- 6
  theta <- c(0.1, 0.5)
  V0 <- matrix(c(1, 0.4, 0.4, 2), 2) / 10
  ci <- c(1, 2, 0.5, 1.5, 3, 0.8)
  w <- (1 / ci) / sum(1 / ci)
  V0h <- gtmeta:::.sqrtm(V0)
  fq <- qf(0.95, df1 = m, df2 = K - m)
  nrep <- 3000L
  cover <- logical(nrep)
  for (r in seq_len(nrep)) {
    th <- lapply(ci, function(c) {
      theta + sqrt(c) * as.numeric(V0h %*% rnorm(m))
    })
    tb <- Reduce(`+`, Map(`*`, th, w))
    Om <- Reduce(`+`, Map(function(t, wi) wi * tcrossprod(t - tb), th, w))
    stat <- (K - m) / m *
      as.numeric(t(tb - theta) %*% solve(Om) %*% (tb - theta))
    cover[r] <- stat <= fq
  }
  expect_equal(mean(cover), 0.95, tolerance = 0.015)

  # and igmi_pivot_ci exposes exactly these ingredients
  studies <- lapply(ci, function(c) igmi_gaussian(theta, c * V0))
  piv <- igmi_pivot_ci(studies, w)
  expect_equal(piv$f_quantile, fq)
  expect_equal(piv$df, c(m, K - m))
})

test_that("Hotelling-F pivot warns outside proportional-covariance conditions", {
  studies <- list(
    igmi_gaussian(c(0, 0), diag(c(1, 2))),
    igmi_gaussian(c(1, 0), diag(c(2, 1))),
    igmi_gaussian(c(0, 1), diag(c(1.5, 3)))
  )
  expect_warning(
    igmi_pivot_ci(studies, weights = c(1, 1, 1)),
    "proportional total covariances"
  )
})

test_that("m = 1 pivot se^2 equals V_loc/(K-1) (Prop 5.16 <-> Cor 5.3)", {
  yi <- c(0.2, -0.5, 0.7, 0.1)
  sei <- c(0.2, 0.3, 0.25, 0.15)
  studies <- igmi_studies(yi = yi, sei = sei)
  w <- igmi_precision_weights(studies)
  piv <- igmi_pivot_ci(studies, w)
  fv <- frechet_variance(studies, w)
  expect_equal(piv$se^2, fv$V_loc / (length(yi) - 1), tolerance = 1e-10)
  # ... which is exactly the UWLS t-interval of benchmark.R
  u <- fit_uwls(yi, sei)
  expect_equal(piv$ci, c(u$ci.lb, u$ci.ub), tolerance = 1e-10)
})

test_that("estimated-weights corrections match MC (Prop 5.17)", {
  set.seed(45)
  K <- 4
  s <- c(0.2, 0.3, 0.4, 0.25)
  nu <- 40
  S <- sum(1 / s^2)
  w <- (1 / s^2) / S
  corr <- igmi_weight_correction(w, nu)

  nrep <- 2e5
  # vectorized: hat s^2 = s^2 chi2_nu / nu; theta_i ~ N(0, s_i^2)
  hs2 <- matrix(rchisq(nrep * K, nu) / nu, nrep) *
    matrix(s^2, nrep, K, byrow = TRUE)
  what <- (1 / hs2) / rowSums(1 / hs2)
  th <- matrix(rnorm(nrep * K), nrep) * matrix(s, nrep, K, byrow = TRUE)
  tbar <- rowSums(what * th)

  # exact unbiasedness (Prop 5.17(1))
  expect_equal(mean(tbar), 0, tolerance = 4 * sd(tbar) / sqrt(nrep))
  # variance inflation 1 + 2 sum w(1-w)/nu (Prop 5.17(3), first display)
  expect_equal(S * var(tbar), corr$variance_inflation, tolerance = 0.01)
  # naive plug-in bias 1 - 2 sum w(1-w)/nu (second display)
  expect_equal(S * mean(1 / rowSums(1 / hs2)), corr$naive_bias_factor,
               tolerance = 0.005)
  # pointwise inflation inequality (Prop 5.17(2))
  expect_true(all(rowSums(what^2 * matrix(s^2, nrep, K, byrow = TRUE))
                  >= 1 / S - 1e-12))

  # nu below the remainder-control threshold warns
  expect_warning(igmi_weight_correction(w, 5), "nu")
})

test_that("bootstrap over studies returns sane percentile intervals", {
  set.seed(46)
  yi <- rnorm(8, 0.5, 0.3)
  sei <- runif(8, 0.15, 0.35)
  studies <- igmi_studies(yi = yi, sei = sei)
  bt <- igmi_bootstrap(studies, B = 200L)
  expect_equal(dim(bt$theta_draws), c(200L, 1L))
  expect_lt(bt$theta_ci[1, 1], bt$theta_ci[1, 2])
  # interval brackets the point estimate
  bar <- bw_barycenter(studies, igmi_precision_weights(studies))
  expect_gt(bar$theta, bt$theta_ci[1, 1])
  expect_lt(bar$theta, bt$theta_ci[1, 2])
  expect_true(all(bt$scale_draws > 0))
})
