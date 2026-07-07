# The named Phase-3 reduction test (TODO 3.B): the scalar IGMI
# Bures-Wasserstein barycenter with precision weights equals the
# metafor fixed-effect estimate exactly (igmi_notes Cor 3.6), with the
# scale rule s_bar = sum w_i s_i (Prop 3.5) -- which is NOT the standard
# error of the pooled mean.

test_that("scalar BW barycenter == metafor fixed effect (Cor 3.6)", {
  yi <- c(0.12, -0.30, 0.25, 0.05, 0.40)
  sei <- c(0.10, 0.25, 0.15, 0.30, 0.20)
  studies <- igmi_studies(yi = yi, sei = sei)
  w <- igmi_precision_weights(studies)
  expect_equal(w, (1 / sei^2) / sum(1 / sei^2), tolerance = 1e-12)

  bar <- bw_barycenter(studies, w)
  fit <- metafor::rma(yi = yi, sei = sei, method = "EE")
  expect_equal(bar$theta, as.numeric(fit$beta), tolerance = 1e-10)

  # scale rule (Prop 3.5): sd of the barycenter = weighted mean of sds
  expect_equal(sqrt(bar$Sigma[1, 1]), sum(w * sei), tolerance = 1e-10)
  # ... and it is NOT the standard error of the pooled mean
  expect_gt(sqrt(bar$Sigma[1, 1]), fit$se)

  # m = 1 converges essentially in one step (Rem 7.6(c)):
  # G(v) = (sum w s)^2 for every v
  bar2 <- bw_barycenter(studies, w, S0 = matrix(7.3, 1, 1))
  expect_equal(bar2$Sigma[1, 1], (sum(w * sei))^2, tolerance = 1e-10)
  expect_lte(bar2$iterations, 3L)
})

test_that("mean decoupling holds for arbitrary weights (Lem 3.3)", {
  set.seed(11)
  yi <- rnorm(4)
  sei <- runif(4, 0.1, 0.5)
  w <- c(0.4, 0.3, 0.2, 0.1)
  bar <- bw_barycenter(igmi_studies(yi = yi, sei = sei), w)
  expect_equal(bar$theta, sum(w * yi), tolerance = 1e-12)
})
