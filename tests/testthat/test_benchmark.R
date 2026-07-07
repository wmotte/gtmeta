# Anchors:
# * FE == metafor rma(method = "EE") and equals the scalar IGMI
#   Bures-Wasserstein barycenter mean under precision weights
#   (igmi_notes Cor 3.6).
# * UWLS == the no-intercept lm of t-values on precisions
#   (Stanley-Doucouliagos), with se = se_FE * sqrt(Q / (K - 1)).
# * The UWLS t_{K-1} interval is exactly the Frechet-scatter pivot
#   interval of igmi_notes Prop 5.16 (m = 1): se_UWLS^2 =
#   Q / (S (K - 1)) = V_F^loc / (K - 1), so under proportional total
#   variances its coverage is exact at every K -- checked by seeded MC.

toy <- function() {
  list(yi = c(0.12, -0.30, 0.25, 0.05, 0.40, -0.10),
       sei = c(0.10, 0.25, 0.15, 0.30, 0.20, 0.12))
}

test_that("fit_fe matches metafor EE", {
  d <- toy()
  f <- fit_fe(d$yi, d$sei)
  m <- metafor::rma(yi = d$yi, sei = d$sei, method = "EE")
  expect_equal(f$estimate, as.numeric(m$beta), tolerance = 1e-10)
  expect_equal(f$se, m$se, tolerance = 1e-10)
  expect_equal(f$logLik, as.numeric(stats::logLik(m)), tolerance = 1e-8)
})

test_that("fit_re matches metafor ML incl. hand log-likelihood", {
  d <- toy()
  f <- fit_re(d$yi, d$sei, method = "ML")
  m <- metafor::rma(yi = d$yi, sei = d$sei, method = "ML")
  expect_equal(f$estimate, as.numeric(m$beta), tolerance = 1e-10)
  expect_equal(f$se, m$se, tolerance = 1e-10)
  expect_equal(unname(f$scale["tau2"]), m$tau2, tolerance = 1e-10)
  expect_equal(f$logLik, as.numeric(stats::logLik(m)), tolerance = 1e-6)
})

test_that("fit_uwls equals the t-on-precision regression", {
  d <- toy()
  f <- fit_uwls(d$yi, d$sei)
  lmfit <- stats::lm(I(d$yi / d$sei) ~ 0 + I(1 / d$sei))
  sm <- summary(lmfit)
  expect_equal(f$estimate, unname(stats::coef(lmfit)), tolerance = 1e-10)
  expect_equal(f$se, unname(sm$coefficients[1, 2]), tolerance = 1e-10)
  # phi_hat = lm residual variance = Q / (K - 1)
  expect_equal(unname(f$scale["phi"]), sm$sigma^2, tolerance = 1e-10)
  q <- sum((d$yi - f$estimate)^2 / d$sei^2)
  expect_equal(unname(f$scale["phi"]), q / (length(d$yi) - 1),
               tolerance = 1e-10)
  # UWLS point estimate == FE point estimate
  expect_equal(f$estimate, fit_fe(d$yi, d$sei)$estimate, tolerance = 1e-12)
  # yi-scale log-likelihood = lm t-scale log-likelihood - Jacobian
  expect_equal(f$logLik,
               as.numeric(stats::logLik(lmfit)) - sum(log(d$sei)),
               tolerance = 1e-8)
})

test_that("benchmark_models tabulates comparable AIC/BIC", {
  d <- toy()
  tab <- benchmark_models(d$yi, d$sei)
  expect_equal(nrow(tab), 3L)
  expect_true(all(is.finite(tab$AIC)), all(is.finite(tab$BIC)))
  k <- length(d$yi)
  expect_equal(tab$AIC, -2 * tab$logLik + 2 * c(1, 2, 2))
  expect_equal(tab$BIC, -2 * tab$logLik + log(k) * c(1, 2, 2))
  expect_equal(sum(tab$aic_best), 1L)
})

test_that("loo_predict FE matches a hand computation", {
  d <- toy()
  sc <- loo_predict(d$yi, d$sei, model = "fe")
  i <- 1L
  f <- fit_fe(d$yi[-i], d$sei[-i])
  expect_equal(sc$pred_mean[i], f$estimate, tolerance = 1e-12)
  expect_equal(sc$pred_sd[i], sqrt(d$sei[i]^2 + f$se^2), tolerance = 1e-12)
  expect_equal(sc$log_score[i],
               stats::dnorm(d$yi[i], f$estimate,
                            sqrt(d$sei[i]^2 + f$se^2), log = TRUE),
               tolerance = 1e-12)
  expect_equal(sc$pit, stats::pnorm(sc$z))
})

test_that("loo_predict keeps scoring when a refit fails", {
  skip_if_not("local_mocked_bindings" %in% getNamespaceExports("testthat"))
  testthat::local_mocked_bindings(
    fit_re = function(...) stop("synthetic refit failure"),
    .package = "gtmeta"
  )
  d <- toy()
  expect_warning(
    sc <- loo_predict(d$yi, d$sei, model = "re"),
    "LOO refit failed"
  )
  expect_equal(nrow(sc), length(d$yi))
  expect_true(all(is.na(sc$pred_mean)))
  expect_true(all(is.na(sc$log_score)))
})

test_that("loo_compare returns one row per model with a winner", {
  d <- toy()
  tab <- loo_compare(d$yi, d$sei)
  expect_equal(tab$model, c("FE", "RE", "UWLS"))
  expect_true(all(is.finite(tab$log_score)))
  expect_equal(sum(tab$loo_best), 1L)
})

test_that("loo_compare chooses among finite LOO scores when one model fails", {
  skip_if_not("local_mocked_bindings" %in% getNamespaceExports("testthat"))
  testthat::local_mocked_bindings(
    fit_re = function(...) stop("synthetic refit failure"),
    .package = "gtmeta"
  )
  d <- toy()
  expect_warning(tab <- loo_compare(d$yi, d$sei), "LOO refit failed")
  expect_equal(tab$model, c("FE", "RE", "UWLS"))
  expect_true(is.na(tab$log_score[tab$model == "RE"]))
  expect_equal(sum(tab$loo_best, na.rm = TRUE), 1L)
  expect_false(tab$loo_best[tab$model == "RE"])
})

test_that("UWLS t-interval has exact coverage (igmi Prop 5.16, m = 1)", {
  # Proportional total variances + precision weights: the pivot
  # (mu_hat - mu) / sqrt(V_F^loc / (K - 1)) is exactly t_{K-1} at every
  # K. fit_uwls's CI is that interval, so its MC coverage must be
  # nominal up to MC error even at K = 5.
  set.seed(42)
  K <- 5L; mu <- 0.3
  sei <- c(0.10, 0.15, 0.20, 0.25, 0.30)   # unequal, tau2 = 0
  nrep <- 4000L
  covered <- logical(nrep)
  for (r in seq_len(nrep)) {
    yi <- mu + stats::rnorm(K, 0, sei)
    f <- fit_uwls(yi, sei)
    covered[r] <- f$ci.lb <= mu && mu <= f$ci.ub
  }
  # binom sd ~ sqrt(.95*.05/4000) ~ 0.0034; allow 4 sd
  expect_equal(mean(covered), 0.95, tolerance = 0.014)

  # equal sei with tau2 > 0: total variances still proportional
  # (equal), so coverage stays exact although weights ignore tau2
  set.seed(43)
  tau <- 0.4
  sei2 <- rep(0.2, K)
  for (r in seq_len(nrep)) {
    theta <- mu + stats::rnorm(K, 0, tau)
    yi <- theta + stats::rnorm(K, 0, sei2)
    f <- fit_uwls(yi, sei2)
    covered[r] <- f$ci.lb <= mu && mu <= f$ci.ub
  }
  expect_equal(mean(covered), 0.95, tolerance = 0.014)
})
