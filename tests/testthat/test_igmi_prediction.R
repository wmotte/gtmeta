# Prediction-model-performance constructor (igmi_notes Section 8.2).
# The headline property -- exact rho-invariance of the pooled operating
# point (Cor 8.5) -- is asserted bit-for-bit and contrasted with mvmeta's
# rho-sensitivity. Comparator agreement (metamisc::valmeta) is opt-in.

cst <- c(0.72, 0.68, 0.80, 0.75, 0.71, 0.77)
cse <- c(0.03, 0.04, 0.02, 0.05, 0.03, 0.04)
oe  <- c(1.10, 0.95, 1.20, 1.05, 0.90, 1.02)
ose <- c(0.08, 0.10, 0.07, 0.12, 0.09, 0.06)

test_that("delta-method transforms match the closed forms (Def 8.4)", {
  sl <- igmi_pred_studies(cst, cse, oe, ose, rho = 0)
  expect_equal(unname(sl[[1]]$theta[1]), qlogis(cst[1]))
  expect_equal(unname(sl[[1]]$theta[2]), log(oe[1]))
  # se(logit C) = se(C) / (C (1-C)); se(log O:E) = se(O:E) / (O:E)
  expect_equal(sl[[1]]$Sigma[1, 1], (cse[1] / (cst[1] * (1 - cst[1])))^2)
  expect_equal(sl[[1]]$Sigma[2, 2], (ose[1] / oe[1])^2)
  expect_equal(sl[[1]]$Sigma[1, 2], 0)               # rho = 0
})

test_that("rho enters Sigma but not the trace (Cor 8.5 mechanism)", {
  s0 <- igmi_pred_studies(cst, cse, oe, ose, rho = 0)[[1]]$Sigma
  s6 <- igmi_pred_studies(cst, cse, oe, ose, rho = 0.6)[[1]]$Sigma
  expect_gt(abs(s6[1, 2]), 0)                         # off-diagonal moves
  expect_equal(sum(diag(s0)), sum(diag(s6)))          # trace is rho-free
})

test_that("BW pooled operating point is exactly rho-invariant (Cor 8.5)", {
  bw_mean <- function(r) {
    g <- lapply(igmi_pred_studies(cst, cse, oe, ose, rho = r),
                function(s) igmi_gaussian(s$theta, s$Sigma))
    bw_barycenter(g, igmi_precision_weights(g))$theta
  }
  m0 <- bw_mean(0); m3 <- bw_mean(0.3); m6 <- bw_mean(0.6)
  expect_equal(m0, m3, tolerance = 1e-12)
  expect_equal(m0, m6, tolerance = 1e-12)
})

test_that("mvmeta pooled estimate IS rho-sensitive (the contrast)", {
  skip_if_not_installed("mvmeta")
  y <- t(vapply(igmi_pred_studies(cst, cse, oe, ose, rho = 0),
                `[[`, numeric(2), "theta"))
  mv_mean <- function(r) {
    S <- lapply(igmi_pred_studies(cst, cse, oe, ose, rho = r),
                `[[`, "Sigma")
    Slong <- t(vapply(S, function(M) c(M[1, 1], M[1, 2], M[2, 2]),
                      numeric(3)))
    fit <- mvmeta::mvmeta(y, S = Slong, method = "reml")
    unname(stats::coef(fit))
  }
  d <- max(abs(mv_mean(0) - mv_mean(0.6)))
  expect_gt(d, 1e-4)          # REML shifts where IGMI-BW does not
})

test_that("IGMI matches metamisc::valmeta pooled C (opt-in)", {
  skip_if_not_installed("metamisc")
  # valmeta pools the c-statistic on the logit scale by REML; the IGMI
  # marginal on discrimination should agree closely.
  g <- lapply(igmi_pred_studies(cst, cse, oe, ose, rho = 0),
              function(s) igmi_gaussian(s$theta, s$Sigma))
  b <- bw_barycenter(g, igmi_precision_weights(g))
  fit <- metamisc::valmeta(cstat = cst, cstat.se = cse)
  expect_equal(unname(b$theta[1]), qlogis(fit$est), tolerance = 0.1)
})
