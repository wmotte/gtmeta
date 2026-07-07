# Diagnostic test accuracy constructor + SROC (igmi_notes Section 8.1, 8.3).
# Identities are checked by hand against the closed-form logit-proportion
# formulas (cf. test_data_adapter.R's log-OR check); the Reitsma
# summary-point agreement (Prop 8.8) is opt-in via mada.

test_that("logit-Se/Sp and variances match the closed forms (Def 8.1)", {
  # TP=80 FN=20 | FP=10 TN=90
  sl <- igmi_dta_studies(80, 20, 10, 90, cc = 0)
  expect_length(sl, 1L)
  th <- sl[[1]]$theta
  expect_equal(unname(th[1]), log(80 / 20))          # logit Se = log 4
  expect_equal(unname(th[2]), log(90 / 10))          # logit Sp = log 9
  Sig <- sl[[1]]$Sigma
  expect_equal(Sig[1, 1], 1 / 80 + 1 / 20)           # delta-method var
  expect_equal(Sig[2, 2], 1 / 90 + 1 / 10)
})

test_that("within-study covariance is exactly diagonal (Lem 8.2)", {
  sl <- igmi_dta_studies(c(80, 60), c(20, 15), c(10, 25), c(90, 75),
                         cc = 0)
  for (s in sl) {
    expect_identical(s$Sigma[1, 2], 0)               # exact zero, not ~0
    expect_identical(s$Sigma[2, 1], 0)
  }
})

test_that("continuity correction rescues a zero cell (Rem 8.3)", {
  # FN = 0 would give an infinite logit; +0.5 to all four cells of that
  # study only.
  sl <- suppressMessages(igmi_dta_studies(30, 0, 8, 92, cc = 0.5))
  expect_length(sl, 1L)
  expect_true(all(is.finite(sl[[1]]$theta)))
  expect_equal(unname(sl[[1]]$theta[1]), log(30.5 / 0.5))
  expect_equal(sl[[1]]$Sigma[1, 1], 1 / 30.5 + 1 / 0.5)
  # a study with no zero cell is left untouched at the same call
  sl2 <- suppressMessages(igmi_dta_studies(c(30, 80), c(0, 20),
                                           c(8, 10), c(92, 90), cc = 0.5))
  expect_equal(unname(sl2[[2]]$theta[1]), log(80 / 20))  # uncorrected
})

test_that("igmi_sroc traces the conditional-mean line through the point", {
  tp <- c(80, 60, 90, 45, 30); fn <- c(20, 15, 10, 5, 4)
  fp <- c(10, 25, 5, 12, 8);   tn <- c(90, 75, 95, 88, 92)
  g  <- lapply(igmi_dta_studies(tp, fn, fp, tn),
               function(s) igmi_gaussian(s$theta, s$Sigma))
  w  <- igmi_precision_weights(g)
  b  <- bw_barycenter(g, w)
  Th <- igmi_mm_cov(g, w)$T_hat
  sr <- igmi_sroc(b, Th)

  slope <- Th[1, 2] / Th[2, 2]
  expect_equal(attr(sr, "slope"), slope)
  # every point satisfies the exact line logit Se = mu_S + slope (logit Sp - mu_C)
  mu_S <- b$theta[1]; mu_C <- b$theta[2]
  expect_equal(sr$logit_se, mu_S + slope * (sr$logit_sp - mu_C))
  # summary operating point is the back-transformed barycenter mean
  sp <- attr(sr, "summary_point")
  expect_equal(unname(sp["sens"]), plogis(mu_S))
  expect_equal(unname(sp["spec"]), plogis(mu_C))
  # curve lives in the unit ROC square
  expect_true(all(sr$sens > 0 & sr$sens < 1 & sr$spec > 0 & sr$spec < 1))
})

test_that("igmi_sroc degenerates to a horizontal curve when logit-Sp variance is zero", {
  mu <- c(logitSe = 0.4, logitSp = 1.1)
  T_hat <- matrix(c(0.2, 0, 0, 0), 2)
  expect_warning(
    sr <- igmi_sroc(mu, T_hat, n = 7L),
    "horizontal line"
  )
  expect_equal(attr(sr, "slope"), 0)
  expect_true(all(sr$logit_se == mu[1]))
  expect_equal(attr(sr, "summary_point"),
               stats::setNames(c(plogis(mu[1]), plogis(mu[2])),
                               c("sens", "spec")))
})

test_that("IGMI summary point is in the same region as mada::reitsma", {
  skip_if_not_installed("mada")
  data("AuditC", package = "mada", envir = environment())
  d <- get("AuditC", envir = environment())
  sl <- suppressMessages(igmi_dta_studies(d$TP, d$FN, d$FP, d$TN))
  g  <- lapply(sl, function(s) igmi_gaussian(s$theta, s$Sigma))
  b  <- bw_barycenter(g, igmi_precision_weights(g))

  fit <- mada::reitsma(d)
  co  <- stats::coef(fit)            # (tsens, tfpr): logit-sens, logit-FPR
  # reitsma is REML random-effects; IGMI-BW is trace-precision (FE-type).
  # Prop 8.8 promises coincidence only under *matched* weights, so on a
  # heterogeneous corpus (AuditC) exact agreement is neither expected nor
  # claimed -- this is the usual FE/RE spread. We assert only same-region
  # agreement on the probability scale: the constructor feeds
  # reitsma-compatible data and IGMI lands in a diagnostically sensible
  # operating region. Tight agreement is reported descriptively in emp/.
  sens_r <- stats::plogis(unname(co[1]))
  spec_r <- stats::plogis(-unname(co[2]))          # Sp = 1 - FPR
  expect_equal(stats::plogis(b$theta[1]), sens_r, tolerance = 0.15)
  expect_equal(stats::plogis(b$theta[2]), spec_r, tolerance = 0.15)
})

test_that("NA cells do not crash the constructor when another study has a zero cell (M2)", {
  # study 1 has an NA cell; study 2 has a genuine zero cell (triggers the
  # continuity-correction subscript-assignment that used to choke on NA)
  tp <- c(50, 0); fn <- c(NA, 20); fp <- c(10, 5); tn <- c(90, 40)
  expect_warning(
    sl <- igmi_dta_studies(tp, fn, fp, tn, cc = 0.5),
    "dropped"
  )
  # the NA study is dropped, the zero-cell study is corrected and kept
  expect_length(sl, 1L)
  expect_true(all(is.finite(sl[[1]]$theta)))
})
