# Anchors: igmi_notes Thm 5.2 (decomposition), Cor 5.3 (V_loc = Q/S),
# Cor 5.10 (I_F^2 = I^2, moment estimator = DL exactly),
# Prop 5.11 (multivariate null law), Prop 5.9 (Q ~ chi^2_{K-1}),
# Rem 5.12 Phase-5.B update (matrix moment estimator of T).

test_that("scalar decomposition: V_loc = Q/S, V_scale explicit (Cor 5.3)", {
  yi <- c(0.12, -0.30, 0.25, 0.05, 0.40)
  sei <- c(0.10, 0.25, 0.15, 0.30, 0.20)
  studies <- igmi_studies(yi = yi, sei = sei)
  w <- igmi_precision_weights(studies)
  fv <- frechet_variance(studies, w)

  fit <- metafor::rma(yi = yi, sei = sei, method = "EE")
  S <- sum(1 / sei^2)
  expect_equal(fv$V_loc, fit$QE / S, tolerance = 1e-10)
  sbar <- sum(w * sei)
  expect_equal(fv$V_scale, sum(w * (sei - sbar)^2), tolerance = 1e-10)
  expect_equal(fv$V_F, fv$V_loc + fv$V_scale, tolerance = 1e-12)
})

test_that("I_F^2 = Higgins-Thompson I^2 and tr_T_hat = DL (Cor 5.10)", {
  yi <- c(0.8, -0.4, 0.6, 0.1, -0.2, 0.9)
  sei <- c(0.15, 0.20, 0.10, 0.30, 0.25, 0.18)
  studies <- igmi_studies(yi = yi, sei = sei)
  w <- igmi_precision_weights(studies)
  het <- igmi_heterogeneity(studies, w)

  K <- length(yi)
  fit <- metafor::rma(yi = yi, sei = sei, method = "EE")
  Q <- fit$QE
  expect_equal(het$I_F2, max(0, (Q - (K - 1)) / Q), tolerance = 1e-10)

  fit_dl <- metafor::rma(yi = yi, sei = sei, method = "DL")
  expect_equal(het$tr_T_hat, fit_dl$tau2, tolerance = 1e-8)

  # E_0[V_loc] = (K-1)/S for precision weights (Prop 5.9)
  expect_equal(het$E0_V_loc, (K - 1) / sum(1 / sei^2), tolerance = 1e-12)

  # homogeneous data: I_F2 truncates to 0 and DL truncates identically
  yi0 <- rep(0.2, 4)
  st0 <- igmi_studies(yi = yi0, sei = c(0.2, 0.3, 0.25, 0.22))
  het0 <- igmi_heterogeneity(st0, igmi_precision_weights(st0))
  expect_equal(het0$I_F2, 0)
  expect_equal(het0$tr_T_hat, 0)
  expect_lt(het0$tr_T_hat_raw, 0)
})

test_that("multivariate null matrix: proportional case (Prop 5.11(3))", {
  # Sigma_i = c_i Sigma_0, w_i prop 1/c_i: M(0) = S^{-1}(I - qq') x Sigma_0,
  # eigenvalues sigma_j/S with multiplicity K-1 each
  Sigma0 <- matrix(c(2, 0.5, 0.5, 1), 2)
  ci <- c(1, 2, 4, 0.5)
  K <- length(ci); S <- sum(1 / ci)
  studies <- lapply(ci, function(c) igmi_gaussian(c(0, 0), c * Sigma0))
  w <- (1 / ci) / S
  lam <- vf_null_matrix(studies, w)$lambda
  sig <- eigen(Sigma0, symmetric = TRUE, only.values = TRUE)$values
  ref <- sort(c(rep(sig / S, each = K - 1), rep(0, 2)), decreasing = TRUE)
  expect_equal(lam, ref, tolerance = 1e-10)

  # E_0[V_loc] = tr M(0) = sum w_i (1 - w_i) tr Sigma_i (Prop 5.11(2))
  expect_equal(sum(lam), vf_null_mean(studies, w), tolerance = 1e-10)
})

test_that("null test reproduces the classical Q-test p-value (Prop 5.9)", {
  set.seed(31)
  yi <- c(0.5, -0.1, 0.3, 0.9, -0.4)
  sei <- c(0.2, 0.3, 0.15, 0.25, 0.35)
  studies <- igmi_studies(yi = yi, sei = sei)
  w <- igmi_precision_weights(studies)
  out <- vf_null_test(studies, w, nsim = 2e5)
  Q <- metafor::rma(yi = yi, sei = sei, method = "EE")$QE
  p_ref <- pchisq(Q, df = length(yi) - 1, lower.tail = FALSE)
  expect_gt(p_ref, 0.01)   # keep the comparison in a testable range
  expect_lt(abs(out$p - p_ref), 0.01)
})

test_that("multivariate null law calibrates in Monte Carlo (Prop 5.11(1))", {
  set.seed(32)
  m <- 2; K <- 4
  Sig <- list(diag(c(1, 2)), matrix(c(2, 0.6, 0.6, 1), 2),
              diag(c(0.5, 0.5)), matrix(c(1, -0.3, -0.3, 3), 2))
  w <- c(0.3, 0.3, 0.2, 0.2)
  studies0 <- lapply(Sig, function(S) igmi_gaussian(c(0, 0), S))
  lam <- vf_null_matrix(studies0, w)$lambda

  nrep <- 4000L
  vloc <- numeric(nrep)
  for (r in seq_len(nrep)) {
    th <- lapply(Sig, function(S) {
      as.numeric(gtmeta:::.sqrtm(S) %*% rnorm(m))
    })
    tb <- Reduce(`+`, Map(`*`, th, w))
    vloc[r] <- sum(vapply(seq_len(K), function(i) {
      w[i] * sum((th[[i]] - tb)^2)
    }, numeric(1)))
  }
  # mean matches tr M(0) and the distribution matches the mixture (KS)
  expect_equal(mean(vloc), sum(lam), tolerance = 0.05)
  mix <- as.numeric(matrix(rnorm(nrep * length(lam))^2, nrep) %*% lam)
  ks <- suppressWarnings(stats::ks.test(vloc, mix))
  expect_gt(ks$p.value, 0.01)
})

test_that("matrix moment estimator: trace identity and m=1 DL (Rem 5.12)", {
  # trace of the raw matrix estimator == tr_T_hat_raw, any weights, m = 2
  Sig <- list(diag(c(1, 2)), matrix(c(2, 0.6, 0.6, 1), 2),
              diag(c(0.5, 0.5)), matrix(c(1, -0.3, -0.3, 3), 2))
  th <- list(c(0.5, -0.2), c(-0.1, 0.4), c(0.3, 0.1), c(-0.4, -0.5))
  studies <- Map(igmi_gaussian, th, Sig)
  w <- c(0.3, 0.3, 0.2, 0.2)
  mm <- igmi_mm_cov(studies, w)
  het <- igmi_heterogeneity(studies, w)
  expect_equal(sum(diag(mm$T_raw)), het$tr_T_hat_raw, tolerance = 1e-12)

  # m = 1, precision weights: T_hat == DerSimonian-Laird tau^2 exactly
  yi <- c(0.8, -0.4, 0.6, 0.1, -0.2, 0.9)
  sei <- c(0.15, 0.20, 0.10, 0.30, 0.25, 0.18)
  st1 <- igmi_studies(yi = yi, sei = sei)
  mm1 <- igmi_mm_cov(st1, igmi_precision_weights(st1))
  fit_dl <- metafor::rma(yi = yi, sei = sei, method = "DL")
  expect_equal(mm1$T_hat[1, 1], fit_dl$tau2, tolerance = 1e-8)

  # homogeneous data: raw estimator indefinite, projection clamps to PSD
  st0 <- igmi_studies(yi = rep(0.2, 4), sei = c(0.2, 0.3, 0.25, 0.22))
  mm0 <- igmi_mm_cov(st0, igmi_precision_weights(st0))
  expect_lt(mm0$T_raw[1, 1], 0)
  expect_equal(mm0$T_hat[1, 1], 0)
})

test_that("matrix moment estimator is unbiased before projection (Rem 5.12)", {
  # E[Omega_hat] = sum w_i(1-w_i) Sigma_i + (1 - sum w_i^2) T, so
  # E[T_raw] = T exactly under fixed weights and known Sigma_i
  set.seed(33)
  m <- 2; K <- 6
  Sig <- lapply(seq_len(K), function(i) diag(c(1, 0.5)) * (0.5 + 0.25 * i))
  Tmat <- matrix(c(0.30, 0.12, 0.12, 0.20), 2)
  w <- rep(1 / K, K)
  Th <- gtmeta:::.sqrtm(Tmat)
  nrep <- 6000L
  acc <- matrix(0, m, m)
  for (r in seq_len(nrep)) {
    studies <- lapply(seq_len(K), function(i) {
      eps <- as.numeric(gtmeta:::.sqrtm(Sig[[i]]) %*% rnorm(m)) +
        as.numeric(Th %*% rnorm(m))
      igmi_gaussian(eps, Sig[[i]])
    })
    acc <- acc + igmi_mm_cov(studies, w)$T_raw
  }
  # entrywise absolute MC bound (~3 standard errors at 6000 reps)
  expect_lt(max(abs(acc / nrep - Tmat)), 0.03)
})

test_that("K = 1 and fully concentrated weights give NA, not NaN/crash (M1)", {
  s1 <- list(igmi_gaussian(0.3, matrix(0.2, 1, 1)))
  # heterogeneity is undefined for one study: NA with a warning, no NaN
  expect_warning(h <- igmi_heterogeneity(s1), "K < 2")
  expect_true(is.na(h$I_F2))
  expect_true(is.na(h$tr_T_hat))
  expect_warning(mm <- igmi_mm_cov(s1), "K < 2")
  expect_true(all(is.na(mm$T_hat)))
  expect_warning(vt <- vf_null_test(s1, nsim = 100), "K < 2")
  expect_true(is.na(vt$p))

  # multivariate single study: igmi_mm_cov must not reach eigen() on 0/0
  s1m <- list(igmi_gaussian(c(0.1, -0.2), diag(c(0.3, 0.4))))
  expect_warning(mmm <- igmi_mm_cov(s1m), "K < 2")
  expect_equal(dim(mmm$T_hat), c(2L, 2L))
  expect_true(all(is.na(mmm$T_hat)))
})

test_that("Monte Carlo null p-value is bounded away from exact 0 (2b)", {
  # a wildly heterogeneous configuration would give 0/nsim without the
  # (r + 1)/(nsim + 1) correction; -log(p) must stay finite
  set.seed(7)
  studies <- lapply(1:6, function(i)
    igmi_gaussian(c(-10, 10)[1 + (i %% 2)], matrix(1e-4, 1, 1)))
  out <- vf_null_test(studies, nsim = 500)
  expect_gt(out$p, 0)
  expect_true(is.finite(-log(out$p)))
})
