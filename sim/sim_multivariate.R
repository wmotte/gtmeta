# sim_multivariate.R -- Phase 4.B, arms A and C (TODO.md)
#
# Arm A (bias / MSE / coverage, m = 2): the IGMI barycenter mean with
# precision weights against mvmeta REML, under known within-study
# covariances with cross-outcome correlation and three between-study
# regimes. Interval estimators compared:
#   igmi_pivot  Frechet-scatter pivot (igmi_notes Prop 5.16):
#               per-component t_{K-1} intervals and the joint Hotelling
#               F_{m,K-m} ellipsoid -- no estimate of T enters
#   igmi_naive  known-Sigma normal theory at T = 0 (Prop 5.13):
#               exact under homogeneity, expected to undercover under
#               heterogeneity (the "classical FE mistake", on purpose)
#   mvmeta      REML random-effects Wald intervals (the standard tool)
# Truth for all: the grand mean mu (all point estimators are unbiased
# for mu by linearity; coverage is the interesting axis).
#
# Arm C (V_F calibration, m = 2): type-I error of the exact null
# mixture test (Prop 5.11 via vf_null_test), its power in tr(Tau), and
# unbiasedness of the raw trace moment estimator tr_T_hat_raw
# (E[V_loc] = sum w_i(1-w_i) tr Sigma_i + (1 - sum w_i^2) tr T).
#
# Output: sim/results/sim_multivariate.rds + printed summary tables.

source(file.path("sim", "common.R"))

NREP    <- as.integer(Sys.getenv("GTMETA_NREP", "1000"))
NREP_VF <- as.integer(Sys.getenv("GTMETA_NREP_VF", "1000"))
LEVEL   <- 0.95

mu <- c(0.3, 0.1)
m  <- 2L
tau_scen <- list(
  none  = matrix(0, 2, 2),
  indep = diag(0.05, 2),
  corr  = 0.05 * matrix(c(1, 0.5, 0.5, 1), 2)
)

## ---------------- arm A ----------------

run_cell_A <- function(K, Tau, rho_w, nrep, seed) {
  set.seed(seed)
  zq <- stats::qnorm(1 - (1 - LEVEL) / 2)
  tq <- stats::qt(1 - (1 - LEVEL) / 2, df = K - 1)
  fq <- stats::qf(LEVEL, m, K - m)
  cq <- stats::qchisq(LEVEL, m)
  cols <- c("est1", "est2",
            "piv_se1", "piv_se2", "piv_c1", "piv_c2", "piv_joint",
            "nai_se1", "nai_se2", "nai_c1", "nai_c2", "nai_joint",
            "mv_est1", "mv_est2", "mv_se1", "mv_se2",
            "mv_c1", "mv_c2", "mv_joint")
  out <- matrix(NA_real_, nrep, length(cols),
                dimnames = list(NULL, cols))
  for (r in seq_len(nrep)) {
    dat <- sim_multivariate(K, mu, Tau, rho_within = rho_w)
    st  <- studies_from_sim(dat$y, dat$Sigma)
    w   <- igmi_precision_weights(st)
    piv <- igmi_pivot_ci(st, w, level = LEVEL)
    nai <- igmi_location_ci(st, w, level = LEVEL)
    est <- as.numeric(piv$theta_bar)
    d   <- est - mu
    piv_se <- sqrt(diag(piv$Omega) / (K - 1))
    nai_se <- sqrt(diag(nai$Sigma_bar_theta))
    out[r, 1:2]  <- est
    out[r, 3:4]  <- piv_se
    out[r, 5:6]  <- abs(d) <= tq * piv_se
    out[r, 7]    <- ((K - m) / m) *
      drop(t(d) %*% solve(piv$Omega, d)) <= fq
    out[r, 8:9]  <- nai_se
    out[r, 10:11] <- abs(d) <= zq * nai_se
    out[r, 12]   <- drop(t(d) %*% solve(nai$Sigma_bar_theta, d)) <= cq
    Smat <- t(vapply(dat$Sigma, function(S) {
      c(S[1, 1], S[2, 1], S[2, 2])
    }, numeric(3)))
    fit <- tryCatch(
      mvmeta::mvmeta(dat$y, S = Smat, method = "reml"),
      error = function(e) NULL, warning = function(w) NULL
    )
    if (!is.null(fit)) {
      mv_est <- as.numeric(coef(fit))
      mv_se  <- sqrt(diag(vcov(fit)))
      dm     <- mv_est - mu
      out[r, 13:14] <- mv_est
      out[r, 15:16] <- mv_se
      out[r, 17:18] <- abs(dm) <= zq * mv_se
      out[r, 19]    <- drop(t(dm) %*% solve(vcov(fit), dm)) <= cq
    }
  }
  out
}

summarize_cell_A <- function(raw) {
  ok <- !is.na(raw[, "mv_est1"])
  row <- function(tag, est, se, cov, j) {
    tab <- metric_table(est, se,
                        lo = est - 1, hi = est + 1,  # placeholders
                        truth = mu[j])
    data.frame(estimator = tag, component = j,
               bias = tab$bias, emp_se = tab$emp_se,
               se_ratio = tab$se_ratio, rmse = tab$rmse,
               coverage = mean(cov))
  }
  comp <- rbind(
    row("igmi_pivot", raw[, "est1"], raw[, "piv_se1"], raw[, "piv_c1"], 1),
    row("igmi_pivot", raw[, "est2"], raw[, "piv_se2"], raw[, "piv_c2"], 2),
    row("igmi_naive", raw[, "est1"], raw[, "nai_se1"], raw[, "nai_c1"], 1),
    row("igmi_naive", raw[, "est2"], raw[, "nai_se2"], raw[, "nai_c2"], 2),
    row("mvmeta", raw[ok, "mv_est1"], raw[ok, "mv_se1"],
        raw[ok, "mv_c1"], 1),
    row("mvmeta", raw[ok, "mv_est2"], raw[ok, "mv_se2"],
        raw[ok, "mv_c2"], 2)
  )
  joint <- data.frame(
    igmi_pivot_joint = mean(raw[, "piv_joint"]),
    igmi_naive_joint = mean(raw[, "nai_joint"]),
    mvmeta_joint = mean(raw[ok, "mv_joint"]),
    mvmeta_fail = mean(!ok)
  )
  list(components = comp, joint = joint)
}

grid <- expand.grid(K = c(5L, 10L, 25L), tau = names(tau_scen),
                    rho_w = c(0, 0.5), stringsAsFactors = FALSE)

cat("Arm A:", nrow(grid), "cells x", NREP, "reps\n")
armA <- vector("list", nrow(grid))
t0 <- proc.time()
for (i in seq_len(nrow(grid))) {
  raw <- run_cell_A(grid$K[i], tau_scen[[grid$tau[i]]], grid$rho_w[i],
                    NREP, seed = 4000L + i)
  armA[[i]] <- c(list(cell = grid[i, ]), summarize_cell_A(raw))
  cat(sprintf("  cell %2d/%d  K=%2d tau=%-5s rho=%.1f  [%.0fs]\n",
              i, nrow(grid), grid$K[i], grid$tau[i], grid$rho_w[i],
              (proc.time() - t0)[3]))
}

armA_components <- do.call(rbind, lapply(armA, function(a) {
  cbind(a$cell, a$components, row.names = NULL)
}))
armA_joint <- do.call(rbind, lapply(armA, function(a) {
  cbind(a$cell, a$joint, row.names = NULL)
}))

print_section("Arm A: per-component coverage (nominal 0.95)")
cov_tab <- reshape(
  armA_components[armA_components$component == 1,
                  c("K", "tau", "rho_w", "estimator", "coverage")],
  idvar = c("K", "tau", "rho_w"), timevar = "estimator",
  direction = "wide"
)
print(cov_tab, digits = 3, row.names = FALSE)
print_section("Arm A: joint 95% region coverage + mvmeta failure rate")
print(armA_joint, digits = 3, row.names = FALSE)

## ---------------- arm C ----------------

run_cell_C <- function(K, tau2, rho_w, nrep, seed) {
  set.seed(seed)
  Tau <- diag(tau2, m)
  out <- matrix(NA_real_, nrep, 3,
                dimnames = list(NULL, c("tr_T_raw", "p", "I_F2")))
  for (r in seq_len(nrep)) {
    dat <- sim_multivariate(K, mu, Tau, rho_within = rho_w)
    st  <- studies_from_sim(dat$y, dat$Sigma)
    w   <- igmi_precision_weights(st)
    het <- igmi_heterogeneity(st, w)
    p   <- vf_null_test(st, w, V_loc = het$V_loc, nsim = 2e4)$p
    out[r, ] <- c(het$tr_T_hat_raw, p, het$I_F2)
  }
  out
}

tau2_grid <- c(0, 0.02, 0.05, 0.1)
cat("\nArm C:", length(tau2_grid), "cells x", NREP_VF, "reps\n")
armC <- do.call(rbind, lapply(seq_along(tau2_grid), function(i) {
  raw <- run_cell_C(K = 10L, tau2 = tau2_grid[i], rho_w = 0.3,
                    nrep = NREP_VF, seed = 6000L + i)
  data.frame(
    tau2 = tau2_grid[i], tr_T_true = m * tau2_grid[i],
    tr_T_raw_mean = mean(raw[, "tr_T_raw"]),
    tr_T_raw_bias = mean(raw[, "tr_T_raw"]) - m * tau2_grid[i],
    reject_05 = detection_power(raw[, "p"], alpha = 0.05),
    I_F2_mean = mean(raw[, "I_F2"])
  )
}))

print_section("Arm C: V_F calibration (K=10, m=2, rho_w=0.3)")
cat("row tau2=0: reject_05 = empirical type-I error (nominal 0.05)\n")
print(armC, digits = 3, row.names = FALSE)

saveRDS(list(grid = grid, armA_components = armA_components,
             armA_joint = armA_joint, armC = armC,
             nrep = NREP, nrep_vf = NREP_VF, mu = mu,
             tau_scen = tau_scen),
        file.path("sim", "results", "sim_multivariate.rds"))
cat("\nsaved sim/results/sim_multivariate.rds\n")
