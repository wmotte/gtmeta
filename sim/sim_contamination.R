# sim_contamination.R -- Phase 4.B, arm B (TODO.md)
#
# Robustness / breakdown of pooled point estimation under contamination
# (m = 1): FE, RE(REML), UWLS versus the WFR one-atom Frechet mean
# wfr_pool() at delta in {0.5, 1, 2}. Since the scalar balanced BW
# barycenter mean *is* the FE estimate (igmi_notes Cor 3.6), any
# robustness gain of WFR over FE here is exactly the "WFR added value
# over balanced BW" that the TODO decision log asks to demonstrate.
#
# Three sub-studies:
#   B1 main grid   bias/RMSE over eps x understate x tau2 (shift = 2)
#   B2 breakdown   bias as a function of the outlier shift at eps = 0.2:
#                  FE/RE/UWLS are dragged linearly; wfr_pool's influence
#                  redescends to exactly zero once the outliers pass the
#                  transport cutoff pi*delta (Prop 4.8(2))
#   B3 coverage    analytic CIs (FE z, RE Wald, UWLS t) vs a
#                  bootstrap-over-studies percentile CI for wfr_pool
#                  (Rem 5.15(a)), clean vs contaminated
#
# Truth: mu (the uncontaminated grand mean). Contaminated studies have
# their yi shifted and (optionally) their sei understated by a factor.
#
# Output: sim/results/sim_contamination.rds + printed summary tables.

source(file.path("sim", "common.R"))

NREP_MAIN <- as.integer(Sys.getenv("GTMETA_NREP_MAIN", "2000"))
NREP_BRK  <- as.integer(Sys.getenv("GTMETA_NREP_BRK", "1000"))
NREP_COV  <- as.integer(Sys.getenv("GTMETA_NREP_COV", "500"))
B_BOOT    <- as.integer(Sys.getenv("GTMETA_B_BOOT", "400"))

MU <- 0.3
K  <- 10L
DELTAS <- c(0.5, 1, 2)

# REML can fail to converge on grossly contaminated data; fall back to
# the closed-form DerSimonian-Laird fit and count how often
fit_re_robust <- function(yi, sei) {
  out <- tryCatch(fit_re(yi, sei, method = "REML"),
                  error = function(e) NULL)
  if (is.null(out)) {
    out <- fit_re(yi, sei, method = "DL")
    out$fallback <- TRUE
  } else {
    out$fallback <- FALSE
  }
  out
}

# The precision-weighted wfr_pool inherits FE's vulnerability to
# *understated* standard errors (contaminated studies claim 9x weight
# at understate = 3 and can win the mode contest); the equal-weight
# variants isolate that mechanism: location-robustness comes from the
# geometry (rejection point pi*delta), weight-robustness must come
# from the lambda choice.
point_estimates <- function(dat) {
  yi <- dat$yi; sei <- dat$sei
  wfr <- vapply(DELTAS, function(d) {
    wfr_pool(yi, sei, delta = d)$estimate
  }, numeric(1))
  wfr_eq <- vapply(c(0.5, 1), function(d) {
    wfr_pool(yi, weights = rep(1, length(yi)), delta = d)$estimate
  }, numeric(1))
  re <- fit_re_robust(yi, sei)
  c(fe = fit_fe(yi, sei)$estimate,
    re = re$estimate,
    uwls = fit_uwls(yi, sei)$estimate,
    stats::setNames(wfr, paste0("wfr", DELTAS)),
    stats::setNames(wfr_eq, c("wfr0.5eq", "wfr1eq")),
    mass1 = wfr_pool(yi, sei, delta = 1)$mass_ratio,
    re_fb = as.numeric(re$fallback))
}

## ---------------- B1: main grid ----------------

gridB <- expand.grid(tau2 = c(0, 0.05), eps = c(0, 0.1, 0.2, 0.3),
                     understate = c(1, 3))
SHIFT <- 2

cat("Arm B1:", nrow(gridB), "cells x", NREP_MAIN, "reps\n")
t0 <- proc.time()
armB1 <- do.call(rbind, lapply(seq_len(nrow(gridB)), function(i) {
  set.seed(7000L + i)
  est <- t(vapply(seq_len(NREP_MAIN), function(r) {
    dat <- sim_contaminated(K, mu = MU, tau2 = gridB$tau2[i],
                            eps = gridB$eps[i], shift = SHIFT,
                            understate = gridB$understate[i])
    point_estimates(dat)
  }, numeric(10)))
  cat(sprintf("  cell %2d/%d [%.0fs]\n", i, nrow(gridB),
              (proc.time() - t0)[3]))
  nm <- setdiff(colnames(est), c("mass1", "re_fb"))
  do.call(rbind, lapply(nm, function(e) {
    data.frame(gridB[i, ], estimator = e,
               bias = metric_bias(est[, e], MU),
               rmse = sqrt(metric_mse(est[, e], MU)),
               mass1 = mean(est[, "mass1"]),
               re_fb = mean(est[, "re_fb"]), row.names = NULL)
  }))
}))

print_section("Arm B1: bias (shift = 2, K = 10)")
print(reshape(armB1[, c("tau2", "eps", "understate", "estimator", "bias")],
              idvar = c("tau2", "eps", "understate"),
              timevar = "estimator", direction = "wide"),
      digits = 3, row.names = FALSE)
print_section("Arm B1: RMSE")
print(reshape(armB1[, c("tau2", "eps", "understate", "estimator", "rmse")],
              idvar = c("tau2", "eps", "understate"),
              timevar = "estimator", direction = "wide"),
      digits = 3, row.names = FALSE)

## ---------------- B2: breakdown curve ----------------

shift_grid <- seq(0, 6, by = 0.75)
cat("\nArm B2:", length(shift_grid), "shifts x", NREP_BRK, "reps\n")
armB2 <- do.call(rbind, lapply(seq_along(shift_grid), function(i) {
  set.seed(8000L + i)
  est <- t(vapply(seq_len(NREP_BRK), function(r) {
    dat <- sim_contaminated(K, mu = MU, tau2 = 0, eps = 0.2,
                            shift = shift_grid[i], understate = 1)
    point_estimates(dat)
  }, numeric(10)))
  nm <- setdiff(colnames(est), c("mass1", "re_fb"))
  data.frame(shift = shift_grid[i],
             t(vapply(nm, function(e) metric_bias(est[, e], MU),
                      numeric(1))),
             mass1 = mean(est[, "mass1"]), row.names = NULL)
}))

print_section("Arm B2: bias vs outlier shift (eps = 0.2, tau2 = 0)")
cat("wfr rejection points: pi*delta =",
    paste(round(pi * DELTAS, 2), collapse = ", "), "\n")
print(armB2, digits = 3, row.names = FALSE)

## ---------------- B3: coverage ----------------

wfr_boot_ci <- function(yi, sei, delta, B, level = 0.95,
                        equal_weights = FALSE) {
  k <- length(yi)
  draws <- vapply(seq_len(B), function(b) {
    idx <- sample.int(k, k, replace = TRUE)
    if (equal_weights) {
      wfr_pool(yi[idx], weights = rep(1, k), delta = delta)$estimate
    } else {
      wfr_pool(yi[idx], sei[idx], delta = delta)$estimate
    }
  }, numeric(1))
  stats::quantile(draws, c((1 - level) / 2, 1 - (1 - level) / 2),
                  names = FALSE)
}

cat("\nArm B3: coverage, eps in {0, 0.2} x", NREP_COV, "reps, B =",
    B_BOOT, "\n")
armB3 <- do.call(rbind, lapply(c(0, 0.2), function(eps) {
  set.seed(9000L + round(100 * eps))
  cov <- t(vapply(seq_len(NREP_COV), function(r) {
    dat <- sim_contaminated(K, mu = MU, tau2 = 0, eps = eps,
                            shift = SHIFT, understate = 3)
    fe <- fit_fe(dat$yi, dat$sei)
    re <- fit_re_robust(dat$yi, dat$sei)
    uw <- fit_uwls(dat$yi, dat$sei)
    bc <- wfr_boot_ci(dat$yi, dat$sei, delta = 1, B = B_BOOT)
    bce <- wfr_boot_ci(dat$yi, dat$sei, delta = 0.5, B = B_BOOT,
                       equal_weights = TRUE)
    c(fe = fe$ci.lb <= MU && MU <= fe$ci.ub,
      re = re$ci.lb <= MU && MU <= re$ci.ub,
      uwls = uw$ci.lb <= MU && MU <= uw$ci.ub,
      wfr1 = bc[1] <= MU && MU <= bc[2],
      wfr0.5eq = bce[1] <= MU && MU <= bce[2])
  }, logical(5)))
  data.frame(eps = eps, t(colMeans(cov)), row.names = NULL)
}))

print_section("Arm B3: 95% coverage (shift = 2, understate = 3)")
print(armB3, digits = 3, row.names = FALSE)

saveRDS(list(armB1 = armB1, armB2 = armB2, armB3 = armB3,
             mu = MU, K = K, deltas = DELTAS, shift = SHIFT,
             nrep = c(main = NREP_MAIN, brk = NREP_BRK,
                      cov = NREP_COV), B_boot = B_BOOT),
        file.path("sim", "results", "sim_contamination.rds"))
cat("\nsaved sim/results/sim_contamination.rds\n")
