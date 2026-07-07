# validate.R -- lean reproducibility artefact for the IGMI manuscript
#
# One always-runnable, dependency-light place that recomputes the
# headline numerical identities of the IGMI framework *live* and prints a
# single PASS/FAIL table with the theorem it is anchored to. Run from the
# repo root:
#
#   RENV_CONFIG_AUTOLOADER_ENABLED=false Rscript validate.R
#
# It sources R/ directly (no installed gtmeta required, exactly as the
# sim/ drivers do) and needs only metafor + stats for the core tiers.
# Exit status is 1 if any check FAILs, so it doubles as a CI gate.
#
# Tier 1  reduction identities        always run (metafor + stats only)
# Tier 2  exact-coverage Monte Carlo  always run (seeded, fast)
# Tier 3  environment-gated checks     SKIP cleanly when deps/files absent
#           - POT round-trip (needs the pinned reticulate virtualenv)
#           - manuscript numbers from committed sim/ and optional emp/ rds
#
# The checks re-express, in one readable place, the assertions that the
# theorem-anchored testthat suite already pins across several files; see
# tests/testthat/{test_scalar_reduction,test_benchmark,
# test_igmi_heterogeneity,test_wfr_limit}.R and ms/figures.R.

if (!dir.exists("R")) {
  stop("run validate.R from the repository root (no R/ directory here)")
}
for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("validate.R needs the 'metafor' package on the library path")
  }
)

# ---- tiny check harness --------------------------------------------------

.results <- list()
.record <- function(tier, name, anchor, status, detail) {
  .results[[length(.results) + 1L]] <<- list(
    tier = tier, name = name, anchor = anchor,
    status = status, detail = detail)
}

# numeric identity: PASS iff max|got - ref| <= tol
num_check <- function(tier, name, anchor, got, ref, tol = 1e-10) {
  d <- suppressWarnings(max(abs(as.numeric(got) - as.numeric(ref))))
  ok <- is.finite(d) && d <= tol
  .record(tier, name, anchor, if (ok) "PASS" else "FAIL",
          sprintf("|d|=%.1e (tol %.0e)", d, tol))
}

# logical predicate with a free-text detail
bool_check <- function(tier, name, anchor, cond, detail = "") {
  .record(tier, name, anchor, if (isTRUE(cond)) "PASS" else "FAIL", detail)
}

skip_check <- function(tier, name, anchor, reason) {
  .record(tier, name, anchor, "SKIP", reason)
}

# ==========================================================================
# Tier 1 -- reduction identities (metafor + stats only)
# ==========================================================================

yi  <- c(0.12, -0.30, 0.25, 0.05, 0.40, -0.10)
sei <- c(0.10,  0.25, 0.15, 0.30, 0.20,  0.12)
K   <- length(yi)
studies <- igmi_studies(yi = yi, sei = sei)
w   <- igmi_precision_weights(studies)

ee <- metafor::rma(yi = yi, sei = sei, method = "EE")
fe <- fit_fe(yi, sei)
bar <- bw_barycenter(studies, w)

# 1. precision weights are the normalised inverse-variance weights
num_check(1, "precision weights = inverse-variance", "Def 1.1",
          w, (1 / sei^2) / sum(1 / sei^2), tol = 1e-12)

# 2. fixed-effect fit reproduces metafor's equal-effects model (point + se)
num_check(1, "fit_fe point = metafor rma(EE)", "Cor 3.6",
          fe$estimate, as.numeric(ee$beta))
num_check(1, "fit_fe se = metafor rma(EE) se", "Cor 3.6",
          fe$se, ee$se)

# 3. the scalar Bures-Wasserstein barycenter mean IS the FE estimate
num_check(1, "scalar BW barycenter mean = FE", "Cor 3.6",
          bar$theta, fe$estimate)

# 4. mean decoupling: theta_bar = sum w_i y_i for arbitrary weights
w_arb <- c(0.30, 0.20, 0.15, 0.15, 0.10, 0.10)
num_check(1, "BW mean decoupling (arbitrary weights)", "Lem 3.3",
          bw_barycenter(studies, w_arb)$theta, sum((w_arb / sum(w_arb)) * yi),
          tol = 1e-12)

# 5. scalar scale rule: sd of the barycenter = weighted mean of the sds
#    (a deterministic design functional, NOT the pooled standard error)
num_check(1, "BW scale rule sd_bar = sum w_i s_i", "Prop 3.5",
          sqrt(bar$Sigma[1, 1]), sum(w * sei))

# 6. UWLS == the no-intercept lm of t-values on precisions; point == FE
uw   <- fit_uwls(yi, sei)
lmf  <- stats::lm(I(yi / sei) ~ 0 + I(1 / sei))
num_check(1, "UWLS point = t-on-precision lm", "Stanley-Doucouliagos",
          uw$estimate, unname(stats::coef(lmf)))
num_check(1, "UWLS se = lm se", "Stanley-Doucouliagos",
          uw$se, unname(summary(lmf)$coefficients[1, 2]))
num_check(1, "UWLS point = FE point", "Cor 3.6",
          uw$estimate, fe$estimate, tol = 1e-12)

# 7-9. heterogeneity calibration against the classical indices
het <- igmi_heterogeneity(studies, w)
Q   <- ee$QE
S   <- sum(1 / sei^2)
num_check(1, "I_F^2 = Higgins-Thompson I^2", "Cor 5.10",
          het$I_F2, max(0, (Q - (K - 1)) / Q))
num_check(1, "V_loc = Q / S", "Cor 5.3",
          het$V_loc, Q / S)
dl <- metafor::rma(yi = yi, sei = sei, method = "DL")
num_check(1, "moment tr(T_hat) = DerSimonian-Laird tau^2", "Cor 5.10",
          het$tr_T_hat, dl$tau2, tol = 1e-8)

# 10. WFR pooling in the delta -> Inf limit collapses onto FE (pure R,
#     no POT): the W_2 limit of the WFR Frechet mean is the FE estimate.
num_check(1, "WFR pool (delta=Inf) = FE", "Thm 4.12 / Cor 3.6",
          wfr_pool(yi, sei, delta = Inf)$estimate, fe$estimate, tol = 1e-12)

# 11. WFR two-atom closed form collapses to the Hellinger anchor when the
#     atoms are co-located (pure R oracle, no POT):
#     4 d^2 (a0 + a1 - 2 sqrt(a0 a1)) = (2 d |sqrt(a0) - sqrt(a1)|)^2
a0 <- 1.0; a1 <- 0.64; dlt <- 0.7
num_check(1, "WFR two-atom oracle -> Hellinger (d=0)", "Prop 4.8",
          wfr_two_atom_oracle(a0, a1, 0, dlt),
          (2 * dlt * abs(sqrt(a0) - sqrt(a1)))^2, tol = 1e-12)

# 12. bivariate BW point estimate is exactly invariant to the plug-in
#     within-study correlation rho (mean decoupling again, Lem 3.3) --
#     the exact rho-invariance that ms/figures.R:138 relies on.
th2  <- rbind(c(0.30, 0.60), c(-0.10, 0.20), c(0.25, 0.50), c(0.05, 0.40))
s1   <- c(0.10, 0.20, 0.15, 0.25); s2 <- c(0.12, 0.18, 0.14, 0.22)
bw_theta_rho <- function(rho) {
  st <- lapply(seq_len(nrow(th2)), function(i) {
    off <- rho * s1[i] * s2[i]
    igmi_gaussian(th2[i, ], matrix(c(s1[i]^2, off, off, s2[i]^2), 2))
  })
  bw_barycenter(st, igmi_precision_weights(st))$theta
}
t0 <- bw_theta_rho(0)
shift_bw <- max(vapply(c(0.3, 0.6), function(r) max(abs(bw_theta_rho(r) - t0)),
                       numeric(1)))
bool_check(1, "bivariate BW point invariant to rho", "Lem 3.3",
           shift_bw == 0, sprintf("max shift = %.1e over rho in {0,.3,.6}",
                                   shift_bw))

# 13. DTA study object (igmi_dta_studies): the within-study covariance is
#     diagonal EXACTLY -- sensitivity and specificity come from disjoint
#     patient groups (igmi_notes Lem 8.2), so unlike the CDSR pairs there
#     is no plug-in rho. logit-Se and its variance also match the closed
#     forms of Def 8.1 (log-odds and 1/TP + 1/FN).
dta <- igmi_dta_studies(80, 20, 10, 90, cc = 0)[[1]]
bool_check(1, "DTA within-study covariance exactly diagonal", "Lem 8.2",
           dta$Sigma[1, 2] == 0 && dta$Sigma[2, 1] == 0,
           "off-diagonal identically zero")
num_check(1, "DTA logit-Se and delta-method variance", "Def 8.1",
          c(dta$theta[1], dta$Sigma[1, 1]),
          c(log(80 / 20), 1 / 80 + 1 / 20))

# 14. Prediction study object (igmi_pred_studies): the trace-precision BW
#     operating point is EXACTLY invariant to the unidentified within-study
#     correlation rho between discrimination and calibration (igmi_notes
#     Cor 8.5) -- the load-bearing rho-invariance case.
pred_theta_rho <- function(rho) {
  st <- lapply(igmi_pred_studies(c(0.72, 0.68, 0.80, 0.75),
                                 c(0.03, 0.04, 0.02, 0.05),
                                 c(1.10, 0.95, 1.20, 1.05),
                                 c(0.08, 0.10, 0.07, 0.12), rho = rho),
               function(s) igmi_gaussian(s$theta, s$Sigma))
  bw_barycenter(st, igmi_precision_weights(st))$theta
}
p0 <- pred_theta_rho(0)
shift_pred <- max(vapply(c(0.3, 0.6),
                         function(r) max(abs(pred_theta_rho(r) - p0)),
                         numeric(1)))
bool_check(1, "prediction BW operating point invariant to rho", "Cor 8.5",
           shift_pred == 0,
           sprintf("max shift = %.1e over rho in {0,.3,.6}", shift_pred))

# ==========================================================================
# Tier 2 -- exact-coverage Monte Carlo (seeded, fast)
# ==========================================================================

# UWLS t_{K-1} interval == the Frechet-scatter pivot interval (Prop 5.16,
# m = 1): under proportional total variances its coverage is exact at
# every K. Seeded MC; nrep overridable via GTMETA_NREP for a longer run.
nrep <- as.integer(Sys.getenv("GTMETA_NREP", "2000"))
set.seed(42)
mu <- 0.3
sei_c <- c(0.10, 0.15, 0.20, 0.25, 0.30)
Kc <- length(sei_c)
covered <- logical(nrep)
for (r in seq_len(nrep)) {
  yr <- mu + stats::rnorm(Kc, 0, sei_c)
  f  <- fit_uwls(yr, sei_c)
  covered[r] <- f$ci.lb <= mu && mu <= f$ci.ub
}
cov_hat <- mean(covered)
# binom sd ~ sqrt(.95*.05/nrep); allow ~4 sd
tol_cov <- 4 * sqrt(0.95 * 0.05 / nrep)
bool_check(2, "UWLS t-interval coverage = 0.95", "Prop 5.16",
           abs(cov_hat - 0.95) <= tol_cov,
           sprintf("cover=%.3f (nrep=%d, tol %.3f)", cov_hat, nrep, tol_cov))

# ==========================================================================
# Tier 3 -- environment-gated checks (SKIP cleanly if unavailable)
# ==========================================================================

# 3a. POT round-trip: the reticulate -> numpy -> POT WFR value reproduces
#     the closed-form two-atom oracle (needs the pinned virtualenv).
if (requireNamespace("reticulate", quietly = TRUE) &&
    isTRUE(suppressMessages(wfr_available()))) {
  got  <- wfr_dist(1.0, 0, 1.0, 0.5, delta = 0.7)$wfr2
  want <- wfr_two_atom_oracle(1.0, 1.0, 0.5, 0.7)
  num_check(3, "POT WFR round-trip = two-atom oracle", "Prop 4.8",
            got, want, tol = 1e-6)
} else {
  skip_check(3, "POT WFR round-trip = two-atom oracle", "Prop 4.8",
             "pinned POT virtualenv not available")
}

# 3b. manuscript breakdown curve (fig_wfr, panel b) from the committed
#     simulation results: FE == UWLS exactly along the curve, and the
#     WFR delta=0.5 estimate rejects a far outlier while FE does not.
f_ct <- "sim/results/sim_contamination.rds"
if (file.exists(f_ct)) {
  b2 <- readRDS(f_ct)$armB2
  num_check(3, "FE == UWLS along breakdown curve", "fig_wfr(b)",
            b2$fe, b2$uwls, tol = 1e-12)
  fe_end  <- b2$fe[nrow(b2)]
  wfr_end <- b2$wfr0.5eq[nrow(b2)]
  bool_check(3, "WFR(0.5) rejects far outlier; FE does not", "Prop 4.8(2)",
             abs(fe_end) > 0.5 && abs(wfr_end) < 0.2,
             sprintf("at shift=%.1f: FE bias=%.2f, WFR=%.2f",
                     b2$shift[nrow(b2)], fe_end, wfr_end))
} else {
  skip_check(3, "FE == UWLS along breakdown curve", "fig_wfr(b)",
             "sim/results/sim_contamination.rds absent")
}

# 3c. univariate corpus LOO winner shares (fig_uni, n = 2445) from the
#     optional empirical results (external Cochrane corpus; gitignored).
f_emp <- "emp/results/emp_univariate.rds"
if (file.exists(f_emp)) {
  u <- readRDS(f_emp)$emp
  meths  <- c("fe", "re", "uwls", "wfr")
  shares <- vapply(meths, function(m) mean(u$loo_win == m, na.rm = TRUE),
                   numeric(1))
  bool_check(3, "corpus LOO winner shares reproduce", "fig_uni",
             all(is.finite(shares)) && abs(sum(shares) - 1) < 1e-8,
             sprintf("FE/RE/UWLS/WFR = %s (n=%d)",
                     paste(sprintf("%.3f", shares), collapse = "/"),
                     nrow(u)))
} else {
  skip_check(3, "corpus LOO winner shares reproduce", "fig_uni",
             "emp/results/emp_univariate.rds absent (set GTMETA_COCHRANE_DIR)")
}

# ==========================================================================
# report
# ==========================================================================

tbl <- do.call(rbind, lapply(.results, function(r)
  data.frame(tier = r$tier, name = r$name, anchor = r$anchor,
             status = r$status, detail = r$detail,
             stringsAsFactors = FALSE)))

mark <- c(PASS = "[PASS]", FAIL = "[FAIL]", SKIP = "[skip]")
cat("\nIGMI reproducibility validation\n")
cat(strrep("=", 78), "\n", sep = "")
for (ti in sort(unique(tbl$tier))) {
  cat(sprintf("\n-- Tier %d %s\n", ti, strrep("-", 66)))
  sub <- tbl[tbl$tier == ti, ]
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %s  %-42s %-18s %s\n",
                mark[[sub$status[i]]], sub$name[i], sub$anchor[i],
                sub$detail[i]))
  }
}

npass <- sum(tbl$status == "PASS")
nfail <- sum(tbl$status == "FAIL")
nskip <- sum(tbl$status == "SKIP")
cat("\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("%d checks: %d PASS, %d FAIL, %d SKIP\n",
            nrow(tbl), npass, nfail, nskip))

if (nfail > 0L) {
  cat("VALIDATION FAILED\n")
  quit(status = 1L, save = "no")
}
cat("ALL IDENTITIES REPRODUCED\n")
quit(status = 0L, save = "no")
