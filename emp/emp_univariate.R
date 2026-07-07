# emp_univariate.R -- Phase 5.B, univariate arm (TODO.md)
#
# IGMI beside FE / RE / UWLS on the real cochrane2025rob corpus,
# mirroring the AIC/BIC model-comparison design of Stanley, Ioannidis,
# Maier, Doucouliagos, Otte & Bartos (2023) and extending it with
# leave-one-study-out predictive scoring. IGMI enters twice:
#
#   IGMI-prop  the proportional-covariance IGMI model. For m = 1 its
#              likelihood, point estimate and t_{K-1} interval are
#              *exactly* UWLS (igmi_notes Prop 5.16; pinned in
#              test_benchmark.R) -- reported as the UWLS row and
#              identified, not duplicated. This is the honest scalar
#              reading of "IGMI as third contender".
#   IGMI-WFR   the genuinely new univariate contender: the equal-weight
#              one-atom WFR Frechet mean wfr_pool() (igmi_notes Rem 4.9
#              updates; Phase 4.B verdict: equal Frechet weights, delta
#              tuned to the effect scale). It has no proper likelihood
#              (exp(-rho) is not integrable for a bounded redescending
#              rho), so it joins the LOO arm only, with predictive
#              N(xhat_{-i}, v_i + tau2DL_{-i} + se_{-i}^2) where se is
#              the sandwich standard error (Rem 4.9, Phase-5.B update).
#
# Corpus: one meta-analysis per review (the largest-K outcome_id; ties ->
# first), K >= K_MIN after FE-pooling subgroup rows within (outcome,
# study) -- the same RevMan-conform rule as igmi_study_list(). One MA per
# review avoids pseudo-replication of reviews with many outcomes.
#
# delta rule: delta = mad(yi) (per meta-analysis), i.e. rejection point
# pi * mad ~ 3.1 robust SDs from the pooled center; falls back to sd(yi)
# when the mad degenerates. Fixed over LOO folds (a tuning constant, not
# a fitted parameter).
#
# Output: emp/results/emp_univariate.rds + printed summary tables.
# Env: GTMETA_COCHRANE_DIR (data tree), GTMETA_EMP_MAX (cap #MAs, 0 =
# all), GTMETA_EMP_KMIN (default 5).

for (f in list.files("R", full.names = TRUE)) source(f)
dir.create(file.path("emp", "results"), showWarnings = FALSE)

DIR   <- Sys.getenv("GTMETA_COCHRANE_DIR",
                    file.path(dirname(getwd()), "cochrane2025rob_data"))
K_MIN <- as.integer(Sys.getenv("GTMETA_EMP_KMIN", "5"))
MAXMA <- as.integer(Sys.getenv("GTMETA_EMP_MAX", "0"))

## ---------------- corpus ----------------

cache <- file.path("emp", "results", "tidy_cache.rds")
if (file.exists(cache)) {
  tidy <- readRDS(cache)
} else {
  tidy <- cochrane_tidy(DIR)
  saveRDS(tidy, cache)
}
cat("tidy corpus:", nrow(tidy), "estimate rows\n")

# FE-pool subgroup duplicates within (outcome_id, study), corpus-wide
key <- paste(tidy$outcome_id, tidy$study, sep = "|")
w <- 1 / tidy$sei^2
num <- tapply(w * tidy$yi, key, sum)
den <- tapply(w, key, sum)
first <- !duplicated(key)
dat <- tidy[first, c("cdid", "review", "outcome_id", "study", "measure",
                     "log_scale", "source", "rob2_overall_bias")]
k1 <- key[first]
dat$yi <- as.numeric(num[k1] / den[k1])
dat$sei <- as.numeric(sqrt(1 / den[k1]))
cat("after subgroup pooling:", nrow(dat), "study-by-outcome estimates\n")

# K per outcome; largest outcome per review
ktab <- table(dat$outcome_id)
dat$K <- as.integer(ktab[dat$outcome_id])
dat <- dat[dat$K >= K_MIN, ]
pick <- do.call(rbind, lapply(split(dat, dat$cdid), function(d) {
  oid <- names(sort(table(d$outcome_id), decreasing = TRUE))[1]
  d[d$outcome_id == oid, ]
}))
mas <- split(pick, pick$outcome_id)
# drop degenerate MAs (no spread)
mas <- mas[vapply(mas, function(d) stats::sd(d$yi) > 0, logical(1))]
if (MAXMA > 0L) mas <- mas[seq_len(min(MAXMA, length(mas)))]
cat("meta-analyses (one per review, K >=", K_MIN, "):", length(mas), "\n")

## ---------------- helpers ----------------

tau2_dl <- function(yi, vi) {
  w <- 1 / vi
  S <- sum(w)
  yb <- sum(w * yi) / S
  Q <- sum(w * (yi - yb)^2)
  cc <- S - sum(w^2) / S
  max(0, (Q - (length(yi) - 1)) / cc)
}

wfr_delta <- function(yi) {
  d <- stats::mad(yi)
  if (d <= 0) d <- stats::sd(yi)
  d
}

# LOO for the WFR arm: robust center from the K-1 remaining studies,
# predictive spread = within-study + DL heterogeneity + sandwich se^2
loo_wfr <- function(yi, sei, delta) {
  k <- length(yi)
  pm <- ps <- numeric(k)
  for (i in seq_len(k)) {
    f <- wfr_pool(yi[-i], weights = rep(1, k - 1), delta = delta)
    pm[i] <- f$estimate
    ps[i] <- sqrt(sei[i]^2 + tau2_dl(yi[-i], sei[-i]^2) + f$se^2)
  }
  z <- (yi - pm) / ps
  data.frame(log_score = stats::dnorm(yi, pm, ps, log = TRUE),
             z = z, cover95 = abs(z) <= stats::qnorm(0.975))
}

## ---------------- main loop ----------------

t0 <- proc.time()
rows <- vector("list", length(mas))
for (j in seq_along(mas)) {
  d <- mas[[j]]
  yi <- d$yi; sei <- d$sei; k <- length(yi)

  bm <- tryCatch(benchmark_models(yi, sei), error = function(e) NULL)
  re_fail <- is.null(bm)
  if (re_fail) {
    fits <- list(fit_fe(yi, sei), fit_uwls(yi, sei))
    bm <- do.call(rbind, lapply(fits, function(f) {
      data.frame(model = f$model, estimate = f$estimate, se = f$se,
                 ci.lb = f$ci.lb, ci.ub = f$ci.ub,
                 scale_par = names(f$scale)[1], scale = unname(f$scale)[1],
                 logLik = f$logLik, AIC = -2 * f$logLik + 2 * f$n_par,
                 BIC = -2 * f$logLik + f$n_par * log(k))
    }))
    bm$aic_best <- bm$AIC == min(bm$AIC)
    bm$bic_best <- bm$BIC == min(bm$BIC)
  }

  delta <- wfr_delta(yi)
  wf <- wfr_pool(yi, weights = rep(1, k), delta = delta)

  loo <- tryCatch({
    cls <- lapply(c("fe", "re", "uwls"), function(m)
      loo_predict(yi, sei, model = m))
    names(cls) <- c("fe", "re", "uwls")
    cls$wfr <- loo_wfr(yi, sei, delta)
    cls
  }, error = function(e) NULL)

  get <- function(model, col) {
    i <- match(model, bm$model)
    if (is.na(i)) NA_real_ else bm[[col]][i]
  }
  loo_sum <- function(m, col) if (is.null(loo)) NA_real_ else sum(loo[[m]][[col]])
  loo_zsd <- function(m) if (is.null(loo)) NA_real_ else stats::sd(loo[[m]]$z)
  loo_cov <- function(m) if (is.null(loo)) NA_real_ else mean(loo[[m]]$cover95)

  ls <- c(fe = loo_sum("fe", "log_score"), re = loo_sum("re", "log_score"),
          uwls = loo_sum("uwls", "log_score"), wfr = loo_sum("wfr", "log_score"))
  rob <- d$rob2_overall_bias
  rows[[j]] <- data.frame(
    cdid = d$cdid[1], outcome_id = d$outcome_id[1], K = k,
    source = d$source[1], measure = d$measure[1],
    fe = get("FE", "estimate"), re = get("RE(ML)", "estimate"),
    uwls = get("UWLS", "estimate"), wfr = wf$estimate,
    se_fe = get("FE", "se"), se_re = get("RE(ML)", "se"),
    se_uwls = get("UWLS", "se"), se_wfr = wf$se,
    tau2 = get("RE(ML)", "scale"), phi = get("UWLS", "scale"),
    delta = delta, mass_ratio = wf$mass_ratio,
    n_killed = k - wf$n_active,
    aic_win = if (any(bm$aic_best)) bm$model[bm$aic_best][1] else NA,
    bic_win = if (any(bm$bic_best)) bm$model[bm$bic_best][1] else NA,
    ls_fe = ls["fe"], ls_re = ls["re"], ls_uwls = ls["uwls"],
    ls_wfr = ls["wfr"],
    loo_win = if (all(is.na(ls))) NA_character_ else names(which.max(ls)),
    zsd_fe = loo_zsd("fe"), zsd_re = loo_zsd("re"),
    zsd_uwls = loo_zsd("uwls"), zsd_wfr = loo_zsd("wfr"),
    cov_fe = loo_cov("fe"), cov_re = loo_cov("re"),
    cov_uwls = loo_cov("uwls"), cov_wfr = loo_cov("wfr"),
    re_fail = re_fail,
    rob2_n = sum(!is.na(rob)),
    rob2_high_frac = if (any(!is.na(rob))) mean(rob == "high", na.rm = TRUE)
                     else NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )
  if (j %% 500L == 0L) {
    cat(sprintf("  %d/%d [%.0fs]\n", j, length(mas), (proc.time() - t0)[3]))
  }
}
emp <- do.call(rbind, rows)
cat(sprintf("fits done: %d MAs in %.0fs (%d RE failures)\n",
            nrow(emp), (proc.time() - t0)[3], sum(emp$re_fail)))

## ---------------- summaries ----------------

print_section <- function(x) cat("\n==", x, "==\n")

print_section("AIC / BIC winner shares (2023-paper design; UWLS == IGMI-prop)")
print(round(rbind(AIC = prop.table(table(emp$aic_win)),
                  BIC = prop.table(table(emp$bic_win))), 3))

print_section("LOO winner shares (incl. IGMI-WFR)")
print(round(prop.table(table(emp$loo_win)), 3))

print_section("LOO calibration: median sd of standardized errors (want ~1)")
print(round(vapply(c("zsd_fe", "zsd_re", "zsd_uwls", "zsd_wfr"),
                   function(c) stats::median(emp[[c]], na.rm = TRUE),
                   numeric(1)), 3))

print_section("LOO 95% predictive coverage (mean over MAs)")
print(round(vapply(c("cov_fe", "cov_re", "cov_uwls", "cov_wfr"),
                   function(c) mean(emp[[c]], na.rm = TRUE), numeric(1)), 3))

print_section("Direct 2023 comparison: RE loses the AIC contest to FE/UWLS")
cat("share of MAs where AIC prefers a fixed-center model over RE:",
    round(mean(emp$aic_win != "RE(ML)", na.rm = TRUE), 3), "\n")

print_section("Median LOO log-score differences vs RE (positive = better)")
# medians: a few MAs have degenerate LOO folds (all remaining yi equal ->
# UWLS-ML predictive sd 0 -> log score -Inf); the median is unaffected
cat("MAs with infinite UWLS log score:", sum(is.infinite(emp$ls_uwls)), "\n")
print(round(c(fe = stats::median(emp$ls_fe - emp$ls_re, na.rm = TRUE),
              uwls = stats::median(emp$ls_uwls - emp$ls_re, na.rm = TRUE),
              wfr = stats::median(emp$ls_wfr - emp$ls_re, na.rm = TRUE)), 4))

print_section("IGMI-WFR LOO win share by K tercile")
kt <- cut(emp$K, stats::quantile(emp$K, c(0, 1/3, 2/3, 1)),
          include.lowest = TRUE)
print(round(tapply(emp$loo_win == "wfr", kt, mean, na.rm = TRUE), 3))

print_section("WFR divergence |wfr - re| / se_re")
dv <- abs(emp$wfr - emp$re) / emp$se_re
print(round(stats::quantile(dv, c(.5, .75, .9, .95, .99), na.rm = TRUE), 3))
emp$wfr_diverges <- dv > 1

print_section("Divergence x RoB2 (covered MAs only)")
cov_ma <- !is.na(emp$rob2_high_frac)
cat("MAs with RoB2:", sum(cov_ma), "\n")
if (sum(cov_ma) > 10) {
  print(round(tapply(emp$rob2_high_frac[cov_ma],
                     emp$wfr_diverges[cov_ma], mean, na.rm = TRUE), 3))
  cat("(mean fraction of high-RoB studies, by WFR-divergence flag)\n")
}

print_section("Killed studies")
cat("MAs with >= 1 killed study:", sum(emp$n_killed > 0),
    sprintf("(%.1f%%); mean retained mass %.3f\n",
            100 * mean(emp$n_killed > 0), mean(emp$mass_ratio)))

saveRDS(list(emp = emp, K_MIN = K_MIN, n_ma = nrow(emp)),
        file.path("emp", "results", "emp_univariate.rds"))
cat("\nsaved emp/results/emp_univariate.rds\n")
