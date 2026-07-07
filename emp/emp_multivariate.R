# emp_multivariate.R -- Phase 5.B, multivariate arm (TODO.md)
#
# The genuine multivariate case: bivariate meta-analyses assembled from
# multi-outcome CDSR reviews, IGMI-BW beside the classical contenders.
# The univariate arm settled the 2023-paper AIC/BIC design; this arm asks
# the multivariate question -- does modelling the joint structure help
# predictively, and how does the exact IGMI inference geometry compare
# with mvmeta's Wald geometry on real data?
#
# Corpus rule (adapter round, TODO Phase 5.B): outcome pairs MUST come
# from different analysis groups -- outcome_ids within one analysis group
# are routinely sensitivity variants of identical data. outcome_id is
# "review:group:number" for both sources, so the group is the middle
# component. One bivariate MA per review: the pair (from different
# groups) with the largest set of studies reporting both outcomes,
# K >= K_MIN complete cases.
#
# Within-study cross-outcome correlation is not recoverable from the
# corpus (adapter header), so it enters as a plug-in rho, run over
# RHOS = {0, 0.3, 0.6} as a sensitivity band; every contender receives
# the same Sigma_i(rho).
#
# Contenders (all leave-one-study-out, joint bivariate predictive):
#   uv_fe   independent per-outcome fixed effect  -- N(fe_-i, v_ij + var)
#   uv_re   independent per-outcome DL random effects
#   mv      mvmeta REML (ML fallback)             -- N(th_-i, Sigma_i +
#           Psi_-i + vcov_-i); the established multivariate benchmark
#   bw      IGMI-BW: barycenter mean = sum w_i theta_i (Lem 3.3) with
#           trace-precision weights, heterogeneity by the matrix moment
#           estimator igmi_mm_cov() (igmi_notes Rem 5.12, Phase-5.B
#           update: E[Omega_hat] = sum w_i(1-w_i) Sigma_i +
#           (1 - sum w_i^2) T), predictive
#           N(theta_bar_-i, Sigma_i + T_-i + sum_j w_j^2 (Sigma_j + T_-i))
#
# Inference geometry (full data, per rho): exact Hotelling region of the
# Frechet-scatter pivot (Prop 5.16) vs the mvmeta Wald ellipse -- area
# ratio; BW fixed-point convergence (Thm 7.5) is verified corpus-wide.
#
# Output: emp/results/emp_multivariate.rds + printed summaries.
# Env: GTMETA_COCHRANE_DIR, GTMETA_EMP_MAX (cap #MAs, 0 = all),
# GTMETA_EMP_KMIN (default 5).

for (f in list.files("R", full.names = TRUE)) source(f)
dir.create(file.path("emp", "results"), showWarnings = FALSE)

DIR   <- Sys.getenv("GTMETA_COCHRANE_DIR",
                    file.path(dirname(getwd()), "cochrane2025rob_data"))
K_MIN <- as.integer(Sys.getenv("GTMETA_EMP_KMIN", "5"))
MAXMA <- as.integer(Sys.getenv("GTMETA_EMP_MAX", "0"))
RHOS  <- c(0, 0.3, 0.6)
LEVEL <- 0.95

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
# (same RevMan-conform rule as igmi_study_list(duplicates = "pool"))
key <- paste(tidy$outcome_id, tidy$study, sep = "|")
w <- 1 / tidy$sei^2
num <- tapply(w * tidy$yi, key, sum)
den <- tapply(w, key, sum)
first <- !duplicated(key)
dat <- tidy[first, c("cdid", "review", "outcome_id", "study", "measure")]
k1 <- key[first]
dat$yi <- as.numeric(num[k1] / den[k1])
dat$sei <- as.numeric(sqrt(1 / den[k1]))
dat <- dat[is.finite(dat$yi) & is.finite(dat$sei) & dat$sei > 0, ]

# analysis group = middle component of "review:group:number"
agrp <- function(oid) sub("^.*:([^:]+):[^:]*$", "\\1", oid)

## ---------------- pair selection: one bivariate MA per review ----------

pick_pair <- function(d) {
  # d: pooled rows of one review; return c(o1, o2) or NULL
  sets <- split(d$study, d$outcome_id)
  sets <- sets[lengths(sets) >= K_MIN]
  if (length(sets) < 2L) return(NULL)
  # bound the pair search in outcome-rich reviews
  if (length(sets) > 40L) {
    sets <- sets[order(lengths(sets), decreasing = TRUE)[1:40]]
  }
  grp <- agrp(names(sets))
  best <- NULL; best_n <- K_MIN - 1L
  ids <- names(sets)
  for (a in seq_along(ids)[-length(ids)]) {
    for (b in seq((a + 1L), length(ids))) {
      if (grp[a] == grp[b]) next
      n <- length(intersect(sets[[a]], sets[[b]]))
      if (n > best_n) { best <- c(ids[a], ids[b]); best_n <- n }
    }
  }
  best
}

cand <- split(dat, dat$cdid)
pairs <- lapply(cand, pick_pair)
pairs <- pairs[!vapply(pairs, is.null, logical(1))]
cat("reviews with an eligible cross-group outcome pair:", length(pairs), "\n")
if (MAXMA > 0L) pairs <- pairs[seq_len(min(MAXMA, length(pairs)))]

## ---------------- helpers ----------------

# bivariate normal log density and Mahalanobis d^2 in closed form
ldmvn2 <- function(x, mu, S) {
  d <- x - mu
  det2 <- S[1, 1] * S[2, 2] - S[1, 2]^2
  q <- (d[1]^2 * S[2, 2] - 2 * d[1] * d[2] * S[1, 2] + d[2]^2 * S[1, 1]) /
    det2
  c(ls = -log(2 * pi) - 0.5 * log(det2) - 0.5 * q, d2 = q)
}

tau2_dl <- function(yi, vi) {
  w <- 1 / vi
  S <- sum(w)
  yb <- sum(w * yi) / S
  Q <- sum(w * (yi - yb)^2)
  cc <- S - sum(w^2) / S
  max(0, (Q - (length(yi) - 1)) / cc)
}

# independent per-outcome predictive (fe: tau2 = 0; re: DL)
uv_predict <- function(Y, V, i, re = TRUE) {
  vapply(1:2, function(j) {
    yi <- Y[-i, j]; vi <- V[-i, j]
    t2 <- if (re) tau2_dl(yi, vi) else 0
    wj <- 1 / (vi + t2)
    est <- sum(wj * yi) / sum(wj)
    c(mu = est, v = V[i, j] + t2 + 1 / sum(wj))
  }, numeric(2))
}

fit_mv <- function(Y, svec) {
  fit <- tryCatch(
    suppressWarnings(mvmeta::mvmeta(Y ~ 1, S = svec, method = "reml")),
    error = function(e) NULL)
  if (is.null(fit)) {
    fit <- tryCatch(
      suppressWarnings(mvmeta::mvmeta(Y ~ 1, S = svec, method = "ml")),
      error = function(e) NULL)
  }
  fit
}

## ---------------- main loop ----------------

t0 <- proc.time()
rows <- list()
for (j in seq_along(pairs)) {
  oids <- pairs[[j]]
  cdid <- names(pairs)[j]
  rev_id <- dat$review[dat$cdid == cdid][1]

  for (rho in RHOS) {
    sl <- tryCatch(
      suppressMessages(igmi_study_list(tidy, rev_id, oids, rho = rho)),
      error = function(e) NULL)
    if (is.null(sl) || length(sl) < K_MIN) next
    sl <- lapply(sl, function(s) igmi_gaussian(s$theta, s$Sigma))
    K <- length(sl)
    Y <- t(vapply(sl, `[[`, numeric(2), "theta"))
    if (any(apply(Y, 2, stats::sd) == 0)) next
    Slist <- lapply(sl, `[[`, "Sigma")
    V <- t(vapply(Slist, diag, numeric(2)))
    svec <- t(vapply(Slist, function(S) c(S[1, 1], S[1, 2], S[2, 2]),
                     numeric(3)))

    # ---- full-data fits ----
    wgt <- igmi_precision_weights(sl)
    piv <- igmi_pivot_ci(sl, wgt, level = LEVEL)
    mm <- igmi_mm_cov(sl, wgt)
    bary <- bw_barycenter(sl, wgt)
    mv <- fit_mv(Y, svec)
    mv_fail <- is.null(mv)

    # exact Hotelling ellipse (Prop 5.16) vs mvmeta Wald ellipse: area and
    # per-component projection half-widths (the area is dominated by the
    # near-rank-1 scatter when the between-study correlation sits at the
    # boundary, so widths are the interpretable comparison)
    c_bw <- (2 / (K - 2)) * piv$f_quantile
    area_bw <- pi * c_bw * sqrt(max(0, det(piv$Omega)))
    hw_bw <- sqrt(c_bw * diag(piv$Omega))
    area_mv <- if (mv_fail) NA_real_ else
      pi * stats::qchisq(LEVEL, 2) * sqrt(max(0, det(stats::vcov(mv))))
    hw_mv <- if (mv_fail) c(NA_real_, NA_real_) else
      sqrt(stats::qchisq(LEVEL, 2) * diag(stats::vcov(mv)))
    th_mv <- if (mv_fail) c(NA_real_, NA_real_) else as.numeric(coef(mv))
    se_mv <- if (mv_fail) c(NA_real_, NA_real_) else
      sqrt(diag(stats::vcov(mv)))
    psi_cor <- if (mv_fail) NA_real_ else {
      P <- mv$Psi
      if (all(diag(P) > 0)) P[1, 2] / sqrt(P[1, 1] * P[2, 2]) else NA_real_
    }
    t_cor <- if (all(diag(mm$T_hat) > 0)) {
      mm$T_hat[1, 2] / sqrt(mm$T_hat[1, 1] * mm$T_hat[2, 2])
    } else NA_real_

    # ---- leave-one-study-out ----
    ls_mat <- matrix(NA_real_, K, 4,
                     dimnames = list(NULL, c("uv_fe", "uv_re", "mv", "bw")))
    d2_mat <- ls_mat
    for (i in seq_len(K)) {
      yi <- Y[i, ]

      for (mod in c("uv_fe", "uv_re")) {
        p <- uv_predict(Y, V, i, re = (mod == "uv_re"))
        z2 <- (yi - p["mu", ])^2 / p["v", ]
        ls_mat[i, mod] <- sum(stats::dnorm(yi, p["mu", ], sqrt(p["v", ]),
                                           log = TRUE))
        d2_mat[i, mod] <- sum(z2)
      }

      mvf <- fit_mv(Y[-i, , drop = FALSE], svec[-i, , drop = FALSE])
      if (!is.null(mvf)) {
        Sp <- Slist[[i]] + mvf$Psi + stats::vcov(mvf)
        r <- ldmvn2(yi, as.numeric(coef(mvf)), Sp)
        ls_mat[i, "mv"] <- r["ls"]; d2_mat[i, "mv"] <- r["d2"]
      }

      sub <- sl[-i]
      wr <- igmi_precision_weights(sub)
      tb <- Reduce(`+`, Map(function(s, wi) wi * s$theta, sub, wr))
      Tr <- igmi_mm_cov(sub, wr)$T_hat
      Vb <- Reduce(`+`, Map(function(s, wi) wi^2 * (s$Sigma + Tr), sub, wr))
      Sp <- Slist[[i]] + Tr + Vb
      r <- ldmvn2(yi, tb, Sp)
      ls_mat[i, "bw"] <- r["ls"]; d2_mat[i, "bw"] <- r["d2"]
    }

    ls <- colSums(ls_mat)
    ok_all <- stats::complete.cases(ls_mat)
    ls_cmp <- colSums(ls_mat[ok_all, , drop = FALSE])   # common folds only
    chi_crit <- stats::qchisq(LEVEL, 2)

    rows[[length(rows) + 1L]] <- data.frame(
      cdid = cdid, o1 = oids[1], o2 = oids[2], K = K, rho = rho,
      bw1 = piv$theta_bar[1], bw2 = piv$theta_bar[2],
      mv1 = th_mv[1], mv2 = th_mv[2],
      se_mv1 = se_mv[1], se_mv2 = se_mv[2],
      max_std_diff = max(abs(piv$theta_bar - th_mv) / se_mv),
      area_bw = area_bw, area_mv = area_mv,
      wr1 = hw_bw[1] / hw_mv[1], wr2 = hw_bw[2] / hw_mv[2],
      t_cor = t_cor, psi_cor = psi_cor,
      tr_T = sum(diag(mm$T_hat)), tr_Psi = if (mv_fail) NA_real_ else
        sum(diag(mv$Psi)),
      bw_converged = bary$converged, bw_iter = bary$iterations,
      mv_fail = mv_fail, mv_loo_fail = sum(is.na(ls_mat[, "mv"])),
      ls_uv_fe = ls["uv_fe"], ls_uv_re = ls["uv_re"],
      ls_mv = ls["mv"], ls_bw = ls["bw"],
      loo_win = if (all(is.na(ls_cmp))) NA_character_ else
        names(which.max(ls_cmp)),
      d2h_uv_fe = mean(d2_mat[, "uv_fe"]) / 2,
      d2h_uv_re = mean(d2_mat[, "uv_re"]) / 2,
      d2h_mv = mean(d2_mat[, "mv"], na.rm = TRUE) / 2,
      d2h_bw = mean(d2_mat[, "bw"]) / 2,
      cov_uv_fe = mean(d2_mat[, "uv_fe"] <= chi_crit),
      cov_uv_re = mean(d2_mat[, "uv_re"] <= chi_crit),
      cov_mv = mean(d2_mat[, "mv"] <= chi_crit, na.rm = TRUE),
      cov_bw = mean(d2_mat[, "bw"] <= chi_crit),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }
  if (j %% 100L == 0L) {
    cat(sprintf("  %d/%d reviews [%.0fs]\n", j, length(pairs),
                (proc.time() - t0)[3]))
  }
}
emp <- do.call(rbind, rows)
cat(sprintf("fits done: %d rows (%d reviews x rho) in %.0fs\n",
            nrow(emp), length(unique(emp$cdid)), (proc.time() - t0)[3]))

## ---------------- summaries ----------------

print_section <- function(x) cat("\n==", x, "==\n")

for (rho in RHOS) {
  e <- emp[emp$rho == rho, ]
  print_section(sprintf("rho = %.1f  (%d bivariate MAs)", rho, nrow(e)))

  cat("LOO winner shares (common folds):\n")
  print(round(prop.table(table(e$loo_win)), 3))

  cat("\nLOO joint calibration, median mean(d^2)/2 (want ~1):\n")
  print(round(vapply(c("d2h_uv_fe", "d2h_uv_re", "d2h_mv", "d2h_bw"),
                     function(c) stats::median(e[[c]], na.rm = TRUE),
                     numeric(1)), 3))

  cat("\nLOO joint 95% predictive coverage (mean over MAs):\n")
  print(round(vapply(c("cov_uv_fe", "cov_uv_re", "cov_mv", "cov_bw"),
                     function(c) mean(e[[c]], na.rm = TRUE), numeric(1)), 3))

  cat("\nMedian LOO log-score difference vs mvmeta (positive = better):\n")
  print(round(c(uv_re = stats::median(e$ls_uv_re - e$ls_mv, na.rm = TRUE),
                bw = stats::median(e$ls_bw - e$ls_mv, na.rm = TRUE)), 4))

  cat("\nCI geometry: median area ratio Hotelling(Prop 5.16) / mvmeta-Wald:",
      round(stats::median(e$area_bw / e$area_mv, na.rm = TRUE), 3),
      "(area is det-dominated; see widths)\n")
  cat("Per-component projection width ratio Hotelling/Wald, median:",
      round(stats::median(c(e$wr1, e$wr2), na.rm = TRUE), 3), "\n")
  cat("Point agreement: median / p90 of max std diff (bw vs mv):",
      round(stats::quantile(e$max_std_diff, c(.5, .9), na.rm = TRUE), 3), "\n")
  bnd_t <- mean(abs(e$t_cor) > 0.99, na.rm = TRUE)
  bnd_p <- mean(abs(e$psi_cor) > 0.99, na.rm = TRUE)
  cat("Between-study correlation at the +/-1 boundary (rank-1 after",
      "truncation): T_hat", round(bnd_t, 3), "; mvmeta Psi",
      round(bnd_p, 3), "\n")
  cat("Interior T_hat correlation, median:",
      round(stats::median(e$t_cor[abs(e$t_cor) <= 0.99], na.rm = TRUE), 3),
      "\n")
  cat("mvmeta failures (full fit):", sum(e$mv_fail),
      "; MAs with failed LOO folds:", sum(e$mv_loo_fail > 0), "\n")
  cat("BW fixed point: all converged:", all(e$bw_converged),
      "; median iterations:", stats::median(e$bw_iter), "\n")
}

print_section("rho sensitivity of the point estimates")
# The IGMI-BW estimate is rho-invariant BY CONSTRUCTION: the mean
# decouples linearly (Lem 3.3) and trace-precision weights depend on
# Sigma_i(rho) only through the rho-free diagonal. mvmeta's GLS weights
# do use the off-diagonal, so its estimate moves with the unknown rho.
by_ma <- split(emp, emp$cdid)
shift_bw <- vapply(by_ma, function(g)
  max(abs(c(diff(range(g$bw1)), diff(range(g$bw2))))), numeric(1))
cat("max abs shift of IGMI-BW over rho band (should be 0):",
    max(shift_bw, na.rm = TRUE), "\n")
shift_mv <- vapply(by_ma, function(g) {
  if (sum(!is.na(g$mv1)) < 2L) return(NA_real_)
  se <- c(g$se_mv1[g$rho == 0][1], g$se_mv2[g$rho == 0][1])
  max(abs(c(diff(range(g$mv1, na.rm = TRUE)),
            diff(range(g$mv2, na.rm = TRUE)))) / se)
}, numeric(1))
cat("mvmeta shift across rho, standardized by se_mv(rho=0):\n")
print(round(stats::quantile(shift_mv, c(.5, .75, .9, .95), na.rm = TRUE), 3))

print_section("LOO winner stability across rho")
win_tab <- tapply(emp$loo_win, list(emp$cdid), function(x)
  length(unique(x[!is.na(x)])))
cat("MAs whose LOO winner is identical at all rho:",
    round(mean(win_tab == 1, na.rm = TRUE), 3), "\n")

saveRDS(list(emp = emp, K_MIN = K_MIN, RHOS = RHOS,
             n_ma = length(unique(emp$cdid))),
        file.path("emp", "results", "emp_multivariate.rds"))
cat("\nsaved emp/results/emp_multivariate.rds\n")
