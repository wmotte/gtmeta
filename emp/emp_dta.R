# emp_dta.R -- Phase 7, diagnostic test accuracy arm (TODO.md)
#
# The canonical bivariate setting: paired sensitivity/specificity from
# benchmark DTA corpora (mada). IGMI-BW beside the field standard
# (Reitsma bivariate REML) and mvmeta, with joint leave-one-study-out
# predictive scoring and the IGMI summary ROC curve (igmi_notes Def 8.7,
# Prop 8.8). Unlike the generic multivariate arm there is NO plug-in rho:
# the within-study covariance is diagonal exactly (Lemma 8.2).
#
# Data: mada's built-in DTA datasets (each with TP, FN, FP, TN). No
# Cochrane corpus and no network needed.
# Output: emp/results/emp_dta.rds (+ per-dataset SROC curves for the
# manuscript figure) and printed summaries.
# Env: GTMETA_DTA_KMIN (default 5).

for (f in list.files("R", full.names = TRUE)) source(f)
dir.create(file.path("emp", "results"), showWarnings = FALSE)
stopifnot(requireNamespace("mada", quietly = TRUE))
have_mv <- requireNamespace("mvmeta", quietly = TRUE)

K_MIN <- as.integer(Sys.getenv("GTMETA_DTA_KMIN", "5"))
LEVEL <- 0.95

## ---------------- helpers (shared with emp_multivariate.R) -------------

ldmvn2 <- function(x, mu, S) {
  d <- x - mu
  det2 <- S[1, 1] * S[2, 2] - S[1, 2]^2
  q <- (d[1]^2 * S[2, 2] - 2 * d[1] * d[2] * S[1, 2] + d[2]^2 * S[1, 1]) /
    det2
  c(ls = -log(2 * pi) - 0.5 * log(det2) - 0.5 * q, d2 = q)
}

tau2_dl <- function(yi, vi) {
  w <- 1 / vi; S <- sum(w); yb <- sum(w * yi) / S
  Q <- sum(w * (yi - yb)^2); cc <- S - sum(w^2) / S
  max(0, (Q - (length(yi) - 1)) / cc)
}

# independent per-margin predictive (re: DL heterogeneity per outcome)
uv_predict <- function(Y, V, i, re = TRUE) {
  vapply(1:2, function(j) {
    yi <- Y[-i, j]; vi <- V[-i, j]
    t2 <- if (re) tau2_dl(yi, vi) else 0
    wj <- 1 / (vi + t2)
    c(mu = sum(wj * yi) / sum(wj), v = V[i, j] + t2 + 1 / sum(wj))
  }, numeric(2))
}

fit_mv <- function(Y, svec) {
  if (!have_mv) return(NULL)
  fit <- tryCatch(suppressWarnings(
    mvmeta::mvmeta(Y ~ 1, S = svec, method = "reml")), error = function(e) NULL)
  if (is.null(fit)) fit <- tryCatch(suppressWarnings(
    mvmeta::mvmeta(Y ~ 1, S = svec, method = "ml")), error = function(e) NULL)
  fit
}

# BW leave-one-out predictive at study i (matrix moment heterogeneity)
bw_predict <- function(sl, Slist, i) {
  sub <- sl[-i]
  wr <- igmi_precision_weights(sub)
  tb <- Reduce(`+`, Map(function(s, wi) wi * s$theta, sub, wr))
  Tr <- igmi_mm_cov(sub, wr)$T_hat
  Vb <- Reduce(`+`, Map(function(s, wi) wi^2 * (s$Sigma + Tr), sub, wr))
  list(mu = tb, S = Slist[[i]] + Tr + Vb)
}

# Reitsma LOO predictive in its own (logit-sens, logit-FPR) coordinates
reitsma_predict <- function(d, i) {
  fit <- tryCatch(suppressWarnings(mada::reitsma(d[-i, , drop = FALSE])),
                  error = function(e) NULL)
  if (is.null(fit) || is.null(fit$Psi)) return(NULL)
  mu <- as.numeric(stats::coef(fit))                 # (tsens, tfpr)
  # study i observed point and within-study cov on the same scale
  cc <- 0.5
  row <- as.numeric(d[i, c("TP", "FN", "FP", "TN")])
  if (any(row == 0)) row <- row + cc
  TP <- row[1]; FN <- row[2]; FP <- row[3]; TN <- row[4]
  yR <- c(stats::qlogis(TP / (TP + FN)), stats::qlogis(FP / (FP + TN)))
  SiR <- diag(c(1 / TP + 1 / FN, 1 / FP + 1 / TN))
  list(y = yR, mu = mu, S = fit$Psi + SiR)
}

## ---------------- datasets ----------------

get_dta <- function(nm) {
  e <- new.env()
  utils::data(list = nm, package = "mada", envir = e)
  d <- e[[nm]]
  need <- c("TP", "FN", "FP", "TN")
  if (!all(need %in% names(d))) return(NULL)
  d <- d[stats::complete.cases(d[, need]), need, drop = FALSE]
  if (nrow(d) < K_MIN) return(NULL)
  d
}

DATASETS <- c("AuditC", "Dementia", "IAQ", "SAQ", "skin_tests", "smoking")

## ---------------- main loop ----------------

rows <- list(); sroc_store <- list()
for (nm in DATASETS) {
  d <- get_dta(nm); if (is.null(d)) next
  sl_raw <- suppressMessages(
    igmi_dta_studies(d$TP, d$FN, d$FP, d$TN, study = seq_len(nrow(d))))
  if (length(sl_raw) < K_MIN) next
  sl <- lapply(sl_raw, function(s) igmi_gaussian(s$theta, s$Sigma))
  K <- length(sl)
  Y <- t(vapply(sl, `[[`, numeric(2), "theta"))        # (logitSe, logitSp)
  Slist <- lapply(sl, `[[`, "Sigma")
  V <- t(vapply(Slist, diag, numeric(2)))
  svec <- t(vapply(Slist, function(S) c(S[1, 1], S[1, 2], S[2, 2]),
                   numeric(3)))

  # ---- full-data IGMI + comparators ----
  wgt <- igmi_precision_weights(sl)
  bary <- bw_barycenter(sl, wgt)
  Th <- igmi_mm_cov(sl, wgt)$T_hat
  fit_r <- tryCatch(suppressWarnings(mada::reitsma(d)), error = function(e) NULL)
  co_r <- if (is.null(fit_r)) c(NA, NA) else as.numeric(stats::coef(fit_r))
  sens_bw <- stats::plogis(bary$theta[1]); spec_bw <- stats::plogis(bary$theta[2])
  sens_r <- stats::plogis(co_r[1]); spec_r <- stats::plogis(-co_r[2])

  # SROC curves for the figure. Both are the Def 8.7 conditional-mean line
  # (Prop 8.8): IGMI at the moment covariance Th, Reitsma at the REML
  # covariance fit_r$Psi, so they share functional form and differ only
  # through the covariance estimator. mada stores Psi in (tsens, tfpr) =
  # (logit Se, logit FPR = -logit Sp) coordinates; convert to (Se, Sp)
  # order. mada::sroc() defaults to the Rutter-Gatsonis curve (a different
  # functional form) and is kept only as reitsma_sroc_rg for reference.
  reitsma_cm <- if (is.null(fit_r)) NULL else tryCatch({
    Psi <- fit_r$Psi
    mu_r <- c(co_r[1], -co_r[2])                         # (logit Se, logit Sp)
    T_r <- matrix(c(Psi[1, 1], -Psi[1, 2],
                    -Psi[1, 2], Psi[2, 2]), 2, 2)        # cov in (Se, Sp) order
    igmi_sroc(mu_r, T_r)
  }, error = function(e) NULL)
  sroc_store[[nm]] <- list(
    points = data.frame(fpr = 1 - stats::plogis(Y[, 2]),
                        sens = stats::plogis(Y[, 1])),
    igmi = igmi_sroc(bary, Th),
    reitsma_cm = reitsma_cm,
    reitsma_sroc_rg = if (is.null(fit_r)) NULL else
      tryCatch(as.data.frame(mada::sroc(fit_r)), error = function(e) NULL),
    summary_igmi = c(sens = sens_bw, spec = spec_bw),
    summary_reitsma = c(sens = sens_r, spec = spec_r))

  # ---- leave-one-study-out joint predictive ----
  cn <- c("uv_re", "mv", "bw", "reitsma")
  ls_mat <- matrix(NA_real_, K, 4, dimnames = list(NULL, cn)); d2_mat <- ls_mat
  for (i in seq_len(K)) {
    yi <- Y[i, ]
    p <- uv_predict(Y, V, i, re = TRUE)
    ls_mat[i, "uv_re"] <- sum(stats::dnorm(yi, p["mu", ], sqrt(p["v", ]), log = TRUE))
    d2_mat[i, "uv_re"] <- sum((yi - p["mu", ])^2 / p["v", ])

    mvf <- fit_mv(Y[-i, , drop = FALSE], svec[-i, , drop = FALSE])
    if (!is.null(mvf)) {
      Sp <- Slist[[i]] + mvf$Psi + stats::vcov(mvf)
      r <- ldmvn2(yi, as.numeric(stats::coef(mvf)), Sp)
      ls_mat[i, "mv"] <- r["ls"]; d2_mat[i, "mv"] <- r["d2"]
    }
    bp <- bw_predict(sl, Slist, i)
    r <- ldmvn2(yi, bp$mu, bp$S); ls_mat[i, "bw"] <- r["ls"]; d2_mat[i, "bw"] <- r["d2"]

    rp <- reitsma_predict(d, i)
    if (!is.null(rp)) {
      r <- ldmvn2(rp$y, rp$mu, rp$S)
      ls_mat[i, "reitsma"] <- r["ls"]; d2_mat[i, "reitsma"] <- r["d2"]
    }
  }
  ok <- stats::complete.cases(ls_mat)
  ls_cmp <- colSums(ls_mat[ok, , drop = FALSE])
  chi_crit <- stats::qchisq(LEVEL, 2)

  rows[[length(rows) + 1L]] <- data.frame(
    dataset = nm, K = K,
    sens_bw = sens_bw, spec_bw = spec_bw,
    sens_reitsma = sens_r, spec_reitsma = spec_r,
    d_sens = sens_bw - sens_r, d_spec = spec_bw - spec_r,
    sroc_slope = attr(sroc_store[[nm]]$igmi, "slope"),
    bw_converged = bary$converged, bw_iter = bary$iterations,
    loo_win = if (all(is.na(ls_cmp))) NA_character_ else names(which.max(ls_cmp)),
    ls_uv_re = ls_cmp["uv_re"], ls_mv = ls_cmp["mv"],
    ls_bw = ls_cmp["bw"], ls_reitsma = ls_cmp["reitsma"],
    cal_bw = mean(d2_mat[, "bw"]) / 2,
    cal_reitsma = mean(d2_mat[, "reitsma"], na.rm = TRUE) / 2,
    cov_bw = mean(d2_mat[, "bw"] <= chi_crit),
    cov_reitsma = mean(d2_mat[, "reitsma"] <= chi_crit, na.rm = TRUE),
    row.names = NULL, stringsAsFactors = FALSE)
  cat(sprintf("  %-11s K=%2d  IGMI (Se %.3f, Sp %.3f) vs Reitsma (%.3f, %.3f)  LOO win: %s\n",
              nm, K, sens_bw, spec_bw, sens_r, spec_r,
              rows[[length(rows)]]$loo_win))
}
emp <- do.call(rbind, rows)

## ---------------- summaries ----------------

cat("\n== DTA summary (", nrow(emp), " datasets ) ==\n", sep = "")
cat("BW fixed point: all converged:", all(emp$bw_converged),
    "; median iterations:", stats::median(emp$bw_iter), "\n")
cat("IGMI-vs-Reitsma summary point |Delta sens| median:",
    round(stats::median(abs(emp$d_sens)), 3),
    " |Delta spec| median:", round(stats::median(abs(emp$d_spec)), 3), "\n")
cat("LOO joint winner shares:\n"); print(round(prop.table(table(emp$loo_win)), 3))
cat("Median LOO log-score vs Reitsma (positive = IGMI better):",
    round(stats::median(emp$ls_bw - emp$ls_reitsma, na.rm = TRUE), 3), "\n")
cat("Joint 95% predictive coverage: IGMI-BW",
    round(mean(emp$cov_bw), 3), " Reitsma", round(mean(emp$cov_reitsma, na.rm = TRUE), 3), "\n")
cat("Joint calibration mean(d^2)/2 (want ~1): IGMI-BW",
    round(stats::median(emp$cal_bw), 3), " Reitsma",
    round(stats::median(emp$cal_reitsma, na.rm = TRUE), 3), "\n")

saveRDS(list(emp = emp, sroc = sroc_store, K_MIN = K_MIN),
        file.path("emp", "results", "emp_dta.rds"))
cat("\nsaved emp/results/emp_dta.rds\n")
