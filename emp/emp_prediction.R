# emp_prediction.R -- Phase 7, prediction-model-performance arm (TODO.md)
#
# Meta-analysis of external-validation studies reporting discrimination
# (c-statistic) and calibration (observed:expected ratio). This is the
# home where IGMI's rho-invariance is load-bearing: the within-study
# correlation between discrimination and calibration is genuinely unknown
# (Def 8.4), yet the trace-precision BW point estimate is exactly
# invariant to it (Cor 8.5), while mvmeta-REML shifts. IGMI-BW beside
# mvmeta and metamisc::valmeta, with joint leave-one-study-out scoring and
# the rho-sensitivity band (the direct analogue of the multivariate arm's
# rho figure).
#
# Data: metamisc benchmark validation datasets (EuroSCORE). Transforms
# use metamisc::ccalc / oecalc (logit-C, log-O:E), matching valmeta.
# Output: emp/results/emp_prediction.rds + printed summaries.
# Env: GTMETA_PRED_KMIN (default 5).

for (f in list.files("R", full.names = TRUE)) source(f)
dir.create(file.path("emp", "results"), showWarnings = FALSE)
stopifnot(requireNamespace("metamisc", quietly = TRUE))
have_mv <- requireNamespace("mvmeta", quietly = TRUE)

K_MIN <- as.integer(Sys.getenv("GTMETA_PRED_KMIN", "5"))
RHOS  <- c(0, 0.3, 0.6)
LEVEL <- 0.95

ldmvn2 <- function(x, mu, S) {
  d <- x - mu; det2 <- S[1, 1] * S[2, 2] - S[1, 2]^2
  q <- (d[1]^2 * S[2, 2] - 2 * d[1] * d[2] * S[1, 2] + d[2]^2 * S[1, 1]) / det2
  c(ls = -log(2 * pi) - 0.5 * log(det2) - 0.5 * q, d2 = q)
}
fit_mv <- function(Y, svec) {
  if (!have_mv) return(NULL)
  fit <- tryCatch(suppressWarnings(
    mvmeta::mvmeta(Y ~ 1, S = svec, method = "reml")), error = function(e) NULL)
  if (is.null(fit)) fit <- tryCatch(suppressWarnings(
    mvmeta::mvmeta(Y ~ 1, S = svec, method = "ml")), error = function(e) NULL)
  fit
}

## ---------------- prepare transformed data --------------------------

prep <- function(nm) {
  e <- new.env(); utils::data(list = nm, package = "metamisc", envir = e)
  d <- e[[nm]]
  cc <- tryCatch(metamisc::ccalc(
    cstat = d$c.index, cstat.se = d$se.c.index,
    cstat.cilb = d$c.index.95CIl, cstat.ciub = d$c.index.95CIu,
    N = d$n, O = d$n.events, slab = d$Study, g = "log(cstat/(1-cstat))"),
    error = function(e) NULL)
  oe <- tryCatch(metamisc::oecalc(
    O = d$n.events, E = d$e.events, N = d$n, slab = d$Study, g = "log(OE)"),
    error = function(e) NULL)
  if (is.null(cc) || is.null(oe)) return(NULL)
  df <- data.frame(logitC = cc$theta, seC = cc$theta.se,
                   logOE = oe$theta, seOE = oe$theta.se)
  keep <- stats::complete.cases(df) & df$seC > 0 & df$seOE > 0
  structure(df[keep, , drop = FALSE], raw = d[keep, , drop = FALSE])
}

DATASETS <- c("EuroSCORE")

## ---------------- main ----------------

rows <- list(); store <- list()
for (nm in DATASETS) {
  df <- prep(nm); if (is.null(df) || nrow(df) < K_MIN) next
  K <- nrow(df)

  bw_pt <- mv_pt <- matrix(NA_real_, length(RHOS), 2)
  mv_se <- matrix(NA_real_, length(RHOS), 2)
  conv <- logical(length(RHOS))
  for (r in seq_along(RHOS)) {
    rho <- RHOS[r]
    sl_raw <- igmi_pred_studies(df$logitC, df$seC, df$logOE, df$seOE,
                                rho = rho, transform = FALSE)
    sl <- lapply(sl_raw, function(s) igmi_gaussian(s$theta, s$Sigma))
    wgt <- igmi_precision_weights(sl)
    bary <- bw_barycenter(sl, wgt)
    bw_pt[r, ] <- bary$theta; conv[r] <- bary$converged
    Y <- t(vapply(sl, `[[`, numeric(2), "theta"))
    svec <- t(vapply(sl, function(s) c(s$Sigma[1, 1], s$Sigma[1, 2],
                                       s$Sigma[2, 2]), numeric(3)))
    mv <- fit_mv(Y, svec)
    if (!is.null(mv)) {
      mv_pt[r, ] <- as.numeric(stats::coef(mv))
      mv_se[r, ] <- sqrt(diag(stats::vcov(mv)))
    }
  }

  # rho-sensitivity: IGMI-BW should be flat; mvmeta shifts (std by se@rho0)
  bw_shift <- max(abs(sweep(bw_pt, 2, bw_pt[1, ])))
  mv_shift <- if (all(is.na(mv_pt))) NA_real_ else
    max(abs(sweep(mv_pt, 2, mv_pt[1, ])) / matrix(mv_se[1, ], length(RHOS), 2,
                                                  byrow = TRUE), na.rm = TRUE)

  # valmeta marginal on discrimination (REML), on the c-statistic scale,
  # over the same complete-case rows; compare to IGMI pooled C
  draw <- attr(df, "raw")
  vm <- tryCatch(metamisc::valmeta(
    measure = "cstat", cstat = draw$c.index, cstat.se = draw$se.c.index,
    cstat.cilb = draw$c.index.95CIl, cstat.ciub = draw$c.index.95CIu,
    N = draw$n, O = draw$n.events), error = function(e) NULL)
  val_C <- if (is.null(vm)) NA_real_ else as.numeric(vm$est)   # C scale
  bw_C <- stats::plogis(bw_pt[1, 1])

  # joint LOO at rho = 0 (BW vs mvmeta)
  rho0 <- igmi_pred_studies(df$logitC, df$seC, df$logOE, df$seOE,
                            rho = 0, transform = FALSE)
  sl0 <- lapply(rho0, function(s) igmi_gaussian(s$theta, s$Sigma))
  Slist <- lapply(sl0, `[[`, "Sigma")
  Y0 <- t(vapply(sl0, `[[`, numeric(2), "theta"))
  svec0 <- t(vapply(Slist, function(S) c(S[1, 1], S[1, 2], S[2, 2]), numeric(3)))
  ls_bw <- ls_mv <- rep(NA_real_, K)
  for (i in seq_len(K)) {
    sub <- sl0[-i]; wr <- igmi_precision_weights(sub)
    tb <- Reduce(`+`, Map(function(s, wi) wi * s$theta, sub, wr))
    Tr <- igmi_mm_cov(sub, wr)$T_hat
    Vb <- Reduce(`+`, Map(function(s, wi) wi^2 * (s$Sigma + Tr), sub, wr))
    ls_bw[i] <- ldmvn2(Y0[i, ], tb, Slist[[i]] + Tr + Vb)["ls"]
    mvf <- fit_mv(Y0[-i, , drop = FALSE], svec0[-i, , drop = FALSE])
    if (!is.null(mvf)) ls_mv[i] <- ldmvn2(Y0[i, ], as.numeric(stats::coef(mvf)),
      Slist[[i]] + mvf$Psi + stats::vcov(mvf))["ls"]
  }

  store[[nm]] <- list(rhos = RHOS, bw_pt = bw_pt, mv_pt = mv_pt, mv_se = mv_se)
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = nm, K = K,
    bw_logitC = bw_pt[1, 1], bw_logOE = bw_pt[1, 2],
    bw_C = bw_C, val_C = val_C,
    bw_shift = bw_shift, mv_shift = mv_shift,
    bw_all_converged = all(conv),
    loo_bw = sum(ls_bw, na.rm = TRUE), loo_mv = sum(ls_mv, na.rm = TRUE),
    row.names = NULL, stringsAsFactors = FALSE)
  cat(sprintf("  %-10s K=%d  IGMI C=%.3f  valmeta C=%.3f  BW rho-shift %.2e  mvmeta rho-shift %.2f se\n",
              nm, K, bw_C, val_C, bw_shift, mv_shift))
}
emp <- do.call(rbind, rows)

## ---------------- summaries ----------------

cat("\n== Prediction summary (", nrow(emp), " datasets ) ==\n", sep = "")
cat("BW fixed point all converged (all rho):", all(emp$bw_all_converged), "\n")
cat("rho-sensitivity -- IGMI-BW max shift over rho band (want 0):",
    format(max(emp$bw_shift), scientific = TRUE), "\n")
cat("rho-sensitivity -- mvmeta max shift (in se units of the rho=0 fit):",
    round(max(emp$mv_shift, na.rm = TRUE), 3), "\n")
cat("IGMI-BW vs valmeta pooled C-statistic (median abs diff):",
    round(stats::median(abs(emp$bw_C - emp$val_C), na.rm = TRUE), 4), "\n")
cat("Joint LOO log-score: IGMI-BW", round(sum(emp$loo_bw), 2),
    " mvmeta", round(sum(emp$loo_mv), 2), "\n")

saveRDS(list(emp = emp, store = store, RHOS = RHOS, K_MIN = K_MIN),
        file.path("emp", "results", "emp_prediction.rds"))
cat("\nsaved emp/results/emp_prediction.rds\n")
