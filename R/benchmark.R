# benchmark.R -- FE / RE / UWLS benchmark with AIC/BIC and LOO scoring
#
# Mirrors the model-comparison design of Stanley, Ioannidis, Maier,
# Doucouliagos, Otte & Bartos (2023, J Clin Epidemiol 157:53-58): compare
# fixed effect (FE), random effects (RE) and unrestricted weighted least
# squares (UWLS) per meta-analysis by AIC/BIC, extended here with
# leave-one-study-out predictive scoring.
#
# All three log-likelihoods are computed on the *yi scale* with ML
# variance parameters, so AIC/BIC are directly comparable across models:
#   FE:   yi ~ N(mu, vi)             1 parameter
#   RE:   yi ~ N(mu, vi + tau2)      2 parameters (ML via metafor)
#   UWLS: yi ~ N(mu, phi * vi)       2 parameters
# Note that logLik() of the classical UWLS lm-fit (t-values on precisions)
# lives on the t = yi/sei scale and differs by the constant Jacobian
# -sum(log sei); tests pin this identity.
#
# Anchors in math/: the FE estimate is exactly the scalar IGMI
# Bures-Wasserstein barycenter mean under precision weights
# (igmi_notes Cor 3.6) -- test_scalar_reduction.R (Phase 3) points here.

.check_meta_input <- function(yi, sei) {
  stopifnot(is.numeric(yi), is.numeric(sei), length(yi) == length(sei),
            length(yi) >= 2L, all(is.finite(yi)), all(is.finite(sei)),
            all(sei > 0))
}

#' Fixed-effect (common-effect) fit
#'
#' Closed-form inverse-variance pooling: `mu = sum(wi yi) / sum(wi)` with
#' `wi = 1 / sei^2` and `se = 1 / sqrt(sum(wi))`. The log-likelihood is
#' the ML value of the model `yi ~ N(mu, vi)`.
#'
#' @param yi,sei Effect estimates and their standard errors.
#' @param level Confidence level.
#' @return List with `estimate`, `se`, `ci.lb`, `ci.ub`, `logLik`,
#'   `n_par`, `k`, `model`.
#' @export
fit_fe <- function(yi, sei, level = 0.95) {
  .check_meta_input(yi, sei)
  vi <- sei^2
  w <- 1 / vi
  mu <- sum(w * yi) / sum(w)
  se <- 1 / sqrt(sum(w))
  z <- stats::qnorm(1 - (1 - level) / 2)
  ll <- sum(stats::dnorm(yi, mu, sei, log = TRUE))
  list(estimate = mu, se = se,
       ci.lb = mu - z * se, ci.ub = mu + z * se,
       scale = c(phi = 1), logLik = ll, n_par = 1L,
       k = length(yi), model = "FE")
}

#' Random-effects fit (via metafor)
#'
#' Wraps [metafor::rma()] with `method = "ML"` by default so that the
#' returned log-likelihood is comparable with [fit_fe()] and
#' [fit_uwls()] in AIC/BIC. Pass `method = "REML"` for estimation
#' purposes, but do not mix REML likelihoods into the AIC/BIC comparison.
#'
#' @inheritParams fit_fe
#' @param method Heterogeneity estimator handed to [metafor::rma()].
#' @return List as in [fit_fe()], with `scale = c(tau2 = ...)`.
#' @export
fit_re <- function(yi, sei, level = 0.95, method = "ML") {
  .check_meta_input(yi, sei)
  fit <- metafor::rma(yi = yi, sei = sei, method = method,
                      level = 100 * level)
  ll <- sum(stats::dnorm(yi, as.numeric(fit$beta),
                         sqrt(sei^2 + fit$tau2), log = TRUE))
  list(estimate = as.numeric(fit$beta), se = fit$se,
       ci.lb = fit$ci.lb, ci.ub = fit$ci.ub,
       scale = c(tau2 = fit$tau2), logLik = ll, n_par = 2L,
       k = length(yi), model = paste0("RE(", method, ")"))
}

#' Unrestricted weighted least squares (UWLS)
#'
#' The Stanley-Doucouliagos estimator: WLS of `yi` on an intercept with
#' weights `1/vi`, equivalently the no-intercept regression of t-values
#' `yi/sei` on precisions `1/sei`. The point estimate equals the FE
#' estimate; the standard error is the FE standard error times
#' `sqrt(phi_hat)` with `phi_hat = RSS_w / (K - 1)` the unbiased
#' multiplicative dispersion (note `RSS_w = Q`, Cochran's statistic, so
#' `phi_hat = Q / (K - 1)`). The CI uses `t_{K-1}` quantiles, matching
#' the lm-based practice of the 2023 paper. The reported log-likelihood
#' uses the ML dispersion `phi_ml = RSS_w / K` on the yi scale.
#'
#' @inheritParams fit_fe
#' @return List as in [fit_fe()], with `scale = c(phi = phi_hat)`.
#' @export
fit_uwls <- function(yi, sei, level = 0.95) {
  .check_meta_input(yi, sei)
  k <- length(yi)
  vi <- sei^2
  w <- 1 / vi
  mu <- sum(w * yi) / sum(w)
  se_fe <- 1 / sqrt(sum(w))
  rss <- sum(w * (yi - mu)^2)
  phi_hat <- rss / (k - 1)
  phi_ml <- rss / k
  se <- se_fe * sqrt(phi_hat)
  tq <- stats::qt(1 - (1 - level) / 2, df = k - 1)
  ll <- sum(stats::dnorm(yi, mu, sqrt(phi_ml) * sei, log = TRUE))
  list(estimate = mu, se = se,
       ci.lb = mu - tq * se, ci.ub = mu + tq * se,
       scale = c(phi = phi_hat), logLik = ll, n_par = 2L,
       k = k, model = "UWLS")
}

#' Fit FE, RE and UWLS and tabulate AIC/BIC
#'
#' @inheritParams fit_fe
#' @return Data frame with one row per model: estimate, se, CI, scale
#'   parameter, logLik, AIC, BIC and the AIC/BIC winner flags.
#' @export
benchmark_models <- function(yi, sei, level = 0.95) {
  fits <- list(fit_fe(yi, sei, level),
               fit_re(yi, sei, level, method = "ML"),
               fit_uwls(yi, sei, level))
  k <- length(yi)
  out <- do.call(rbind, lapply(fits, function(f) {
    data.frame(model = f$model, estimate = f$estimate, se = f$se,
               ci.lb = f$ci.lb, ci.ub = f$ci.ub,
               scale_par = names(f$scale)[1], scale = unname(f$scale)[1],
               logLik = f$logLik,
               AIC = -2 * f$logLik + 2 * f$n_par,
               BIC = -2 * f$logLik + f$n_par * log(k),
               stringsAsFactors = FALSE)
  }))
  out$aic_best <- out$AIC == min(out$AIC)
  out$bic_best <- out$BIC == min(out$BIC)
  out
}

#' Leave-one-study-out predictive scores for one model
#'
#' For each study i the model is refit on the remaining K-1 studies and
#' the left-out `yi` is scored against its normal predictive distribution:
#'   FE:   N(mu_-i, vi + se_-i^2)
#'   RE:   N(mu_-i, vi + tau2_-i + se_-i^2)
#'   UWLS: N(mu_-i, phi_-i * vi + se_-i^2)
#' Returned per study: predictive mean/sd, log score, raw error,
#' standardized error `z`, PIT value and 95% predictive-interval coverage
#' indicator. Under a well-specified model the PIT values are U(0,1) and
#' `z` is standard normal -- these feed [metric_calibration()].
#'
#' @inheritParams fit_fe
#' @param model One of `"fe"`, `"re"`, `"uwls"`.
#' @param method RE heterogeneity estimator (see [fit_re()]).
#' @return Data frame with K rows.
#' @export
loo_predict <- function(yi, sei, model = c("fe", "re", "uwls"),
                        method = "ML") {
  .check_meta_input(yi, sei)
  model <- match.arg(model)
  stopifnot(length(yi) >= 3L)
  k <- length(yi)
  out <- data.frame(study = seq_len(k), pred_mean = NA_real_,
                    pred_sd = NA_real_)
  refit_failures <- character(0)
  for (i in seq_len(k)) {
    f <- tryCatch(
      switch(model,
        fe = fit_fe(yi[-i], sei[-i]),
        re = fit_re(yi[-i], sei[-i], method = method),
        uwls = fit_uwls(yi[-i], sei[-i])),
      error = function(e) {
        refit_failures <<- c(
          refit_failures,
          paste0("study ", i, " under ", toupper(model), ": ",
                 conditionMessage(e))
        )
        NULL
      }
    )
    if (is.null(f)) next
    extra <- switch(model,
      fe = 0,
      re = unname(f$scale["tau2"]),
      uwls = (unname(f$scale["phi"]) - 1) * sei[i]^2)
    out$pred_mean[i] <- f$estimate
    out$pred_sd[i] <- sqrt(sei[i]^2 + extra + f$se^2)
  }
  if (length(refit_failures) > 0) {
    warning("LOO refit failed for ", length(refit_failures),
            " study/studies (", paste(refit_failures, collapse = "; "), ")")
  }
  out$log_score <- stats::dnorm(yi, out$pred_mean, out$pred_sd, log = TRUE)
  out$err <- yi - out$pred_mean
  out$z <- out$err / out$pred_sd
  out$pit <- stats::pnorm(out$z)
  out$cover95 <- abs(out$z) <= stats::qnorm(0.975)
  out
}

#' Compare FE, RE and UWLS by LOO predictive performance
#'
#' @inheritParams fit_fe
#' @return Data frame with one row per model: summed log score, RMSE of
#'   predictive errors, mean standardized error, sd of standardized
#'   errors (calibration: should be near 1), and 95% predictive coverage.
#' @export
loo_compare <- function(yi, sei) {
  models <- c("fe", "re", "uwls")
  out <- do.call(rbind, lapply(models, function(m) {
    sc <- loo_predict(yi, sei, model = m)
    data.frame(model = toupper(m),
               log_score = sum(sc$log_score),
               rmse = sqrt(mean(sc$err^2)),
               z_mean = mean(sc$z),
               z_sd = stats::sd(sc$z),
               cover95 = mean(sc$cover95),
               stringsAsFactors = FALSE)
  }))
  finite <- is.finite(out$log_score)
  out$loo_best <- FALSE
  if (any(finite)) {
    out$loo_best[finite] <- out$log_score[finite] == max(out$log_score[finite])
  }
  out
}
