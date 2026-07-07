# data_adapter.R -- cochrane2025rob outputs -> tidy estimate-level data
#
# The upstream pipeline (github.com/wmotte/cochrane2025rob) produces four
# flat files consumed here:
#
#   out.03.rm5s/rm5_data.csv.gz   estimate-level rows parsed from RevMan5
#                                 files via meta::read.rm5(); raw arm data
#                                 (event.e, n.e, ..., mean.e, sd.e, ...) plus
#                                 the summary measure `sm`; TE/seTE columns
#                                 are dropped upstream, so effects are
#                                 *recomputed* here with metafor::escalc().
#   out.03.zips/zip_data.csv.gz   estimate-level rows from RevMan-web zip
#                                 exports; no `sm` column, but generic
#                                 inverse-variance columns (GIV.Mean, GIV.SE)
#                                 and/or the displayed effect with its CI
#                                 (Mean, CI.start, CI.end).
#   out.03.zips/RoB2.csv          estimate-level RoB2 judgements.
#   out.03.zips/RoB2b.csv         study-level RoB2 judgements (fallback).
#
# The adapter returns one tidy data frame with columns
#   review, cdid, study, subgroup, outcome_id, network_id, yi, vi, sei,
#   measure, log_scale, yi_source, source, rob2_* (six domains), rob2_level
# where yi is always on the analysis (log where applicable) scale.
# network_id is NA for this corpus: cochrane2025rob is pairwise; looped
# networks are a tracked open item (TODO.md Phase 2).
#
# Verified against the real corpus (2026-07-03): the rm5 and zip corpora
# are disjoint at review level (no double counting on row-bind), and the
# RoB2 tables only cover (a subset of) zip reviews, so rm5 rows carry NA
# RoB2 columns by data availability, not by key mismatch.

ROB2_DOMAINS <- c(
  "bias_arising_from_the_randomization_process",
  "bias_due_to_deviations_from_intended_interventions",
  "bias_due_to_missing_outcome_data",
  "bias_in_measurement_of_the_outcome",
  "bias_in_selection_of_the_reported_result",
  "overall_bias"
)

# summary measures whose analysis scale is the log scale
LOG_MEASURES <- c("OR", "RR", "PETO", "HR", "IRR")

# severity order used to break ties when collapsing estimate-level RoB2
# judgements to study level (more severe wins a tie)
ROB2_SEVERITY <- c(low = 1L, some_concerns = 2L, high = 3L)

# shared join-key normalization for the RoB2 joins: lowercase, trimmed,
# empty/NA mapped to "." so they compare equal across sources
rob2_norm_key <- function(...) {
  parts <- lapply(list(...), function(x) {
    x <- tolower(trimws(as.character(x)))
    x[is.na(x) | x == "" | x == "na"] <- "."
    x
  })
  do.call(paste, c(parts, sep = "|"))
}

#' Canonical Cochrane review id
#'
#' Review identifiers differ per source: rm5 filenames yield
#' `14651858.CD000009.pub4` while zip exports and the RoB2 tables use the
#' bare `CD000009`. This strips the DOI prefix and `.pubN` suffix so ids
#' can be compared across sources.
#'
#' @param x Character vector of review identifiers.
#' @return Character vector of canonical `CDxxxxxx` ids (uppercase).
#' @export
norm_review_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^14651858\\.", "", x)
  sub("\\.PUB\\d+$", "", x)
}

#' Read rm5-derived estimates and recompute effect sizes
#'
#' Reads `rm5_data.csv.gz` as written by `08_proccess_estimates.R` of
#' cochrane2025rob and recomputes `yi`/`vi` from the raw arm data with
#' [metafor::escalc()], using the RevMan summary measure `sm` per row.
#' Upstream drops the RevMan `TE`/`seTE` columns, so rows without usable
#' raw arm data (e.g. generic inverse-variance outcomes) are dropped with
#' a message.
#'
#' @param path Path to `rm5_data.csv.gz`.
#' @return Tidy data frame (one row per study-by-outcome estimate).
#' @export
read_rm5_estimates <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("studlab", "sm", "infile")
  missing_cols <- setdiff(need, names(df))
  if (length(missing_cols) > 0L) {
    stop("rm5 input lacks expected columns: ",
         paste(missing_cols, collapse = ", "))
  }
  for (col in c("event.e", "n.e", "event.c", "n.c",
                "mean.e", "sd.e", "mean.c", "sd.c",
                "comp.no", "outcome.no", "group.no")) {
    if (is.null(df[[col]])) df[[col]] <- NA_real_
  }

  sm <- toupper(trimws(df$sm))
  yi <- rep(NA_real_, nrow(df))
  vi <- rep(NA_real_, nrow(df))

  bin_measure <- c(OR = "OR", RR = "RR", RD = "RD", PETO = "PETO")
  con_measure <- c(MD = "MD", SMD = "SMD")

  # real exports contain rows with zero-sized arms, on which escalc() errors
  pos_n <- !is.na(df$n.e) & df$n.e > 0 & !is.na(df$n.c) & df$n.c > 0
  is_bin <- sm %in% names(bin_measure) & pos_n &
    stats::complete.cases(df$event.e, df$event.c)
  is_con <- sm %in% names(con_measure) & pos_n &
    stats::complete.cases(df$mean.e, df$sd.e, df$mean.c, df$sd.c)

  for (m in names(bin_measure)) {
    idx <- is_bin & sm == m
    if (!any(idx)) next
    es <- metafor::escalc(
      measure = bin_measure[[m]],
      ai = df$event.e[idx], n1i = df$n.e[idx],
      ci = df$event.c[idx], n2i = df$n.c[idx]
    )
    yi[idx] <- as.numeric(es$yi)
    vi[idx] <- as.numeric(es$vi)
  }
  for (m in names(con_measure)) {
    idx <- is_con & sm == m
    if (!any(idx)) next
    es <- metafor::escalc(
      measure = con_measure[[m]],
      m1i = df$mean.e[idx], sd1i = df$sd.e[idx], n1i = df$n.e[idx],
      m2i = df$mean.c[idx], sd2i = df$sd.c[idx], n2i = df$n.c[idx]
    )
    yi[idx] <- as.numeric(es$yi)
    vi[idx] <- as.numeric(es$vi)
  }

  review <- sub("_data\\.rm5$", "", df$infile)
  out <- data.frame(
    review = review,
    cdid = norm_review_id(review),
    study = df$studlab,
    study_year = NA_integer_,
    subgroup = as.character(df$group.no),
    outcome_id = paste(review, df$comp.no, df$outcome.no, sep = ":"),
    network_id = NA_character_,
    yi = yi,
    vi = vi,
    sei = sqrt(vi),
    measure = sm,
    log_scale = sm %in% LOG_MEASURES,
    yi_source = "raw_arms",
    source = "rm5",
    stringsAsFactors = FALSE
  )

  usable <- is.finite(out$yi) & is.finite(out$vi) & out$vi > 0
  if (any(!usable)) {
    message(sum(!usable), " of ", nrow(out),
            " rm5 rows dropped (no usable raw arm data for sm, ",
            "or zero/infinite variance).")
  }
  out[usable, , drop = FALSE]
}

#' Infer the reporting scale of a (Mean, CI) triple
#'
#' RevMan-web zip exports report the displayed effect `Mean` with
#' `CI.start`/`CI.end` but not the summary measure. A Wald interval is
#' symmetric on its analysis scale, so the scale is inferred by comparing
#' the relative asymmetry of the interval around `Mean` on the identity
#' scale versus the log scale; the log scale is only eligible when all
#' three quantities are positive. Returns `"identity"`, `"log"`, or `NA`
#' (undecidable / malformed).
#'
#' @param m,lo,hi Numeric vectors: point estimate and CI bounds.
#' @param tol Maximum relative asymmetry accepted as "symmetric".
#' @return Character vector over `c("identity", "log", NA)`.
#' @export
infer_ci_scale <- function(m, lo, hi, tol = 0.05) {
  n <- length(m)
  out <- rep(NA_character_, n)
  ok <- is.finite(m) & is.finite(lo) & is.finite(hi) & lo < hi
  asym <- function(m, lo, hi) abs((hi - m) - (m - lo)) / (hi - lo)
  a_id <- rep(Inf, n)
  a_id[ok] <- asym(m[ok], lo[ok], hi[ok])
  pos <- ok & m > 0 & lo > 0
  a_log <- rep(Inf, n)
  a_log[pos] <- asym(log(m[pos]), log(lo[pos]), log(hi[pos]))
  out[ok & a_id <= tol & a_id <= a_log] <- "identity"
  out[ok & a_log <= tol & a_log < a_id] <- "log"
  out
}

#' Read zip-derived estimates
#'
#' Reads `zip_data.csv.gz` as written by `08_proccess_estimates.R` of
#' cochrane2025rob. Effects are taken, in order of preference:
#' generic inverse-variance columns (`GIV.Mean`, `GIV.SE`; already on the
#' analysis scale), else reconstructed from the displayed `Mean` +
#' `CI.start`/`CI.end` with the scale inferred by [infer_ci_scale()]
#' (Wald: `sei = (hi - lo) / (2 * qnorm(1 - (1 - level) / 2))` on the
#' inferred scale). Rows with neither are dropped with a message.
#'
#' @param path Path to `zip_data.csv.gz`.
#' @param level Confidence level of the reported CIs (RevMan default 0.95).
#' @return Tidy data frame (one row per study-by-outcome estimate).
#' @export
read_zip_estimates <- function(path, level = 0.95) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (is.null(df$review) || is.null(df$Study)) {
    stop("zip input lacks expected columns: review, Study")
  }
  for (col in c("GIV.Mean", "GIV.SE", "Mean", "CI.start", "CI.end",
                "Analysis.group", "Analysis.number", "Study.year")) {
    if (is.null(df[[col]])) df[[col]] <- NA
  }
  subgroup <- if (is.null(df$Subgroup)) NA_character_ else as.character(df$Subgroup)

  z <- stats::qnorm(1 - (1 - level) / 2)
  n <- nrow(df)
  yi <- rep(NA_real_, n)
  sei <- rep(NA_real_, n)
  yi_source <- rep(NA_character_, n)
  log_scale <- rep(NA, n)

  giv <- is.finite(df$GIV.Mean) & is.finite(df$GIV.SE) & df$GIV.SE > 0
  yi[giv] <- df$GIV.Mean[giv]
  sei[giv] <- df$GIV.SE[giv]
  yi_source[giv] <- "giv"
  # GIV rows enter RevMan on the analysis scale already; whether that scale
  # is a log scale is not recoverable from the export, so leave it NA.

  rest <- !giv
  scale <- infer_ci_scale(df$Mean, df$CI.start, df$CI.end)
  id <- rest & !is.na(scale) & scale == "identity"
  lg <- rest & !is.na(scale) & scale == "log"
  yi[id] <- df$Mean[id]
  sei[id] <- (df$CI.end[id] - df$CI.start[id]) / (2 * z)
  yi_source[id] <- "ci_identity"
  log_scale[id] <- FALSE
  yi[lg] <- log(df$Mean[lg])
  sei[lg] <- (log(df$CI.end[lg]) - log(df$CI.start[lg])) / (2 * z)
  yi_source[lg] <- "ci_log"
  log_scale[lg] <- TRUE

  out <- data.frame(
    review = df$review,
    cdid = norm_review_id(df$review),
    study = df$Study,
    study_year = suppressWarnings(as.integer(df$Study.year)),
    subgroup = subgroup,
    analysis_group = as.character(df$Analysis.group),
    analysis_number = as.character(df$Analysis.number),
    outcome_id = paste(df$review, df$Analysis.group, df$Analysis.number,
                       sep = ":"),
    network_id = NA_character_,
    yi = yi,
    vi = sei^2,
    sei = sei,
    measure = NA_character_,
    log_scale = log_scale,
    yi_source = yi_source,
    source = "zip",
    stringsAsFactors = FALSE
  )

  usable <- is.finite(out$yi) & is.finite(out$sei) & out$sei > 0
  if (any(!usable)) {
    message(sum(!usable), " of ", nrow(out),
            " zip rows dropped (no GIV columns and no symmetric CI ",
            "on either scale).")
  }
  out[usable, , drop = FALSE]
}

#' Read RoB2 judgements (estimate-level plus study-level fallback)
#'
#' Reads `RoB2.csv` (estimate-level) and optionally `RoB2b.csv`
#' (study-level) as written by `07_process_ROB2.R` of cochrane2025rob.
#' Judgement strings are normalised to `"low"`, `"some_concerns"`,
#' `"high"` (anything else, including empty, becomes `NA`).
#'
#' @param path_estimate Path to `RoB2.csv`.
#' @param path_study Optional path to `RoB2b.csv`.
#' @return List with data frames `estimate` and `study` (the latter may be
#'   `NULL`); both carry the six `rob2_*` domain columns.
#' @export
read_rob2 <- function(path_estimate, path_study = NULL) {
  clean <- function(df) {
    found <- intersect(ROB2_DOMAINS, names(df))
    if (length(found) == 0L) {
      stop("RoB2 input lacks all six domain columns")
    }
    for (d in ROB2_DOMAINS) {
      x <- if (is.null(df[[d]])) rep(NA_character_, nrow(df)) else tolower(trimws(df[[d]]))
      x[grepl("^low", x)] <- "low"
      x[grepl("some concern", x, fixed = TRUE)] <- "some_concerns"
      x[grepl("^high", x)] <- "high"
      x[!x %in% c("low", "some_concerns", "high")] <- NA_character_
      df[[paste0("rob2_", d)]] <- x
      df[[d]] <- NULL
    }
    df
  }
  est <- clean(utils::read.csv(path_estimate, stringsAsFactors = FALSE))
  stu <- if (!is.null(path_study)) {
    clean(utils::read.csv(path_study, stringsAsFactors = FALSE))
  }
  list(estimate = est, study = stu)
}

#' Collapse estimate-level RoB2 judgements to study level
#'
#' Aggregates the cleaned estimate-level RoB2 table to one row per
#' `(review, study)`, taking per domain the modal judgement with ties
#' broken toward the more severe category (`low < some_concerns < high`).
#' On the real corpus most studies carry one judgement across their
#' estimates (519 of 660 review-study pairs), so the mode is usually the
#' unanimous value; the severe tie-break is the conservative choice for
#' the rest.
#'
#' @param est_tab Cleaned estimate-level table (element `estimate` of
#'   [read_rob2()]).
#' @return Data frame with a `key` column (`rob2_norm_key()` over
#'   `(cdid, study)`) and the six `rob2_*` domain columns.
#' @keywords internal
collapse_rob2_study <- function(est_tab) {
  rob_cols <- paste0("rob2_", ROB2_DOMAINS)
  key <- rob2_norm_key(norm_review_id(est_tab$Review), est_tab$Study)
  out <- data.frame(key = unique(key), stringsAsFactors = FALSE)
  for (col in rob_cols) {
    val <- tapply(est_tab[[col]], key, function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0L) return(NA_character_)
      tab <- table(x)
      cand <- names(tab)[tab == max(tab)]
      cand[which.max(ROB2_SEVERITY[cand])]
    })
    out[[col]] <- as.character(val[out$key])
  }
  out
}

#' Attach RoB2 judgements to tidy estimates
#'
#' Three-stage join, each stage only filling rows still unmatched:
#' (1) estimate-level on
#' `(review, analysis_group, analysis_number, subgroup, study)` where
#' those keys exist (zip-derived rows); (2) the study-level table on
#' `(review, study)`; (3) the estimate-level table collapsed to study
#' level ([collapse_rob2_study()]) on `(review, study)`. Later stages
#' fill *per domain*: a row whose estimate-level RoB2 row exists but
#' carries empty judgements (30% of matched rows on the real corpus)
#' still gets its missing domains backfilled from the study level.
#' Review ids are canonicalized with [norm_review_id()] on both sides.
#' The `rob2_level` column records the finest stage whose *key* matched
#' (`"estimate"`, `"study"`, `"study_from_estimates"`, or `NA`);
#' individual domain values may come from a coarser stage via backfill.
#'
#' @param estimates Tidy estimates from [read_rm5_estimates()] /
#'   [read_zip_estimates()] (or their row-bind).
#' @param rob List from [read_rob2()].
#' @return `estimates` with six `rob2_*` columns and `rob2_level` appended.
#' @export
attach_rob2 <- function(estimates, rob) {
  rob_cols <- paste0("rob2_", ROB2_DOMAINS)
  for (col in rob_cols) estimates[[col]] <- NA_character_
  estimates$rob2_level <- NA_character_

  fill <- function(estimates, tab, hit, level) {
    has <- !is.na(hit)
    for (col in rob_cols) {
      take <- has & is.na(estimates[[col]])
      estimates[[col]][take] <- tab[[col]][hit[take]]
    }
    new <- has & is.na(estimates$rob2_level)
    estimates$rob2_level[new] <- level
    estimates
  }

  est_tab <- rob$estimate
  if (!is.null(est_tab) &&
      all(c("Review", "Study") %in% names(est_tab)) &&
      all(c("analysis_group", "analysis_number") %in% names(estimates))) {
    # The estimate-level key needs the finer join columns on both sides;
    # a NULL column (e.g. Analysis.group after read.csv dot-conversion)
    # would silently produce zero matches, so require them explicitly and
    # report a skip instead of falling back to study-level unannounced.
    rob_join_cols <- c("Analysis_group", "Analysis_number", "Subgroup")
    est_join_cols <- c("subgroup", "study")
    miss_rob <- setdiff(rob_join_cols, names(est_tab))
    miss_est <- setdiff(est_join_cols, names(estimates))
    if (length(miss_rob) == 0L && length(miss_est) == 0L) {
      key_rob <- rob2_norm_key(norm_review_id(est_tab$Review), est_tab$Analysis_group,
                          est_tab$Analysis_number, est_tab$Subgroup,
                          est_tab$Study)
      key_est <- rob2_norm_key(norm_review_id(estimates$review),
                          estimates$analysis_group,
                          estimates$analysis_number, estimates$subgroup,
                          estimates$study)
      estimates <- fill(estimates, est_tab, match(key_est, key_rob), "estimate")
    } else {
      warning("attach_rob2: estimate-level RoB2 join skipped (missing ",
              "column(s): ",
              paste(c(if (length(miss_rob)) paste0("rob$estimate$", miss_rob),
                      if (length(miss_est)) paste0("estimates$", miss_est)),
                    collapse = ", "),
              "); falling back to study-level RoB2.", call. = FALSE)
    }
  }

  key_est <- rob2_norm_key(norm_review_id(estimates$review), estimates$study)
  stu_tab <- rob$study
  if (!is.null(stu_tab) && all(c("Review", "Study") %in% names(stu_tab))) {
    key_rob <- rob2_norm_key(norm_review_id(stu_tab$Review), stu_tab$Study)
    estimates <- fill(estimates, stu_tab, match(key_est, key_rob), "study")
  }

  if (!is.null(est_tab) && all(c("Review", "Study") %in% names(est_tab))) {
    coll <- collapse_rob2_study(est_tab)
    estimates <- fill(estimates, coll, match(key_est, coll$key),
                      "study_from_estimates")
  }
  estimates
}

#' One-call adapter: cochrane2025rob output directory -> tidy estimates
#'
#' Convenience orchestration of [read_rm5_estimates()],
#' [read_zip_estimates()], [read_rob2()] and [attach_rob2()]. All inputs
#' are optional; whatever is present in `dir` (with the upstream layout,
#' see file header) is read and row-bound on the shared tidy columns.
#'
#' @param dir Root of a cochrane2025rob checkout/output tree.
#' @param level Confidence level of reported CIs in the zip export.
#' @return Tidy estimate-level data frame with RoB2 columns.
#' @export
cochrane_tidy <- function(dir, level = 0.95) {
  p_rm5 <- file.path(dir, "out.03.rm5s", "rm5_data.csv.gz")
  p_zip <- file.path(dir, "out.03.zips", "zip_data.csv.gz")
  p_rob <- file.path(dir, "out.03.zips", "RoB2.csv")
  p_robb <- file.path(dir, "out.03.zips", "RoB2b.csv")

  parts <- list()
  if (file.exists(p_rm5)) parts$rm5 <- read_rm5_estimates(p_rm5)
  if (file.exists(p_zip)) parts$zip <- read_zip_estimates(p_zip, level = level)
  if (length(parts) == 0L) {
    stop("no cochrane2025rob estimate files found under ", dir)
  }

  all_cols <- Reduce(union, lapply(parts, names))
  parts <- lapply(parts, function(df) {
    for (col in setdiff(all_cols, names(df))) df[[col]] <- NA
    df[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, c(parts, list(make.row.names = FALSE)))

  if (file.exists(p_rob)) {
    rob <- read_rob2(p_rob, if (file.exists(p_robb)) p_robb)
    out <- attach_rob2(out, rob)
  }
  out
}

#' Group tidy estimates into per-study Gaussian inputs for IGMI
#'
#' Collects, within one review, the outcomes listed in `outcome_ids` per
#' study into a mean vector and covariance matrix, giving the
#' `N(theta_i, Sigma_i)` study objects of the IGMI notes (igmi_notes
#' section 1). Within-study cross-outcome correlations are not recoverable
#' from the corpus, so `Sigma_i` is built from the reported variances with
#' a plug-in common correlation `rho` (default 0, i.e. diagonal);
#' sensitivity over `rho` is the intended use. Only studies observed on
#' *all* requested outcomes are returned.
#'
#' On the real corpus 14% of rows are subgroup duplicates within
#' `(study, outcome_id)` (RevMan subgroup rows). `duplicates` selects the
#' policy: `"pool"` (default) combines them by fixed-effect
#' inverse-variance pooling - exactly what RevMan does when totalling
#' subgroups of one study; `"first"` keeps the first row; `"error"`
#' refuses.
#'
#' @param tidy Tidy data frame from [cochrane_tidy()].
#' @param review_id Review to extract.
#' @param outcome_ids Character vector (length m >= 1) of `outcome_id`s.
#' @param rho Plug-in within-study correlation between outcomes.
#' @param duplicates Policy for subgroup rows duplicated within
#'   `(study, outcome_id)`: `"pool"`, `"first"`, or `"error"`.
#' @return List of study objects `list(study =, theta =, Sigma =)`.
#' @export
igmi_study_list <- function(tidy, review_id, outcome_ids, rho = 0,
                            duplicates = c("pool", "first", "error")) {
  duplicates <- match.arg(duplicates)
  m <- length(outcome_ids)
  stopifnot(m >= 1L, rho > -1 / max(1, m - 1), rho < 1)
  sub <- tidy[tidy$review == review_id & tidy$outcome_id %in% outcome_ids, ]
  sub <- sub[is.finite(sub$yi) & is.finite(sub$sei) & sub$sei > 0, ]

  key <- paste(sub$study, sub$outcome_id, sep = "|")
  if (anyDuplicated(key)) {
    if (duplicates == "error") {
      stop("duplicated (study, outcome_id) rows (subgroups); set ",
           "`duplicates` to \"pool\" or \"first\"")
    }
    if (duplicates == "first") {
      sub <- sub[!duplicated(key), ]
    } else {
      message(sum(duplicated(key)), " subgroup rows pooled ",
              "(fixed-effect) within (study, outcome_id).")
      w <- 1 / sub$sei^2
      num <- tapply(w * sub$yi, key, sum)
      den <- tapply(w, key, sum)
      keep <- !duplicated(key)
      k1 <- key[keep]
      sub <- sub[keep, ]
      sub$yi <- as.numeric(num[k1] / den[k1])
      sub$sei <- as.numeric(sqrt(1 / den[k1]))
    }
  }
  out <- list()
  for (s in unique(sub$study)) {
    rows <- sub[sub$study == s, ]
    rows <- rows[match(outcome_ids, rows$outcome_id), ]
    if (anyNA(rows$yi)) next
    sei <- rows$sei
    R <- matrix(rho, m, m)
    diag(R) <- 1
    out[[length(out) + 1L]] <- list(
      study = s,
      theta = stats::setNames(rows$yi, outcome_ids),
      Sigma = diag(sei, m) %*% R %*% diag(sei, m)
    )
  }
  out
}
