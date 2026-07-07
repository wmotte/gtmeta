# Fixtures mimic the exact column layout written by cochrane2025rob's
# 08_proccess_estimates.R and 07_process_ROB2.R (see R/data_adapter.R
# header); effect-size recomputation is checked by hand against the
# closed-form log-OR / MD formulas.

make_rm5_fixture <- function(path) {
  df <- data.frame(
    comp.no = c(1, 1, 2, 3),
    outcome.no = c(1, 1, 1, 1),
    group.no = c(1, 1, 1, 1),
    studlab = c("Alpha 2001", "Beta 2003", "Alpha 2001", "Gamma 2010"),
    sm = c("OR", "OR", "MD", "HR"),
    event.e = c(10, 5, NA, NA),
    n.e = c(100, 50, 40, NA),
    event.c = c(20, 8, NA, NA),
    n.c = c(100, 50, 40, NA),
    mean.e = c(NA, NA, 1.2, NA),
    sd.e = c(NA, NA, 0.5, NA),
    mean.c = c(NA, NA, 1.0, NA),
    sd.c = c(NA, NA, 0.4, NA),
    infile = "14651858.CD000001.pub2_data.rm5"
  )
  utils::write.csv(df, path, row.names = FALSE)
  df
}

make_zip_fixture <- function(path) {
  df <- data.frame(
    Analysis.group = c(1, 1, 1, 1, 1),
    Analysis.number = c(1, 1, 2, 3, 1),
    Analysis.name = "outcome",
    Subgroup = NA,
    Study = c("Delta 2015", "Epsilon 2018", "Delta 2015", "Zeta 2020",
              "Eta 2021"),
    Study.year = c(2015, 2018, 2015, 2020, 2021),
    GIV.Mean = c(0.30, NA, NA, NA, 0.10),
    GIV.SE = c(0.10, NA, NA, NA, 0.20),
    Mean = c(NA, 2.0, 1.0, 3.0, NA),
    CI.start = c(NA, 1.0, 0.0, 1.0, NA),
    CI.end = c(NA, 4.0, 2.0, 4.0, NA),
    review = "CD000002"
  )
  utils::write.csv(df, path, row.names = FALSE)
  df
}

make_rob2_fixtures <- function(path_est, path_stu) {
  est <- data.frame(
    Review = "CD000002",
    Analysis_group = c(1, 1, 2), Analysis_number = c(1, 1, 1),
    Subgroup = NA,
    Study = c("Delta 2015", "Eta 2021", "Eta 2021"),
    bias_arising_from_the_randomization_process =
      c("Low risk of bias", "", "High"),
    bias_due_to_deviations_from_intended_interventions =
      c("Some concerns", "", "High"),
    bias_due_to_missing_outcome_data = c("High risk of bias", "", "High"),
    bias_in_measurement_of_the_outcome = c("Low", "", "High"),
    bias_in_selection_of_the_reported_result = c("", "", "High"),
    overall_bias = c("Some concerns", "", "High")
  )
  utils::write.csv(est, path_est, row.names = FALSE)
  stu <- data.frame(
    Review = c("CD000002", "14651858.CD000001.pub2"),
    Study = c("Epsilon 2018", "Alpha 2001"),
    bias_arising_from_the_randomization_process = "High",
    bias_due_to_deviations_from_intended_interventions = "Low",
    bias_due_to_missing_outcome_data = "Low",
    bias_in_measurement_of_the_outcome = "Low",
    bias_in_selection_of_the_reported_result = "Low",
    overall_bias = "High"
  )
  utils::write.csv(stu, path_stu, row.names = FALSE)
  list(est = est, stu = stu)
}

test_that("read_rm5_estimates recomputes effects and drops unusable rows", {
  path <- withr::local_tempfile(fileext = ".csv")
  make_rm5_fixture(path)
  expect_message(out <- read_rm5_estimates(path), "1 of 4")
  expect_equal(nrow(out), 3L)

  # hand-checked log OR for Alpha 2001: a=10, b=90, c=20, d=80
  expect_equal(out$yi[1], log((10 * 80) / (90 * 20)), tolerance = 1e-12)
  expect_equal(out$vi[1], 1 / 10 + 1 / 90 + 1 / 20 + 1 / 80,
               tolerance = 1e-12)
  expect_true(out$log_scale[1])

  # hand-checked MD row
  md <- out[out$measure == "MD", ]
  expect_equal(md$yi, 1.2 - 1.0, tolerance = 1e-12)
  expect_equal(md$vi, 0.5^2 / 40 + 0.4^2 / 40, tolerance = 1e-12)
  expect_false(md$log_scale)

  expect_equal(unique(out$review), "14651858.CD000001.pub2")
  expect_equal(out$outcome_id[1], "14651858.CD000001.pub2:1:1")
  expect_true(all(is.na(out$network_id)))
})

test_that("infer_ci_scale distinguishes identity, log, and junk", {
  expect_equal(infer_ci_scale(2, 1, 4), "log")          # symmetric in log
  expect_equal(infer_ci_scale(1, 0, 2), "identity")     # symmetric, lo = 0
  expect_true(is.na(infer_ci_scale(3, 1, 4)))           # neither
  expect_true(is.na(infer_ci_scale(NA, 1, 4)))
  expect_true(is.na(infer_ci_scale(2, 4, 1)))           # malformed
})

test_that("read_zip_estimates prefers GIV, reconstructs from CI", {
  path <- withr::local_tempfile(fileext = ".csv")
  make_zip_fixture(path)
  expect_message(out <- read_zip_estimates(path), "1 of 5")
  expect_equal(nrow(out), 4L)

  giv <- out[out$yi_source == "giv", ]
  expect_equal(giv$yi, c(0.30, 0.10))
  expect_equal(giv$sei, c(0.10, 0.20))

  z <- qnorm(0.975)
  lg <- out[out$yi_source == "ci_log", ]
  expect_equal(lg$yi, log(2))
  expect_equal(lg$sei, (log(4) - log(1)) / (2 * z), tolerance = 1e-12)
  expect_true(lg$log_scale)

  id <- out[out$yi_source == "ci_identity", ]
  expect_equal(id$yi, 1.0)
  expect_equal(id$sei, 2 / (2 * z), tolerance = 1e-12)
  expect_false(id$log_scale)
})

test_that("RoB2 attaches at estimate level with study-level fallback", {
  d <- withr::local_tempdir()
  p_zip <- file.path(d, "zip.csv"); make_zip_fixture(p_zip)
  p_rm5 <- file.path(d, "rm5.csv"); make_rm5_fixture(p_rm5)
  p_est <- file.path(d, "RoB2.csv"); p_stu <- file.path(d, "RoB2b.csv")
  make_rob2_fixtures(p_est, p_stu)

  zip <- suppressMessages(read_zip_estimates(p_zip))
  rm5 <- suppressMessages(read_rm5_estimates(p_rm5))
  rob <- read_rob2(p_est, p_stu)

  z <- attach_rob2(zip, rob)
  delta <- z[z$study == "Delta 2015" & z$analysis_number == "1", ]
  expect_equal(delta$rob2_bias_arising_from_the_randomization_process, "low")
  expect_equal(
    delta$rob2_bias_due_to_deviations_from_intended_interventions,
    "some_concerns")
  expect_equal(delta$rob2_bias_due_to_missing_outcome_data, "high")
  expect_true(is.na(delta$rob2_bias_in_selection_of_the_reported_result))

  expect_equal(delta$rob2_level, "estimate")

  # Epsilon has no estimate-level row: falls back to study level
  eps <- z[z$study == "Epsilon 2018", ]
  expect_equal(eps$rob2_overall_bias, "high")
  expect_equal(eps$rob2_level, "study")

  # Delta's *other* analysis (analysis_number 2) misses the estimate-level
  # key but the same (review, study) exists in the estimate table: filled
  # by the collapsed study-level stage
  delta2 <- z[z$study == "Delta 2015" & z$analysis_number == "2", ]
  expect_equal(delta2$rob2_overall_bias, "some_concerns")
  expect_equal(delta2$rob2_level, "study_from_estimates")

  # Eta's estimate-level row exists but carries empty judgements: the key
  # match wins the level, the values are backfilled from the collapsed
  # study stage (Eta's other analysis says "high")
  eta <- z[z$study == "Eta 2021", ]
  expect_equal(eta$rob2_level, "estimate")
  expect_equal(eta$rob2_overall_bias, "high")

  # rm5 rows only ever match at study level; review ids are canonicalized
  # on both sides (14651858.CD000001.pub2 == CD000001)
  r <- attach_rob2(rm5, rob)
  alpha <- r[r$study == "Alpha 2001", ]
  expect_true(all(alpha$rob2_overall_bias == "high"))
  expect_true(all(alpha$rob2_level == "study"))
  expect_true(all(is.na(r[r$study == "Beta 2003", ]$rob2_overall_bias)))
})

test_that("norm_review_id canonicalizes both id formats", {
  expect_equal(norm_review_id("14651858.CD000009.pub4"), "CD000009")
  expect_equal(norm_review_id("CD000031"), "CD000031")
  expect_equal(norm_review_id(" cd000031 "), "CD000031")
})

test_that("collapse_rob2_study takes the mode, ties broken toward severe", {
  est <- data.frame(
    Review = "CD9", Study = c("S1", "S1", "S1", "S2", "S2"),
    stringsAsFactors = FALSE
  )
  for (d in paste0("rob2_", gtmeta:::ROB2_DOMAINS)) {
    est[[d]] <- c("low", "low", "high", "low", "high")
  }
  coll <- gtmeta:::collapse_rob2_study(est)
  s1 <- coll[coll$key == gtmeta:::rob2_norm_key("CD9", "S1"), ]
  s2 <- coll[coll$key == gtmeta:::rob2_norm_key("CD9", "S2"), ]
  expect_equal(s1$rob2_overall_bias, "low")    # 2x low beats 1x high
  expect_equal(s2$rob2_overall_bias, "high")   # 1-1 tie -> more severe
})

test_that("cochrane_tidy orchestrates a full directory tree", {
  d <- withr::local_tempdir()
  dir.create(file.path(d, "out.03.rm5s"))
  dir.create(file.path(d, "out.03.zips"))
  make_rm5_fixture(file.path(d, "out.03.rm5s", "rm5_data.csv.gz"))
  make_zip_fixture(file.path(d, "out.03.zips", "zip_data.csv.gz"))
  make_rob2_fixtures(file.path(d, "out.03.zips", "RoB2.csv"),
                     file.path(d, "out.03.zips", "RoB2b.csv"))

  out <- suppressMessages(cochrane_tidy(d))
  expect_equal(nrow(out), 7L)
  expect_setequal(unique(out$source), c("rm5", "zip"))
  expect_true(all(c("yi", "sei", "outcome_id", "network_id",
                    "rob2_overall_bias") %in% names(out)))
  expect_true(all(is.finite(out$yi)))
})

test_that("igmi_study_list builds N(theta, Sigma) with plug-in correlation", {
  tidy <- data.frame(
    review = "CD1",
    study = rep(c("S1", "S2", "S3"), each = 2),
    outcome_id = rep(c("CD1:1:1", "CD1:2:1"), 3),
    yi = c(0.1, 0.5, 0.2, 0.6, NA, 0.7),
    sei = c(0.10, 0.20, 0.15, 0.25, 0.10, 0.30)
  )
  out <- igmi_study_list(tidy, "CD1", c("CD1:1:1", "CD1:2:1"), rho = 0.3)
  # S3 lacks the first outcome -> only S1, S2 usable
  expect_length(out, 2L)
  expect_equal(unname(out[[1]]$theta), c(0.1, 0.5))
  expect_equal(out[[1]]$Sigma[1, 1], 0.10^2, tolerance = 1e-12)
  expect_equal(out[[1]]$Sigma[1, 2], 0.3 * 0.10 * 0.20, tolerance = 1e-12)
  expect_true(isSymmetric(out[[1]]$Sigma))

  # rho = 0 gives a diagonal Sigma
  out0 <- igmi_study_list(tidy, "CD1", c("CD1:1:1", "CD1:2:1"), rho = 0)
  expect_equal(out0[[1]]$Sigma[1, 2], 0)
})

test_that("igmi_study_list pools subgroup duplicates by fixed effect", {
  # S1 has two subgroup rows on outcome 1: FE pool = classical
  # inverse-variance combination (hand-checked)
  tidy <- data.frame(
    review = "CD1",
    study = c("S1", "S1", "S1", "S2", "S2"),
    outcome_id = c("CD1:1:1", "CD1:1:1", "CD1:2:1", "CD1:1:1", "CD1:2:1"),
    yi = c(0.2, 0.6, 0.5, 0.1, 0.4),
    sei = c(0.1, 0.2, 0.2, 0.1, 0.3)
  )
  expect_message(
    out <- igmi_study_list(tidy, "CD1", c("CD1:1:1", "CD1:2:1"), rho = 0),
    "1 subgroup rows pooled")
  expect_length(out, 2L)
  w <- c(1 / 0.1^2, 1 / 0.2^2)
  expect_equal(unname(out[[1]]$theta[1]), sum(w * c(0.2, 0.6)) / sum(w),
               tolerance = 1e-12)
  expect_equal(out[[1]]$Sigma[1, 1], 1 / sum(w), tolerance = 1e-12)
  # the non-duplicated outcome is untouched
  expect_equal(unname(out[[1]]$theta[2]), 0.5)

  # "first" keeps the first row; "error" refuses
  out_first <- igmi_study_list(tidy, "CD1", c("CD1:1:1", "CD1:2:1"),
                               duplicates = "first")
  expect_equal(unname(out_first[[1]]$theta[1]), 0.2)
  expect_error(
    igmi_study_list(tidy, "CD1", c("CD1:1:1", "CD1:2:1"),
                    duplicates = "error"),
    "duplicated")
})

test_that("attach_rob2 warns instead of silently skipping the estimate-level join (M4)", {
  d <- withr::local_tempdir()
  p_zip <- file.path(d, "zip.csv"); make_zip_fixture(p_zip)
  p_est <- file.path(d, "RoB2.csv"); p_stu <- file.path(d, "RoB2b.csv")
  make_rob2_fixtures(p_est, p_stu)

  zip <- suppressMessages(read_zip_estimates(p_zip))
  rob <- read_rob2(p_est, p_stu)

  # simulate a read.csv dot-conversion that renames Analysis_group ->
  # Analysis.group: the estimate-level key column is now missing
  broken <- rob
  names(broken$estimate)[names(broken$estimate) == "Analysis_group"] <-
    "Analysis.group"

  expect_warning(z <- attach_rob2(zip, broken),
                 "estimate-level RoB2 join skipped")
  # study-level fallback still populates rows that have a study match
  expect_true(any(z$rob2_level %in% c("study", "study_from_estimates")))
  # and no row is wrongly labelled as an estimate-level match
  expect_false(any(z$rob2_level %in% "estimate", na.rm = TRUE))
})
