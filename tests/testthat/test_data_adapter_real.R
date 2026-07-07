# Opt-in real-data smoke test against an actual cochrane2025rob output
# tree (the LFS files of github.com/wmotte/cochrane2025rob). Point
# GTMETA_COCHRANE_DIR at a directory containing out.03.rm5s/ and
# out.03.zips/; skipped otherwise (CI has no data checkout).
#
# The assertions pin the corpus facts established 2026-07-03:
# disjoint rm5/zip corpora, RoB2 coverage limited to zip reviews, and a
# 100% RoB2 match rate within covered reviews after the three-stage join.

test_that("cochrane_tidy digests the real corpus", {
  dir <- Sys.getenv("GTMETA_COCHRANE_DIR", "")
  skip_if(dir == "" || !dir.exists(dir),
          "GTMETA_COCHRANE_DIR not set to a cochrane2025rob checkout")

  tidy <- suppressMessages(suppressWarnings(cochrane_tidy(dir)))

  expect_gt(nrow(tidy), 5e5)
  expect_true(all(is.finite(tidy$yi)))
  expect_true(all(tidy$sei > 0))
  expect_setequal(unique(tidy$source), c("rm5", "zip"))

  # canonical ids: every cdid is a bare CDxxxxxx
  expect_true(all(grepl("^CD[0-9]+$", tidy$cdid)))

  # the two corpora are disjoint at review level (row-bind = no dupes)
  expect_length(intersect(tidy$cdid[tidy$source == "rm5"],
                          tidy$cdid[tidy$source == "zip"]), 0L)

  # RoB2: within covered reviews every row matches at some level
  rob <- read_rob2(file.path(dir, "out.03.zips", "RoB2.csv"))
  covered <- unique(norm_review_id(rob$estimate$Review))
  zt <- tidy[tidy$source == "zip" & tidy$cdid %in% covered, ]
  expect_gt(nrow(zt), 1000L)
  # every covered row key-matches at estimate level; empty judgements are
  # backfilled from the collapsed study stage to > 90% coverage (the rest
  # are studies whose every RoB2 row is empty)
  expect_true(all(zt$rob2_level == "estimate"))
  expect_gt(mean(!is.na(zt$rob2_overall_bias)), 0.9)
  expect_true(all(zt$rob2_overall_bias %in%
                    c("low", "some_concerns", "high", NA_character_)))
  expect_true(all(tidy$rob2_level %in%
                    c("estimate", "study", "study_from_estimates",
                      NA_character_)))

  # rm5 rows carry no RoB2 (coverage, not key mismatch)
  expect_true(all(is.na(tidy$rob2_level[tidy$source == "rm5"])))
})
