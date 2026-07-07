test_that("scalar metrics compute known values", {
  est <- c(1.1, 0.9, 1.0, 1.2)
  expect_equal(metric_bias(est, 1), 0.05)
  expect_equal(metric_mse(est, 1), mean(c(0.01, 0.01, 0, 0.04)))
  expect_equal(metric_coverage(c(0, 0, 2), c(2, 0.5, 3), 1), 1 / 3)
  expect_equal(metric_ci_width(c(0, 1), c(2, 4)), 2.5)
})

test_that("metric_calibration is happy with uniform PIT", {
  pit <- (1:999) / 1000
  cal <- metric_calibration(pit)
  expect_lt(cal$ks, 0.01)
  expect_gt(cal$ks_p, 0.9)
  expect_equal(cal$z_mean, 0, tolerance = 1e-10)
  expect_equal(cal$z_sd, 1, tolerance = 0.05)

  # grossly overdispersed PIT is flagged (coarser grid keeps the
  # transformed values strictly inside (0, 1) in double precision)
  bad <- stats::pnorm(stats::qnorm((1:199) / 200) * 3)
  cal_bad <- metric_calibration(bad)
  expect_gt(cal_bad$z_sd, 2.5)
  expect_lt(cal_bad$ks_p, 0.05)
})

test_that("detection_power is the rejection rate", {
  expect_equal(detection_power(c(0.01, 0.2, 0.04, 0.6)), 0.5)
  expect_equal(detection_power(c(0.01, 0.2), alpha = 0.001), 0)
})

test_that("localization_accuracy ranks the true candidate", {
  scores <- rbind(c(0.1, 0.9, 0.3),   # picks 2 (true)
                  c(0.8, 0.2, 0.1),   # picks 1 (true = 2 -> rank 2)
                  c(0.2, 0.5, 0.9))   # picks 3 (true = 2 -> rank 2)
  loc <- localization_accuracy(scores, true_index = 2)
  expect_equal(loc$top1, 1 / 3)
  expect_equal(loc$mean_rank, mean(c(1, 2, 2)))

  # vector input = single replication
  loc1 <- localization_accuracy(c(3, 1, 2), true_index = 1)
  expect_equal(loc1$top1, 1)
  expect_equal(loc1$mean_rank, 1)
})

test_that("metric_table aggregates over replications", {
  set.seed(5)
  n <- 5000L
  est <- stats::rnorm(n, 1, 0.1)
  se <- rep(0.1, n)
  tab <- metric_table(est, se, est - 1.96 * se, est + 1.96 * se, truth = 1)
  expect_equal(tab$bias, 0, tolerance = 0.01)
  expect_equal(tab$se_ratio, 1, tolerance = 0.03)
  expect_equal(tab$coverage, 0.95, tolerance = 0.01)
  expect_equal(tab$rmse, 0.1, tolerance = 0.02)
})
