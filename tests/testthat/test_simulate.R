test_that("sim_meta reproduces its variance decomposition", {
  set.seed(1)
  k <- 20000L
  d <- sim_meta(k, mu = 0.5, tau2 = 0.09, sei = rep(0.2, k))
  expect_equal(nrow(d), k)
  expect_equal(mean(d$yi), 0.5, tolerance = 0.01)
  expect_equal(stats::var(d$yi), 0.09 + 0.04, tolerance = 0.005)
  expect_equal(attr(d, "tau2"), 0.09)

  # tau2 = 0: theta degenerate at mu
  d0 <- sim_meta(10, mu = 1, tau2 = 0)
  expect_true(all(d0$theta == 1))
})

test_that("sim_multivariate builds the requested Sigma structure", {
  set.seed(2)
  mu <- c(0.2, -0.1, 0.4)
  Tau <- diag(c(0.04, 0.01, 0.09))
  s <- sim_multivariate(5000, mu, Tau, rho_within = 0.5,
                        sei_range = c(0.2, 0.2))
  expect_equal(dim(s$y), c(5000L, 3L))
  # every within-study covariance has the plug-in equicorrelation
  S1 <- s$Sigma[[1]]
  expect_equal(S1[1, 2] / sqrt(S1[1, 1] * S1[2, 2]), 0.5,
               tolerance = 1e-12)
  # marginal covariance of y = Tau + E[Sigma_i] (absolute tolerance:
  # off-diagonal targets are small)
  ESig <- Reduce(`+`, s$Sigma) / length(s$Sigma)
  expect_lt(max(abs(stats::cov(s$y) - (Tau + ESig))), 0.01)
  expect_equal(colMeans(s$y), mu, tolerance = 0.02)
})

test_that("sim_contaminated flags and perturbs the right fraction", {
  set.seed(3)
  k <- 1000L
  d <- sim_contaminated(k, mu = 0, tau2 = 0, sei = rep(0.1, k),
                        eps = 0.2, shift = 5, understate = 2)
  expect_equal(sum(d$is_outlier), 200L)
  # outliers are shifted by ~5
  expect_equal(mean(d$yi[d$is_outlier]), 5, tolerance = 0.1)
  # and report half their true sd
  expect_true(all(d$sei[d$is_outlier] == 0.05))
  expect_true(all(d$sei[!d$is_outlier] == 0.1))
})

test_that("sim_contaminated contaminates at least one study when eps is positive", {
  set.seed(31)
  d <- sim_contaminated(k = 5L, mu = 0, tau2 = 0, sei = rep(0.1, 5),
                        eps = 0.1, shift = 5, understate = 2)
  expect_equal(sum(d$is_outlier), 1L)
  expect_true(all(d$sei[d$is_outlier] == 0.05))
})

test_that("sim_network is cycle-consistent iff omega vanishes", {
  edges <- rbind(c("A", "B"), c("B", "C"), c("A", "C"))
  d <- c(A = 0, B = 0.5, C = 1.2)

  set.seed(4)
  net0 <- sim_network(edges, d, omega = 0, n_per_edge = 1,
                      sei_range = c(1e-9, 1e-9))
  # signed cycle sum AB + BC - AC of the (noise-free) contrasts is zero
  cyc <- net0$yi[1] + net0$yi[2] - net0$yi[3]
  expect_equal(cyc, 0, tolerance = 1e-6)
  expect_equal(net0$yi[1], 0.5, tolerance = 1e-6)

  # inject omega = 0.3 on the BC edge: cycle sum becomes 0.3
  net1 <- sim_network(edges, d, omega = c(0, 0.3, 0), n_per_edge = 1,
                      sei_range = c(1e-9, 1e-9))
  cyc1 <- net1$yi[1] + net1$yi[2] - net1$yi[3]
  expect_equal(cyc1, 0.3, tolerance = 1e-6)
  expect_equal(unname(attr(net1, "omega")["B:C"]), 0.3)

  # recycling of n_per_edge
  net2 <- sim_network(edges, d, n_per_edge = c(2, 3, 4))
  expect_equal(nrow(net2), 9L)
})

test_that("sim_multivariate handles a singular (rank-deficient) PSD Tau (M3)", {
  # perfectly correlated between-study effects: Tau is PSD, rank 1, so
  # chol() would fail; the eigen square root must still draw effects that
  # lie (up to within-study noise) on the degenerate direction.
  set.seed(11)
  Tau <- matrix(c(0.2, 0.2, 0.2, 0.2), 2)   # rank 1
  d <- sim_multivariate(k = 5000L, mu = c(0, 0), Tau = Tau, rho_within = 0,
                        sei_range = c(1e-4, 1e-4))
  expect_equal(dim(d$theta), c(5000L, 2L))
  # true effects collinear: theta[,1] == theta[,2] before within-study noise
  expect_lt(max(abs(d$theta[, 1] - d$theta[, 2])), 1e-8)
  expect_equal(stats::var(d$theta[, 1]), 0.2, tolerance = 0.02)
})
