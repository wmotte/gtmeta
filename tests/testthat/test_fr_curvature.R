# Promotes (part of) the scratch verification of igmi_notes
# Rem 6.22(c): the invariant sectional-curvature formula eq. (6b) at
# N(0, I) and its anchors from Thm 6.19,
#   K = (1/16) [ 8 Gamma - 2 ||[A1,A2] + u1 u2' - u2 u1'||_F^2
#                - 4 ||A2 u1 - A1 u2||^2 ]      (g_FR-orthonormal pairs)
# with Gamma = |u1|^2 |u2|^2 - (u1.u2)^2. Anchors: every m = 1 plane has
# K = -1/2; mean-mean planes have K = +1/4; the explicit maximizer of
# Thm 6.19(3) attains exactly kappa_m = 2/7; and no random plane exceeds
# 2/7. (The finite-difference comparison against the Fisher metric
# remains a scratch-level check -- it needs numerical Christoffel
# machinery that does not belong in the package.)

K_fr <- function(u1, A1, u2, A2) {
  ip <- function(ua, Aa, ub, Ab) sum(ua * ub) + 0.5 * sum(Aa * Ab)
  g11 <- ip(u1, A1, u1, A1)
  g22 <- ip(u2, A2, u2, A2)
  g12 <- ip(u1, A1, u2, A2)
  gam <- sum(u1^2) * sum(u2^2) - sum(u1 * u2)^2
  K1 <- A1 %*% A2 - A2 %*% A1 + outer(u1, u2) - outer(u2, u1)
  uterm <- A2 %*% u1 - A1 %*% u2
  (8 * gam - 2 * sum(K1^2) - 4 * sum(uterm^2)) / 16 /
    (g11 * g22 - g12^2)
}

test_that("m = 1 planes have constant curvature -1/2 (Prop 2.2(1))", {
  expect_equal(K_fr(1, matrix(0), 0, matrix(sqrt(2))), -0.5,
               tolerance = 1e-12)
  # any basis of the same plane gives the same value
  expect_equal(K_fr(2, matrix(1), 1, matrix(3)), -0.5, tolerance = 1e-12)
})

test_that("mean-mean planes have curvature +1/4 (Thm 6.19(2))", {
  Z <- matrix(0, 2, 2)
  expect_equal(K_fr(c(1, 0), Z, c(0, 1), Z), 0.25, tolerance = 1e-12)
  Z3 <- matrix(0, 3, 3)
  expect_equal(K_fr(c(1, 1, 0) / sqrt(2), Z3, c(0, 0, 2), Z3), 0.25,
               tolerance = 1e-12)
})

test_that("the explicit maximizer attains kappa_m = 2/7 (Thm 6.19(3))", {
  a <- 1 / sqrt(7); r <- sqrt(6 / 7)
  u1 <- r * c(1, 0); A1 <- a * diag(c(1, -1))
  u2 <- r * c(0, 1); A2 <- -a * matrix(c(0, 1, 1, 0), 2)
  expect_equal(K_fr(u1, A1, u2, A2), 2 / 7, tolerance = 1e-14)
  # the same block embedded in m = 3 still attains 2/7 (dimension-free)
  pad <- function(u, A) {
    list(u = c(u, 0), A = rbind(cbind(A, 0), 0))
  }
  p1 <- pad(u1, A1); p2 <- pad(u2, A2)
  expect_equal(K_fr(p1$u, p1$A, p2$u, p2$A), 2 / 7, tolerance = 1e-14)
})

test_that("random planes never exceed 2/7 (Thm 6.19(3), scan)", {
  set.seed(51)
  for (m in 2:3) {
    worst <- -Inf
    for (t in seq_len(4000)) {
      u1 <- rnorm(m); u2 <- rnorm(m)
      A1 <- matrix(rnorm(m * m), m); A1 <- (A1 + t(A1)) / 2
      A2 <- matrix(rnorm(m * m), m); A2 <- (A2 + t(A2)) / 2
      worst <- max(worst, K_fr(u1, A1, u2, A2))
    }
    expect_lte(worst, 2 / 7 + 1e-12)
    expect_gt(worst, 0.2)   # the scan does approach the supremum
  }
})
