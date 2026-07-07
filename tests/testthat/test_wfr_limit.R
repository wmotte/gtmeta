# The named Phase-3 WFR tests (TODO 3.B): oracles from igmi_notes
# Prop 4.8 (two-atom closed form; Hellinger anchor; cutoff pi*delta)
# and Cor 4.13 (monotone bracket WFR_{d1} <= WFR_{d2} <= W_2 with the
# closed-form Gaussian ceiling). Environmental: they run only where the
# pinned POT virtualenv exists (wfr_setup()).

skip_if_no_wfr <- function() {
  skip_if_not_installed("reticulate")
  skip_if(!wfr_available(), "pinned WFR virtualenv not set up")
}

test_that("two-atom closed form (Prop 4.8) is reproduced by POT", {
  skip_if_no_wfr()
  delta <- 0.7
  cases <- list(
    list(a0 = 1.0, a1 = 1.0, d = 0.5),   # transport regime
    list(a0 = 0.6, a1 = 1.3, d = 1.0),   # unbalanced, transport
    list(a0 = 1.0, a1 = 0.8, d = 0.0),   # co-located: Hellinger anchor
    list(a0 = 0.9, a1 = 1.1, d = 10.0)   # beyond cutoff pi*delta
  )
  for (cs in cases) {
    got <- wfr_dist(cs$a0, 0, cs$a1, cs$d, delta = delta)
    want <- wfr_two_atom_oracle(cs$a0, cs$a1, cs$d, delta)
    expect_equal(got$wfr2, want, tolerance = 1e-6)
  }
  # unit-atom birth cost is exactly (2 delta)^2 (Prop 4.8(1))
  got0 <- wfr_dist(1, 0, 1e-12, 0, delta = delta)
  expect_equal(got0$wfr2, 4 * delta^2, tolerance = 1e-5)
})

test_that("beyond the cutoff the value is pure kill-and-create", {
  skip_if_no_wfr()
  delta <- 0.5
  got <- wfr_dist(c(0.5, 0.5), c(-10, 10), 1, 0, delta = delta)
  expect_equal(got$wfr2, 4 * delta^2 * (1 + 1), tolerance = 1e-6)
  # and the plan carries no mass (all pairs are capped)
  expect_lt(max(got$plan), 1e-10)
})

test_that("monotone bracket WFR_delta <= W2 for Gaussian pairs (Cor 4.13)", {
  skip_if_no_wfr()
  g0 <- wfr_gauss_grid(0, 1, n = 60, span = 5)
  g1 <- wfr_gauss_grid(1.2, 1.5, n = 60, span = 5)

  w2 <- w2_dist2_discrete(g0$a, g0$x, g1$a, g1$x)
  # discrete W2 close to the continuous closed form (Thm 2.4, m = 1):
  # ||mu0-mu1||^2 + (s0-s1)^2
  expect_equal(w2, 1.2^2 + 0.5^2, tolerance = 0.02)

  deltas <- c(0.5, 1, 2, 4, 12)
  vals <- vapply(deltas, function(d) {
    wfr_dist(g0$a, g0$x, g1$a, g1$x, delta = d,
             max_iter = 100000L)$wfr2
  }, numeric(1))
  # nondecreasing in delta (Lem 4.5(4) pointwise / Thm 4.12(2))
  expect_true(all(diff(vals) >= -1e-8))
  # bracketed by the balanced W2 on the same discretization; the MM
  # solver certifies an upper bound of the true WFR^2, so allow a small
  # solver slack on top of the exact inequality WFR <= W2
  expect_true(all(vals <= w2 * (1 + 2e-3)))
  # and by delta = 12 the gap has essentially closed (the O(delta^-2)
  # tail correction scales with the fourth moment of the displacement,
  # so moderate deltas still show a visible gap)
  expect_equal(vals[length(vals)], w2, tolerance = 0.02)
})

test_that("equal-mass atom pair: gap to W2 closes at rate O(delta^-2)", {
  skip_if_no_wfr()
  # Prop 4.8(3): for a0 = a1 = a, WFR^2 = 8 delta^2 a (1 - cos(eps/(2 delta)))
  # = a eps^2 (1 + O(eps^2/delta^2))
  a <- 1; eps <- 0.8
  w2 <- a * eps^2
  gap <- vapply(c(2, 4, 8), function(d) {
    w2 - wfr_two_atom_oracle(a, a, eps, d)
  }, numeric(1))
  # halving 1/delta^2 four-folds the gap closure
  expect_equal(gap[1] / gap[2], 4, tolerance = 0.05)
  expect_equal(gap[2] / gap[3], 4, tolerance = 0.05)
  # and POT agrees with the oracle at one of these deltas
  got <- wfr_dist(a, 0, a, eps, delta = 2)
  expect_equal(got$wfr2, wfr_two_atom_oracle(a, a, eps, 2),
               tolerance = 1e-6)
})
