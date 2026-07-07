# The bridge tests are environmental: they run only where the pinned
# virtualenv exists (wfr_setup() has been run). CI and fresh checkouts
# skip them; the WFR *oracle* tests (two-atom closed form, monotone
# bracket) are Phase 3, in test_wfr_limit.R.

test_that("pinned dependency table is well-formed", {
  deps <- wfr_py_deps()
  expect_named(deps, c("numpy", "scipy", "POT"))
  expect_true(all(grepl("^[0-9]+\\.[0-9]+", deps)))
})

test_that("bridge smoke test passes where the env is present", {
  skip_if_not_installed("reticulate")
  skip_if(!wfr_available(), "pinned WFR virtualenv not set up")
  expect_true(wfr_smoke_test())
})
