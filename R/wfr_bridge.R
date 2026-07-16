# wfr_bridge.R -- pinned reticulate -> Python (POT) bridge for WFR
#
# IGMI's Wasserstein-Fisher-Rao barycenter has no mature native R
# implementation; the plan (plans/README.md) is to call Python's POT
# through reticulate. This file only stands up and checks the bridge
# (Phase 2); the WFR computations themselves land in igmi_wfr.R
# (Phase 3), following the pinned POT recipe of igmi_notes Rem 4.9
# (cost `l_delta`, `reg_m = 1`, value `* 4 delta^2`).
#
# Versions are pinned so that simulation and empirical results are
# reproducible; bump them deliberately, in one place, here.

#' Pinned Python dependencies for the WFR bridge
#' @return Named character vector `package = version`.
#' @export
wfr_py_deps <- function() {
  c(numpy = "2.2.6", scipy = "1.18.0", POT = "0.9.5")
}

.wfr_envname <- "r-gtmeta"

#' Install the pinned Python dependencies for the WFR bridge
#'
#' This function installs software on the user's system: it creates a
#' reticulate virtualenv (named `r-gtmeta` by default) and installs the
#' exact Python package versions of [wfr_py_deps()] into it. Requires a
#' system `python3`.
#'
#' Nothing is installed unless the user explicitly consents: in an
#' interactive session the function asks for confirmation first, and in
#' a non-interactive session it stops unless called with
#' `confirm = TRUE`. It is never called by the package itself, its
#' examples, tests, or vignette.
#'
#' @param envname Virtualenv name.
#' @param confirm Set to `TRUE` to consent to the installation without
#'   the interactive prompt (required in non-interactive sessions).
#' @return The environment name, invisibly. Called for its side effect
#'   of creating the virtualenv and installing the pinned Python
#'   packages into it.
#' @export
wfr_setup <- function(envname = .wfr_envname, confirm = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("wfr_setup() needs the 'reticulate' package")
  }
  deps <- wfr_py_deps()
  pkgs <- paste0(names(deps), "==", deps)
  if (!isTRUE(confirm)) {
    if (!interactive()) {
      stop("wfr_setup() installs Python packages (",
           paste(pkgs, collapse = ", "),
           ") into the '", envname, "' virtualenv. ",
           "Call wfr_setup(confirm = TRUE) to consent.")
    }
    ok <- utils::askYesNo(paste0(
      "Create virtualenv '", envname, "' and install ",
      paste(pkgs, collapse = ", "), "?"))
    if (!isTRUE(ok)) {
      message("Installation cancelled.")
      return(invisible(envname))
    }
  }
  if (!reticulate::virtualenv_exists(envname)) {
    reticulate::virtualenv_create(envname)
  }
  reticulate::virtualenv_install(envname, packages = pkgs)
  invisible(envname)
}

#' Is the pinned WFR bridge usable?
#'
#' @param envname Virtualenv name.
#' @return `TRUE` if the virtualenv exists and `ot` imports at the pinned
#'   version, else `FALSE` (with a message).
#' @export
wfr_available <- function(envname = .wfr_envname) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  if (!reticulate::virtualenv_exists(envname)) return(FALSE)
  ok <- tryCatch({
    reticulate::use_virtualenv(envname, required = TRUE)
    ot <- reticulate::import("ot")
    ver <- as.character(ot$`__version__`)
    if (!identical(ver, unname(wfr_py_deps()["POT"]))) {
      message("POT version ", ver, " differs from pinned ",
              wfr_py_deps()["POT"])
    }
    TRUE
  }, error = function(e) {
    message("WFR bridge not usable: ", conditionMessage(e))
    FALSE
  })
  ok
}

#' Smoke test of the WFR bridge
#'
#' Runs a balanced 2 x 2 exact OT problem through `ot.emd2` and checks
#' the value against the obvious by-hand optimum. This verifies the
#' end-to-end R -> reticulate -> numpy -> POT round trip, not WFR itself
#' (the WFR oracle tests, e.g. the two-atom closed form of igmi_notes
#' Prop 4.8 and the monotone bracket of Cor 4.13, belong to Phase 3's
#' test_wfr_limit.R).
#'
#' @param envname Virtualenv name.
#' @return `TRUE` on success; errors otherwise.
#' @export
wfr_smoke_test <- function(envname = .wfr_envname) {
  stopifnot(wfr_available(envname))
  ot <- reticulate::import("ot")
  np <- reticulate::import("numpy")
  a <- np$array(c(0.5, 0.5))
  b <- np$array(c(0.5, 0.5))
  M <- np$array(matrix(c(0, 1, 1, 0), 2, 2))
  val <- ot$emd2(a, b, M)
  stopifnot(abs(val - 0) < 1e-12)
  M2 <- np$array(matrix(c(1, 0, 0, 1), 2, 2))
  val2 <- ot$emd2(a, b, M2)
  stopifnot(abs(val2 - 0) < 1e-12)
  M3 <- np$array(matrix(c(2, 1, 1, 2), 2, 2))
  val3 <- ot$emd2(a, b, M3)
  stopifnot(abs(val3 - 1) < 1e-12)
  TRUE
}
