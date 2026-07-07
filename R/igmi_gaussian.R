# igmi_gaussian.R -- study objects N(theta_i, Sigma_i) for IGMI
#
# Constructors and validation for the study objects of igmi_notes
# Definition 1.1: each study is the point N(theta_i, Sigma_i) on the
# manifold of nondegenerate m-variate Gaussians, with weights a modelling
# decision separate from the geometry.

.check_spd <- function(S, name = "Sigma") {
  if (!is.matrix(S) || nrow(S) != ncol(S)) {
    stop(name, " must be a square matrix")
  }
  if (max(abs(S - t(S))) > 1e-8 * max(1, max(abs(S)))) {
    stop(name, " must be symmetric")
  }
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) {
    stop(name, " must be positive definite (min eigenvalue ",
         format(min(ev)), ")")
  }
  invisible((S + t(S)) / 2)
}

#' Construct one IGMI study object
#'
#' The study object `N(theta, Sigma)` of igmi_notes Definition 1.1.
#' Scalars are accepted for `m = 1`: `igmi_gaussian(0.3, sei = 0.1)`.
#'
#' @param theta Effect estimate (length-m numeric).
#' @param Sigma Sampling covariance (m x m positive definite). For
#'   `m = 1` you may instead pass `sei`.
#' @param sei Standard error (only for `m = 1`, instead of `Sigma`).
#' @return Object of class `igmi_gaussian`: list `(theta, Sigma, m)`.
#' @export
igmi_gaussian <- function(theta, Sigma = NULL, sei = NULL) {
  theta <- as.numeric(theta)
  m <- length(theta)
  if (is.null(Sigma)) {
    if (is.null(sei) || m != 1L) {
      stop("provide Sigma (or, for m = 1, sei)")
    }
    stopifnot(is.numeric(sei), length(sei) == 1L, sei > 0)
    Sigma <- matrix(sei^2, 1, 1)
  }
  if (!is.matrix(Sigma)) Sigma <- as.matrix(Sigma)
  stopifnot(nrow(Sigma) == m)
  Sigma <- .check_spd(Sigma)
  structure(list(theta = theta, Sigma = Sigma, m = m),
            class = "igmi_gaussian")
}

#' Construct a list of IGMI study objects
#'
#' Univariate: pass `yi` and `sei`. Multivariate: pass `theta` (K x m
#' matrix) and `Sigma` (list of K covariance matrices) -- the output
#' format of [sim_multivariate()] and [igmi_study_list()].
#'
#' @param yi,sei Univariate effect estimates and standard errors.
#' @param theta K x m matrix of multivariate estimates.
#' @param Sigma List of K m x m covariances.
#' @return List of `igmi_gaussian` objects.
#' @export
igmi_studies <- function(yi = NULL, sei = NULL, theta = NULL, Sigma = NULL) {
  if (!is.null(yi)) {
    stopifnot(!is.null(sei), length(yi) == length(sei))
    return(lapply(seq_along(yi), function(i) {
      igmi_gaussian(yi[i], sei = sei[i])
    }))
  }
  stopifnot(is.matrix(theta), is.list(Sigma), nrow(theta) == length(Sigma))
  lapply(seq_len(nrow(theta)), function(i) {
    igmi_gaussian(theta[i, ], Sigma[[i]])
  })
}

.check_studies <- function(studies) {
  stopifnot(length(studies) >= 1L,
            all(vapply(studies, inherits, logical(1), "igmi_gaussian")))
  m <- studies[[1]]$m
  stopifnot(all(vapply(studies, `[[`, integer(1), "m") == m))
  m
}

.check_weights <- function(weights, K) {
  if (is.null(weights)) weights <- rep(1 / K, K)
  stopifnot(length(weights) == K, all(weights > 0))
  weights / sum(weights)
}

#' Precision weights for a list of studies
#'
#' `w_i` proportional to `1 / tr(Sigma_i)`, normalised to sum to one.
#' For `m = 1` this is the classical inverse-variance weight, so the
#' Bures-Wasserstein barycenter mean under these weights is the
#' fixed-effect estimate (igmi_notes Cor 3.6). The weight choice is a
#' modelling decision (igmi_notes Def 1.1), not part of the geometry.
#'
#' @param studies List of `igmi_gaussian` objects.
#' @return Normalised weight vector.
#' @export
igmi_precision_weights <- function(studies) {
  .check_studies(studies)
  w <- vapply(studies, function(s) 1 / sum(diag(s$Sigma)), numeric(1))
  w / sum(w)
}
