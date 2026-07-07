# gtmeta 0.1.0

First release.

* Study constructors for univariate (`igmi_studies()`), multivariate
  (`igmi_gaussian()`), diagnostic test accuracy (`igmi_dta_studies()`), and
  prediction-model performance study objects.
* Bures-Wasserstein barycenter pooling (`bw_barycenter()`) with the exact
  fixed-point iteration, descent/trace invariants as runtime checks, and
  reduction to the inverse-variance fixed-effect estimate in the scalar
  precision-weighted case.
* Fisher-Rao distances and means with the Karcher-ball admissibility check.
* Frechet-variance heterogeneity measures that recover classical `Q`, `I2`,
  and the DerSimonian-Laird moment estimator exactly in the scalar limit.
* Exact location inference (chi-square ellipsoids), Frechet-scatter
  (HKSJ-type) inference, and bootstrap-over-studies defaults.
* Robust Wasserstein-Fisher-Rao pooling (`wfr_pool()`), a redescending
  Andrews-sine M-estimator with transport cutoff `pi * delta`, plus an
  optional Python/POT bridge for general WFR distances.
* Classical fixed-effect, random-effects, and UWLS benchmarks with AIC/BIC
  and leave-one-study-out predictive scoring.
* Simulation generators, evaluation metrics, and adapters for Cochrane-style
  extracted estimates.
* Theorem-anchored test suite: every estimator is tested against a numbered
  result in the accompanying mathematical notes.
