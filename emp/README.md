# emp/ — Phase-5 empirical validation (Cochrane)

Empirical drivers over the real cochrane2025rob corpus. Not part of the
built package (`.Rbuildignore`); `emp/results/` is gitignored — outputs
are reproducible from these scripts plus the corpus.

## Requirements

A cochrane2025rob output tree (the four LFS files, `out.03.*` layout);
see `R/data_adapter.R`. Point `GTMETA_COCHRANE_DIR` at it (default:
`../cochrane2025rob_data` next to the repo).

## Scripts

- `emp_univariate.R` — Phase 5.B univariate arm: IGMI beside FE/RE/UWLS
  across one meta-analysis per review (largest outcome, K ≥ 5,
  subgroup rows FE-pooled), mirroring the AIC/BIC design of Stanley et
  al. 2023 plus LOO predictive scoring. IGMI-prop ≡ UWLS exactly
  (igmi_notes Prop 5.16) and is identified, not duplicated; IGMI-WFR =
  equal-weight `wfr_pool()` at `delta = mad(yi)` with the sandwich
  standard error (igmi_notes Rem 4.9 updates), LOO-scored with
  predictive `N(x_-i, v_i + tau2_DL,-i + se_-i^2)`.

  ```sh
  RENV_CONFIG_AUTOLOADER_ENABLED=false Rscript emp/emp_univariate.R
  ```

  Env overrides: `GTMETA_EMP_KMIN` (default 5), `GTMETA_EMP_MAX`
  (cap the number of meta-analyses; 0 = all), `GTMETA_COCHRANE_DIR`.
  Full corpus ≈ 100 s. Results: `emp/results/emp_univariate.rds`
  (+ `.log`); the tidy corpus is cached in
  `emp/results/tidy_cache.rds` (delete to re-run the adapter).

- `emp_multivariate.R` — Phase 5.B multivariate arm: bivariate MAs from
  multi-outcome reviews, one per review, outcome pairs from *different*
  analysis groups (same-group outcome_ids are sensitivity variants),
  K ≥ 5 complete cases. Contenders under joint LOO predictive scoring:
  independent per-outcome FE/DL, `mvmeta` (REML, ML fallback), and
  IGMI-BW (linear barycenter mean, trace-precision weights, matrix
  moment estimator `igmi_mm_cov()` — igmi_notes Rem 5.12 Phase-5.B
  update). Within-study cross-outcome correlation is a plug-in
  sensitivity parameter, `rho ∈ {0, 0.3, 0.6}`; the IGMI-BW point
  estimate is rho-invariant by construction. Also compares the exact
  Hotelling region (Prop 5.16) with the mvmeta Wald ellipse and
  verifies BW fixed-point convergence (Thm 7.5) corpus-wide.

  ```sh
  RENV_CONFIG_AUTOLOADER_ENABLED=false Rscript emp/emp_multivariate.R
  ```

  Same env overrides. Results: `emp/results/emp_multivariate.rds`
  (+ `.log`).
