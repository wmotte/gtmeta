# sim/ — Phase 4.B simulation studies (IGMI)

Drivers for the TODO.md Phase-4.B evaluation. They source `R/` directly
(no installed `gtmeta` needed) and write results + logs to
`sim/results/`. Run from the repo root:

```sh
RENV_CONFIG_AUTOLOADER_ENABLED=false Rscript sim/sim_multivariate.R
RENV_CONFIG_AUTOLOADER_ENABLED=false Rscript sim/sim_contamination.R
```

Replication counts are env-overridable (`GTMETA_NREP`,
`GTMETA_NREP_VF`, `GTMETA_NREP_MAIN`, `GTMETA_NREP_BRK`,
`GTMETA_NREP_COV`, `GTMETA_B_BOOT`); every cell seeds its own RNG, so
results are reproducible cell by cell.

## sim_multivariate.R — arms A and C

* **Arm A** (m = 2): bias / MSE / coverage of the IGMI barycenter mean
  (precision weights) with the Fréchet-scatter pivot intervals
  (`igmi_pivot_ci`, Prop 5.16) and the known-Σ normal theory at T = 0
  (`igmi_location_ci`, Prop 5.13), against `mvmeta` REML Wald
  intervals. Grid: K ∈ {5, 10, 25} × between-study T ∈ {0, 0.05·I,
  0.05·equicorr(0.5)} × within-study correlation ∈ {0, 0.5}.
* **Arm C**: calibration of the V_F machinery — type-I error and power
  of `vf_null_test` (Prop 5.11 mixture), unbiasedness of the raw trace
  moment estimator `tr_T_hat_raw`.

## sim_contamination.R — arm B

Robustness/breakdown (m = 1) of FE / RE(REML, DL fallback) / UWLS
versus the WFR one-atom Fréchet mean `wfr_pool()` (igmi_notes Rem 4.9
Update; rejection point = transport cutoff πδ, Prop 4.8(2)). Since the
scalar balanced BW barycenter mean equals FE (Cor 3.6), the WFR gain
over FE under contamination is exactly the decision-log question "WFR
added value over balanced BW". Sub-studies: B1 bias/RMSE grid
(ε × understatement × τ²), B2 breakdown curve in the outlier shift, B3
coverage (analytic CIs vs bootstrap-over-studies for WFR).
