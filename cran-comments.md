# cran-comments

## Resubmission (0.1.1)

This is a resubmission addressing the review comments of 2026-07-16
(Konstanze Lauseker) on gtmeta 0.1.0:

* "Please always explain all acronyms in the description text. -> FE,
  RE, UWLS"
  The Description now spells these out: "classical fixed-effect (FE),
  random-effects (RE), and unrestricted weighted least squares (UWLS)
  benchmarks". IGMI was already expanded.

* "If there are references describing the methods in your package,
  please add these in the description field"
  Added, in the requested format: Agueh and Carlier (2011)
  <doi:10.1137/100805741>, Chizat, Peyre, Schmitzer, and Vialard (2018)
  <doi:10.1007/s10208-016-9331-y>, DerSimonian and Laird (1986)
  <doi:10.1016/0197-2456(86)90046-2>, and Stanley and Doucouliagos
  (2015) <doi:10.1002/sim.6481>.

* "Missing Rd-tags: metric_ci_width.Rd: \value, metric_mse.Rd: \value"
  Both help pages now document the return value (a single numeric
  value in each case). All exported functions now have \value tags.

* "Please do not install packages in your functions, examples or
  vignette."
  The only installing code is `wfr_setup()`, a dedicated opt-in
  installer for the optional Python bridge (in the spirit of
  `reticulate::install_python()`). It is never called by the package
  itself, its examples, tests, or the vignette, so nothing is installed
  during checks. It now additionally requires explicit user consent: an
  interactive confirmation prompt, or `confirm = TRUE` in
  non-interactive sessions; otherwise it stops without installing
  anything. Its documentation states clearly that it installs Python
  packages into a virtualenv.

## R CMD check results

0 errors | 0 warnings | 0 notes (new submission NOTE excepted).
