# M2FPCA 0.1.1

* Requires SGCTools >= 1.1.4. That release fixes `SGCTools::getLatentPreds()`, which
  failed with "'upper' not specified or contains NA" whenever a sparsely observed
  ordinal column did not show every level (typical of EMA sampling). The latent
  predictions used by `mpfpca()` (ps-M2FPCA) and `mpfpca_scores()`, and any joint
  prediction with an M2FPCA covariance, are unchanged wherever the old code ran.
* `pfpca_crosscov()` (used by `mpfpca.dir()`): ordinal cutoffs are built over the
  full level set taken from all columns. The number of levels was read from
  column 1 only and the per-column cutoffs from the levels observed in that
  column, which failed or misaligned when a sparse column missed a level.
  Identical results whenever every column shows every level.
* `mpfpca()` (ps-M2FPCA) now requests imputed latent trajectories
  (`impute.missing = TRUE`), as in the research code behind the paper. The port
  used `FALSE`, which leaves `NA` at unobserved cells and breaks the partially
  separable step on irregularly sampled data. Dense data are unaffected.
* `mpfpca_scores()`: `impute.missing` now defaults to `TRUE` for the same reason.
* Added `tests/testthat` end-to-end checks on sparse, asynchronous mixed data.

# M2FPCA 0.1.0

* Initial release.
