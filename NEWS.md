# M2FPCA 0.1.3

* New `select_df_bic()`: selects the spline dimension of each marginal surface
  from a grid of candidates by the BIC of the paper (Web Appendix C.4),
  `n log|C| + sum_i V_i' C^{-1} V_i + K(K+1)/2 log n`, with the latent BLUP `V_i`.
* `mpfpca.dir()` gains `df_marginal` (per-variable dimensions, or `"bic"`) and
  `df_grid`. The default `df_marginal = NULL` keeps the previous behaviour.
  The fit now also returns `df_marginal` and, for `"bic"`, the BIC table.
* `mpfpca()` (ps-M2FPCA) accepts a different `df` for each variable (a vector of
  length `p`), and the same `df_marginal` / `df_grid` options, including `"bic"`.
  It now also returns `df_marginal` and `bic`.
* README: new section "Choosing the spline dimension".
* `mpfpca.dir()` now passes `...` to every marginal `fpca.sgc.lat()` call, as
  documented; previously it was ignored.
* The spline fits start from random values, so call `set.seed()` for
  reproducible results.

# M2FPCA 0.1.2

* Requires SGCTools >= 1.1.5. That release corrects `SGCTools::Kendall_mixed()`, which
  undercounted concordant and discordant pairs and pulled the sample Kendall's tau
  toward zero when few subjects were pairwise observed (about 15% at 10-15
  subjects, 3-6% at 24-30, under 1% at 100 or more). It now equals Kendall's tau_a
  exactly. `mpfpca()`, `mpfpca.dir()` and `pfpca_crosscov()` use it for every
  surface, so their estimates change slightly, most in sparse cells. SGCTools
  1.1.5 also restores explicit arguments of `bridgeF_tt()`, which 1.1.4 had broken
  for truncated components in `fpca.sgc.lat()`.

# M2FPCA 0.1.1

* `pfpca_crosscov()`, ordinal x truncated pairs: the truncation threshold is now
  `qnorm(P(X = 0))`, as in the other truncated bridges. It was computed as
  `qnorm(1 - mean(x))`, which is not invariant to the marginal transformation
  and biased the bridged correlation (true 0.5 -> 0.464 in a check with
  n = 20000; corrected 0.499). This changes M2FPCA estimates for any pair of an
  ordinal and a truncated variable.

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
