# M2FPCA build notes

This package was assembled from your research code, **without modifying
SGCTools**. Source material:

- `Sim_multivariate_pfpca/mpfpca_functions.R` (regular sampling)
- `Sim_multivariate_pfpca/S2_irregular/mpfpca_functions.R` (irregular sampling — newer)
- scoring logic adapted from `Data_analysis/multivariate/three_var_acti.R`

## What maps to what

| Research code | M2FPCA |
|---|---|
| `mpfpca.dir()` | `mpfpca.dir()` (exported, unchanged numerics; uses the general/irregular code path) |
| `pfpca_crosscov()` | `pfpca_crosscov()` (exported) |
| `cov.from.mpfpca.dir()` | `cov.from.mpfpca.dir()` (exported, **bug fixed** — see below) |
| `mpfpca()` (ps, two-stage) | `mpfpca()` (exported; needs `fgm`) |
| `sigma_from_fpca()` | `sigma_from_fpca()` (exported) |
| (scoring inlined in analysis script) | `mpfpca_scores()` (new exported wrapper) |

## Decisions taken (per our discussion)

1. **Both regular and irregular** sampling are covered by one code path: the
   irregular (more general) version, where regular/dense data is just the case
   with no `NA`s. `min_no_pairs` is the reliability cutoff; the
   `obj_df <- obj_df[complete.cases(obj_df), ]` assignment bug from the regular
   version (which crashed on sparse data) is fixed here.

2. **Self-contained, no fragile `:::` into SGCTools.** The SGCTools internal
   bridging helpers (`bridge_oo/co/bo/ot`, `bridgeF_bb/cc`, `ginv`, `Chat`,
   `make_function`) are **copied** into `R/internal-helpers.R`. SGCTools is still
   imported for the *exported* functions (`fpca.sgc.lat`, `getLatentPreds`,
   `Kendall_mixed`).

3. **External replacements (faithful ports).** `fgm:::.GetEigenAnalysisResults`
   is reproduced **line-for-line** as `eigen_analysis()` from the fgm source
   (https://github.com/javzapata/fgm, `R/pfpca.R`): same `gridSize =
   regGrid[2]-regGrid[1]`, eigenvalues `>= 0`, `maxK` truncation before the FVE
   threshold, trapezoidal eigenfunction normalisation `x/sqrt(int x^2)`, sign
   convention `sum(x * (1:m)) >= 0`, and `lambda = gridSize * eigenvalues`;
   `cumFVE` is on the 0-100 scale as in fgm. `fdapace::trapzRcpp` is reproduced
   as `trapz()` (identical trapezoidal rule). `mpfpca_scores()` uses the same
   projection as fgm's `pfpca` theta, `theta_{ivl} = int (L_vi - mu_v) phi_l`,
   with cross-sectional mean centering. So the package no longer needs
   `fgm`/`fdapace` for the core estimator; `fgm` is only *suggested*, for the
   optional ps-`mpfpca()`.

4. **Scope = core estimator + scoring** (no simulation DGP or multilevel code).

## Bugs fixed during the port

- `cov.from.mpfpca.dir()` referenced `mpf$pfpca_res`, but the field is
  `pfpca_results` — it would have errored. Fixed.
- Covariance averaging (`marginal = "equal"`) now uses
  `Reduce("+", lapply(..., as.matrix))` instead of `Reduce("+", sapply(...))`,
  which only worked by accident when `$cov` happened to be an S4 `Matrix`.
- Solver chatter silenced (`nls.lm.control(nprint=3)` / `trace=TRUE` removed).
- `deriv.b2` renamed to `deriv_b2` (avoids roxygen2 mistaking it for an S3
  method of `stats::deriv`, the same issue we hit in SGCTools).

## The four `mixedCCA:::` calls

`bridgeF_bt`, `bridgeF_bc`, `bridgeF_ct`, `bridgeF_tt` (binary–truncated,
binary–continuous, continuous–truncated, truncated–truncated bridges) live only
inside `mixedCCA` and have no SGCTools equivalent, so they are wrapped via
`mixedCCA:::` in `R/internal-helpers.R`. `R CMD check` will emit one NOTE about
this; it is intentional and isolated to those four one-line wrappers.

## ✅ Validation status (2026-07-02)

Built and checked under R 4.5.2 (SGCTools from GitHub, `Ddey07/SGCTools`).
`R CMD check` returns **1 NOTE** (the intentional `mixedCCA:::` calls), no
errors or warnings. Fixes applied while finalising:

- Replaced partial `matrix(nc=, nr=)` names with `ncol=`/`nrow=` (clears the
  "partial argument match" NOTE).
- Fixed `m <- ncol(dat_list[[2]])` → `[[1]]` in `mpfpca.dir()`.
- Documented the undocumented `min_no_pairs` argument of `mpfpca()`.
- Corrected an ordinal–ordinal `delta2` slice in `pfpca_crosscov()`
  (`(ncutoff1+2):(ncutoff1+ncutoff2+1)`).

`dev_smoke_test.R` passes end-to-end (dense + 15% MCAR). A full synthetic
mHealth analysis (EMA + actigraphy mimic of the paper) lives in
`inst/synthetic_data_analysis/`; on that run M2FPCA recovers the latent
correlation at `cor(vec)=0.98`, the three dominant diurnal components, and
scores that predict diagnostic group at 60% vs. a 25% baseline.

**Ordinal coding gotcha:** ordinal margins must be coded from **0**
(`0,1,2,…`). `SGCTools::getLatentPreds` indexes category cutoffs as
`cutoff[x+1]`/`cutoff[x+2]`, so 1-based codes yield `NA` in the
truncated-normal step.

## ⚠️ Original porting caveats (addressed above)

These were the open items before validation:

1. **`eigen_analysis()`** is now a faithful copy of fgm's
   `.GetEigenAnalysisResults` (see point 3 above), so `mpf$pfpca_results` should
   match an old `fgm`-based run to numerical precision. The one assumption it
   inherits from fgm is `gridSize = regGrid[2]-regGrid[1]` (uses the first grid
   gap as the spacing) — exact for equally spaced grids. Worth a quick
   `all.equal` check of `lambda`/`phi` against a saved `fgm` result if you have one.

2. **`mpfpca_scores()`** is a new wrapper (the research scoring was inline), but
   it now uses fgm's exact `pfpca` projection `theta_{ivl} = int (L_vi - mu_v)
   phi_l` on the latent predictions. The `ridge` argument is forwarded to
   `getLatentPreds` only if your installed SGCTools accepts it (public 1.1.3 does
   not, so it is silently ignored and `zcap` winsorising is the backstop).

## Build / check locally

`SGCTools` is now pulled from GitHub automatically (DESCRIPTION `Remotes:
Ddey07/SGCTools`); no local `SGCTools_dev` install is needed.

```r
# one-time: install the GitHub dependency (or let install_deps do it)
remotes::install_github("Ddey07/SGCTools")

setwd("path/to/M2FPCA")
devtools::document()      # regenerate NAMESPACE + man/ from roxygen
devtools::load_all()
source("dev_smoke_test.R")
devtools::check()         # remotes::install_deps(dependencies = TRUE) first if needed
```
