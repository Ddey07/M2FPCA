# Synthetic mHealth analysis with **M2FPCA**

An end-to-end, fully reproducible pipeline that **mimics** the mobile-health data
analysis in Dey, Ghosal, Merikangas & Zipunnikov (2026), *Multivariate Functional
Principal Component Analysis for Mixed-Type mHealth Data: An Application to Mood
Disorders* ([arXiv:2603.11385](https://arxiv.org/abs/2603.11385)) — using only
**synthetic data**. No real participant data from the NIMH Family Study is used
or distributed here.

The pipeline exercises every exported function in the package on a mixed-type,
irregularly-sampled, multivariate functional dataset and reproduces the paper's
qualitative findings at a small, fast scale.

## What it mirrors from the paper

| Paper (Section 4)                                             | This synthetic pipeline                                  |
|--------------------------------------------------------------|----------------------------------------------------------|
| 4 domains: Sad Mood, Anxiousness, Energy (ordinal EMA) + Total Log Activity (continuous actigraphy) | Same 4 domains, same measurement types                   |
| 16 one-hour bins over the waking window                      | 12 bins (coarser, for speed)                             |
| Ordinal Likert collapsed to 6 levels                         | Ordinal collapsed to 3 levels (for speed)                |
| Sparse, asynchronous EMA + dense actigraphy                  | ~60% of EMA bins observed, ~96% of activity bins         |
| 307 participants, 4 diagnostic groups                        | 160 synthetic subject-days, 4 groups (40 each)           |
| Partially separable latent process, shared eigenbasis        | Same generative model (`marginal = "equal"`)             |
| 3 dominant components: level / morning–evening / circadian   | Same, recovered from the estimate                        |
| Scores → multinomial-LASSO digital phenotyping               | Scores → multinomial classifier of diagnostic group      |

## Files

- `01_generate_synthetic_data.R` — generates the mixed-type, irregularly sampled
  data from a known partially separable model with group-specific latent shifts.
- `02_run_m2fpca.R` — fits the direct estimator. Equivalent to a single
  `mpfpca.dir()` call (set `SDA_ONECALL=1`); by default it **checkpoints** each
  marginal (`fpca.sgc.lat` via the package) and each pairwise
  `pfpca_crosscov()` to disk so the fit is resumable.
- `03_scores_and_prediction.R` — `mpfpca_scores()` (digital biomarkers),
  `cov.from.mpfpca.dir()` reconstruction vs. analytic truth, and a multinomial
  prediction of diagnostic group with a train/test split.
- `04_figures.R` — renders the six figures below (base graphics only).
- `run_all.R` — runs all four steps in order.
- `analysis.Rmd` / `analysis.html` — a narrative, **step-by-step** report:
  data generation, data exploration, then the M²FPCA analysis (fit →
  eigenfunctions → covariance recovery → scores → prediction) with prose between
  every step. The rendered `analysis.html` is self-contained; open it in a
  browser. Re-knit with `rmarkdown::render("analysis.Rmd")` (the two slow steps
  cache to `analysis_cache/`; delete it to recompute from scratch).

## Run it

```r
# from this directory, with the M2FPCA package installed
Rscript run_all.R
```

Artifacts are written to `./output/*.rds`; figures to `./figures/*.png`.
Runtime is a few minutes (dominated by the ordinal truncated-normal integrals).

## Results (from the shipped run, seed-fixed)

**Shared components / variance explained** — three components capture the joint
diurnal structure, as in the paper:

| Component | Interpretation           | Cumulative FVE |
|-----------|--------------------------|----------------|
| PC1       | average diurnal level    | 49.3%          |
| PC2       | morning–evening contrast | 77.4%          |
| PC3       | circadian return         | 91.7%          |

**Covariance recovery** — the reconstructed 48×48 latent correlation matches the
analytic truth with `cor(vec(true), vec(est)) = 0.98` (Frobenius rel. error 0.23).

**Cross-domain PC1 score correlation** — recovers the generative coupling:
Mood–Anxiousness ≈ 0.52, activity negatively coupled to low mood/anxiety
(≈ −0.35), Energy–activity positive (≈ 0.29).

**Digital-biomarker prediction** — component scores predict diagnostic group at
**60% test accuracy vs. a 25% majority/chance baseline**. As in the paper, the
*activity morning–evening contrast* (TLAC PC2) separates Bipolar I from Bipolar
II, and *energy* components track MDD, while sad-mood components contribute
comparatively little.

## Figures

1. `fig1_example_trajectories.png` — example synthetic diurnal curves by domain and group.
2. `fig2_eigenfunctions.png` — estimated shared eigenfunctions + FVE.
3. `fig3_covariance.png` — true vs. estimated multivariate correlation.
4. `fig4_score_correlation.png` — cross-domain PC1 score correlation.
5. `fig5_biomarkers_by_group.png` — score distributions by diagnosis (digital biomarkers).
6. `fig6_prediction.png` — diagnostic-group confusion matrix.

## Caveats

Numbers are illustrative: sample size and grid are small and the generative
model is known, so recovery is easier than in the real study. The pipeline is a
faithful **software** demonstration of the method, not a clinical result.
