# =============================================================================
# 01_generate_synthetic_data.R
#
# Synthetic mHealth data that MIMICS (does not reproduce) the NIMH Family Study
# of Mood and Affective Spectrum Disorders analysed in Dey, Ghosal, Merikangas &
# Zipunnikov (2026), arXiv:2603.11385.  NO REAL PARTICIPANT DATA IS USED.
#
# Design mirrored from the paper (Section 4):
#   * p = 4 mixed-type functional domains, observed on a common daily grid of
#     m = 16 one-hour bins spanning the 07:00-22:00 waking window:
#         1. Sad Mood     - ordinal  (Likert, collapsed to 3 levels)
#         2. Anxiousness  - ordinal  (Likert, collapsed to 3 levels)
#         3. Energy       - ordinal  (Likert, collapsed to 3 levels)
#         4. TLAC         - continuous (Total Log Activity Count from actigraphy)
#   * Partially separable latent process: the domains SHARE a temporal
#     eigenbasis; cross-domain coupling enters through per-component score
#     covariances.  Three dominant modes, as reported in the paper:
#         phi_1  average diurnal level
#         phi_2  morning-evening contrast (homeostatic)
#         phi_3  morning-to-evening return (circadian)
#   * Four diagnostic groups (Control, MDD, Bipolar I, Bipolar II) whose latent
#     scores are shifted so that subject-level component scores act as
#     recoverable digital biomarkers (mirrors the multinomial-LASSO finding that
#     diurnal activity contrast separates BP-I from BP-II, energy variability
#     tracks MDD, and sad-mood fluctuations add little).
#
# Small sample size on purpose (demonstration, not inference).
# =============================================================================

set.seed(20240701)

## ---- dimensions ----------------------------------------------------------
# The paper uses a 16-bin grid with Likert responses collapsed to 6 levels.
# For a fast, self-contained demonstration we use a coarser 12-bin grid and
# 3-level ordinal responses; the estimator and interpretation are identical.
m       <- 12                       # ~1.25-hour bins, 07:00-22:00
argvals <- seq(0, 1, length = m)    # time-of-day, rescaled to [0, 1]
hours   <- seq(7, 22, length = m)   # clock hours (for plotting)
p       <- 4
domains <- c("Sad Mood", "Anxiousness", "Energy", "TLAC (activity)")
type    <- c("ord", "ord", "ord", "cont")

groups     <- c("Control", "MDD", "Bipolar I", "Bipolar II")
n_per_grp  <- 40                    # subject-days per group (small sample)
n          <- n_per_grp * length(groups)
group      <- factor(rep(groups, each = n_per_grp), levels = groups)

## ---- shared temporal eigenfunctions (partial separability) ---------------
# Orthonormal-ish cosine modes on the grid; nonzero at endpoints so no bin is
# degenerate. Interpreted as level / morning-evening contrast / circadian.
phi <- cbind(
  level     = rep(1, m),
  contrast  = sqrt(2) * cos(1 * pi * argvals),
  circadian = sqrt(2) * cos(2 * pi * argvals),
  minor     = sqrt(2) * cos(3 * pi * argvals)
)
# normalise each mode to unit trapezoidal norm
trapz <- function(x, y) sum((y[-1] + y[-length(y)]) / 2 * diff(x))
for (l in 1:ncol(phi)) phi[, l] <- phi[, l] / sqrt(trapz(argvals, phi[, l]^2))
Lc     <- ncol(phi)
lambda <- c(1.0, 0.55, 0.30, 0.12)   # component variances (decaying)

## ---- cross-domain coupling per component ---------------------------------
# Order of domains: 1 SadMood, 2 Anx, 3 Energy, 4 TLAC.
# Positive mood-anx-energy block; activity negatively coupled to low mood/anx,
# positively to energy -- a plausible affect/behaviour structure.
Rbase <- matrix(c( 1.0,  0.6,  0.4, -0.3,
                   0.6,  1.0,  0.3, -0.3,
                   0.4,  0.3,  1.0,  0.5,
                  -0.3, -0.3,  0.5,  1.0), 4, 4, byrow = TRUE)
# ensure valid correlation matrix
Rbase <- as.matrix(Matrix::nearPD(Rbase, corr = TRUE)$mat)

## ---- group-specific latent-score mean shifts (the digital-biomarker signal)
# Rows = groups, columns = (domain, component) effects we choose to encode.
# Component index: 1 level, 2 contrast, 3 circadian.
# Encode paper-style contrasts:
#   * Activity morning-evening CONTRAST (domain 4, comp 2): separates BP-I/BP-II
#   * Energy CIRCADIAN variability (domain 3, comp 3): elevated in MDD
#   * overall Energy LEVEL (domain 3, comp 1): lower in MDD/BP
shift <- array(0, dim = c(length(groups), p, Lc),
               dimnames = list(groups, domains, colnames(phi)))
# activity contrast: BP-I high, BP-II low
shift["Bipolar I",  "TLAC (activity)", "contrast"]  <-  1.1
shift["Bipolar II", "TLAC (activity)", "contrast"]  <- -1.1
# energy level: reduced in MDD and both bipolar
shift["MDD",        "Energy", "level"]     <- -1.0
shift["Bipolar I",  "Energy", "level"]     <- -0.5
shift["Bipolar II", "Energy", "level"]     <- -0.6
# energy circadian variability: MDD elevated (encoded as a mean circadian shift)
shift["MDD",        "Energy", "circadian"] <-  0.9
# sad-mood level slightly elevated in MDD (kept modest: "adds little" once others in)
shift["MDD",        "Sad Mood", "level"]   <-  0.5

## ---- draw latent scores and build latent trajectories --------------------
# theta[l] is n x p (scores of component l across the 4 domains)
Ldim  <- array(0, c(n, m, p))       # latent Gaussian trajectories
theta_true <- vector("list", Lc)
for (l in 1:Lc) {
  Sig  <- lambda[l] * Rbase
  z    <- mvtnorm::rmvnorm(n, sigma = Sig)          # n x p
  mu   <- t(sapply(seq_len(n), function(i) shift[as.integer(group[i]), , l]))
  th   <- z + mu * sqrt(lambda[l])                  # add group shift
  theta_true[[l]] <- th
  for (v in 1:p) Ldim[, , v] <- Ldim[, , v] + outer(th[, v], phi[, l])
}

## ---- map latent trajectories to observed mixed-type scales ---------------
zcol <- function(M) apply(M, 2, function(col) (col - mean(col)) / stats::sd(col))

to_ordinal <- function(L, levels = 3) {
  # cut standardized latent into ordinal categories coded 0..(levels-1).
  # The SGCTools latent-prediction step keys category cutoffs off a 0-based
  # code (cutoff[x + 1], cutoff[x + 2]), so ordinal margins MUST start at 0.
  Z <- zcol(L)
  qs <- stats::qnorm(seq(0, 1, length = levels + 1))
  apply(Z, 2, function(col) as.integer(cut(col, breaks = qs, include.lowest = TRUE)) - 1L)
}

SadMood <- to_ordinal(Ldim[, , 1], levels = 3)
Anx     <- to_ordinal(Ldim[, , 2], levels = 3)
Energy  <- to_ordinal(Ldim[, , 3], levels = 3)
# TLAC: continuous actigraphy-like signal (log activity counts), keep latent + offset
TLAC    <- 5 + 0.8 * Ldim[, , 4] + matrix(stats::rnorm(n * m, 0, 0.3), n, m)

dat_list <- list(SadMood = SadMood, Anx = Anx, Energy = Energy, TLAC = TLAC)

## ---- introduce EMA-style irregular (asynchronous) sampling ---------------
# Ordinal EMA domains are sparsely observed (few prompts/day); actigraphy dense.
# Mimics the paper's asynchronous mHealth regime. NAs = unobserved bins.
set.seed(99)
ema_keep  <- 0.60     # ~7 of 12 bins observed per EMA domain-day
for (v in 1:3) {
  M <- dat_list[[v]]
  miss <- matrix(stats::runif(n * m) > ema_keep, n, m)
  M[miss] <- NA
  # guarantee every column keeps all categories (estimator keys off col 1 cats)
  dat_list[[v]] <- M
}
# actigraphy: dense, small MCAR
acti_miss <- matrix(stats::runif(n * m) > 0.95, n, m)
dat_list$TLAC[acti_miss] <- NA

## ---- persist ------------------------------------------------------------
outdir <- Sys.getenv("SDA_OUT", unset = "output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(dat_list = dat_list, type = type, argvals = argvals, hours = hours,
             m = m, p = p, n = n, domains = domains, groups = groups,
             group = group, phi = phi, lambda = lambda, Rbase = Rbase,
             theta_true = theta_true, shift = shift),
        file.path(outdir, "synthetic_data.rds"))

cat(sprintf("Generated synthetic mHealth data: n=%d subject-days, p=%d domains, m=%d bins\n",
            n, p, m))
cat("Domains / types:\n"); print(data.frame(domain = domains, type = type))
cat(sprintf("Observed fraction per domain: %s\n",
            paste(sprintf("%s=%.2f", c("Mood","Anx","Energy","TLAC"),
                          sapply(dat_list, function(X) mean(!is.na(X)))), collapse = "  ")))
cat("Saved:", file.path(outdir, "synthetic_data.rds"), "\n")
