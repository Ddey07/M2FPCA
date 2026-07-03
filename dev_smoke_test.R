# dev_smoke_test.R --------------------------------------------------------
# Manual smoke test for the M2FPCA package.
# Run after: devtools::document(); devtools::load_all()
# Not part of the package build (see .Rbuildignore).
# -------------------------------------------------------------------------

set.seed(1)
n <- 150; m <- 10; p <- 3
argvals <- seq(0, 1, length = m)

## Partially separable DGP: shared Fourier eigenfunctions, cross-correlated
## component scores, then mapped to continuous / binary / ordinal margins.
Lc  <- 4
# cosine eigenfunctions: nonzero at the grid endpoints (unlike sin(l*pi*t),
# which is 0 at t = 0 and would make every variable's first time-point degenerate)
phi <- sapply(1:Lc, function(l) sqrt(2) * cos(l * pi * argvals))   # m x Lc
lambda <- c(1, 0.6, 0.3, 0.15)
Rcross <- matrix(0.5, p, p); diag(Rcross) <- 1

V <- array(0, c(n, m, p))
for (l in 1:Lc) {
  th <- mvtnorm::rmvnorm(n, sigma = lambda[l] * Rcross)   # n x p scores
  for (v in 1:p) V[, , v] <- V[, , v] + outer(th[, v], phi[, l])
}

# standardize columns before discretizing so every time-point has spread
zcol <- function(M) apply(M, 2, function(col) (col - mean(col)) / sd(col))
dat_list <- list(
  cont = V[, , 1],
  bin  = (zcol(V[, , 2]) > 0) * 1,
  ord  = apply(zcol(V[, , 3]), 2, function(col) as.numeric(cut(col, c(-Inf, -0.6, 0.6, Inf))) - 1)
)
type <- c("cont", "bin", "ord")
# the ordinal estimator keys the category count off column 1; ensure all
# columns show the full set of 3 categories
stopifnot(all(apply(dat_list$ord, 2, function(c) length(unique(c))) == 3))

cat("=== mpfpca.dir (regular, dense) ===\n")
mpf <- mpfpca.dir(dat_list, type = type, argvals = argvals, df = 5, weights = TRUE)
stopifnot(length(mpf$S) == mpf$L,
          nrow(mpf$raw_est) == m * p,
          all(dim(mpf$pfpca_results[[1]]$phi) == c(m, mpf$L)))
cat("L =", mpf$L, "| raw_est:", paste(dim(mpf$raw_est), collapse = "x"), "OK\n")

cat("\n=== cov.from.mpfpca.dir (reconstruction) ===\n")
Chat <- cov.from.mpfpca.dir(mpf)
stopifnot(all(dim(Chat) == c(m * p, m * p)),
          max(abs(Chat - t(Chat))) < 1e-8)
cat("reconstructed cov:", paste(dim(Chat), collapse = "x"), "symmetric OK\n")

cat("\n=== single cross-covariance surface (cont vs bin) ===\n")
cc <- pfpca_crosscov(list(dat_list$cont, dat_list$bin), type = c("cont", "bin"),
                     argvals = argvals, df = 5)
stopifnot(all(dim(cc$cov) == c(m, m)))
cat("cross-cov surface:", paste(dim(cc$cov), collapse = "x"), "OK\n")

cat("\n=== mpfpca_scores ===\n")
sc <- mpfpca_scores(dat_list, type, mpf, argvals, npc = 3)
stopifnot(length(sc$scores) == 3,
          all(dim(sc$scores[[1]]) == c(n, p)),
          all(dim(sc$scores_wide) == c(n, 3 * p)))
cat("scores_wide:", paste(dim(sc$scores_wide), collapse = "x"), "OK\n")
cat("cross-variable score correlation (component 1):\n")
print(round(cor(sc$scores[[1]]), 2))   # should pick up the 0.5 cross-correlation

cat("\n=== true vs estimated latent correlation ===\n")
# Analytic truth: full latent covariance = Rcross %x% K, K = sum_l lambda_l phi_l phi_l',
# so the true latent CORRELATION is Rcross %x% cov2cor(K).
K_true  <- Reduce("+", lapply(1:Lc, function(l) lambda[l] * outer(phi[, l], phi[, l])))
R_true  <- kronecker(Rcross, cov2cor(K_true))     # (mp x mp) true latent correlation
R_hat   <- cov.from.mpfpca.dir(mpf)               # estimated latent correlation
rel_err <- norm(R_hat - R_true, "F") / norm(R_true, "F")
vcor    <- cor(as.numeric(R_hat), as.numeric(R_true))
cat(sprintf("  relative Frobenius error: %.3f | cor(vec(true), vec(est)): %.3f\n", rel_err, vcor))
stopifnot(vcor > 0.9)
cat("  M2FPCA recovers the true latent correlation.  OK\n")

cat("\n=== irregular sampling (15% MCAR) ===\n")
dat_miss <- lapply(dat_list, function(X) { X[sample(length(X), 0.15 * length(X))] <- NA; X })
mpf2 <- mpfpca.dir(dat_miss, type = type, argvals = argvals, df = 5,
                   weights = TRUE, min_no_pairs = 20)
stopifnot(length(mpf2$S) == mpf2$L)
cat("irregular L =", mpf2$L, "OK\n")

## =========================================================================
## Fidelity check: eigen_analysis() vs the real fgm:::.GetEigenAnalysisResults
## (only runs if fgm is installed; eigen_analysis is internal, hence :::)
## =========================================================================
cat("\n=== fidelity: eigen_analysis() vs fgm:::.GetEigenAnalysisResults ===\n")
if (requireNamespace("fgm", quietly = TRUE)) {
  ea  <- getFromNamespace("eigen_analysis", "M2FPCA")
  fgm_ea <- getFromNamespace(".GetEigenAnalysisResults", "fgm")

  # a positive-definite covariance surface on the same grid
  Cmat <- as.matrix(mpf$cov_marginal[[1]])

  mine <- ea(Cmat, argvals, FVEthreshold = 0.99, maxK = nrow(Cmat))
  theirs <- fgm_ea(Cmat, argvals,
                   optns = list(FVEthreshold = 0.99, maxK = nrow(Cmat), verbose = FALSE))

  K <- length(mine$lambda)
  lam_ok <- isTRUE(all.equal(as.numeric(mine$lambda),
                             as.numeric(theirs$lambda[seq_len(K)]), tolerance = 1e-8))
  # eigenfunctions can differ by an overall sign per component
  phi_ok <- all(sapply(seq_len(K), function(k) {
    a <- mine$phi[, k]; b <- theirs$phi[seq_len(nrow(Cmat)), k]
    isTRUE(all.equal(a, b, tolerance = 1e-6)) ||
      isTRUE(all.equal(a, -b, tolerance = 1e-6))
  }))
  cat(sprintf("  components compared: %d | eigenvalues match: %s | eigenfunctions match: %s\n",
              K, lam_ok, phi_ok))
  stopifnot(lam_ok, phi_ok)
  cat("  eigen_analysis() reproduces fgm exactly.  OK\n")
} else {
  cat("  fgm not installed -- skipping fidelity check (install.packages('fgm') to enable).\n")
}

cat("\nAll M2FPCA smoke tests passed.\n")
