## End-to-end checks on sparse, asynchronous mixed data (M2FPCA 0.1.1).
## They exercise the SGCTools >= 1.1.4 latent-predictor fix: with sparse EMA-like
## sampling an ordinal column can miss a level, which previously made
## getLatentPreds() fail with "'upper' not specified or contains NA".

make_data <- function(n = 300, m = 6, obs = 0.5, seed = 1) {
  set.seed(seed)
  tt <- seq(0, 1, length.out = m)
  C <- exp(-abs(outer(tt, tt, "-")) / 0.5)
  R <- kronecker(matrix(c(1, 0.5, 0.5, 1), 2), C)
  V <- matrix(stats::rnorm(n * 2 * m), n) %*% chol(R)
  ord <- apply(V[, 1:m], 2, function(z) findInterval(z, c(-0.6, 0.1, 0.6)))   # 0..3
  cont <- V[, (m + 1):(2 * m)]
  ema <- matrix(stats::runif(n * m) < obs, n)                                # sparse EMA-like mask
  ord[!ema] <- NA
  ord[ord[, 2] == 1 & !is.na(ord[, 2]), 2] <- NA                             # column 2 misses level 1
  for (i in which(rowSums(!is.na(ord)) == 0)) ord[i, sample.int(m, 1)] <- findInterval(V[i, 1], c(-0.6, 0.1, 0.6))
  cont[matrix(stats::runif(n * m) < 0.05, n)] <- NA                           # dense "activity"
  list(dat = list(ord = ord, cont = cont), type = c("ord", "cont"), argvals = tt, R = R, m = m)
}

test_that("SGCTools >= 1.1.4 is in use", {
  expect_true(utils::packageVersion("SGCTools") >= "1.1.4")
})

test_that("mpfpca.dir + mpfpca_scores run on sparse ordinal data with a missing level", {
  d <- make_data()
  mpf <- suppressWarnings(mpfpca.dir(d$dat, d$type, d$argvals, df = 4, marginal = "equal", min_no_pairs = 10))
  expect_equal(dim(mpf$raw_est), c(2 * d$m, 2 * d$m))
  expect_true(all(is.finite(mpf$raw_est)))
  sc <- mpfpca_scores(d$dat, d$type, mpf, d$argvals, npc = 2)
  expect_length(sc$scores, 2)
  expect_true(all(is.finite(sc$scores_wide)))
})

test_that("joint (stacked) latent prediction works with the M2FPCA covariance", {
  d <- make_data(seed = 2)
  mpf <- suppressWarnings(mpfpca.dir(d$dat, d$type, d$argvals, df = 4, marginal = "diff", min_no_pairs = 10))
  S <- as.matrix(Matrix::nearPD(stats::cov2cor(mpf$raw_est), corr = TRUE, posd.tol = 1e-3)$mat)
  X <- do.call(cbind, d$dat)
  L <- SGCTools::getLatentPreds(X, type = rep(d$type, each = d$m), lat.cov.est = S, impute.missing = TRUE)
  expect_true(all(is.finite(L)))
})

test_that("ps-M2FPCA runs on the same sparse data", {
  skip_if_not_installed("fgm")
  d <- make_data(seed = 3)
  ps <- suppressWarnings(mpfpca(d$dat, d$type, d$argvals, df = 4, min_no_pairs = 10))
  S <- sigma_from_fpca(ps$pf)
  expect_equal(dim(S), c(2 * d$m, 2 * d$m))
  expect_true(all(is.finite(S)))
})
