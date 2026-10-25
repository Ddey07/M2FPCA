## BIC selection of the marginal spline dimension (M2FPCA 0.1.3)

make_bic_data <- function(n = 300, m = 8, seed = 3) {
  set.seed(seed)
  tt <- seq(0, 1, length.out = m)
  C <- exp(-abs(outer(tt, tt, "-")) / 0.3)
  R <- kronecker(matrix(c(1, 0.4, 0.4, 1), 2), C)
  V <- matrix(stats::rnorm(n * 2 * m), n) %*% chol(R)
  ord <- apply(V[, 1:m], 2, function(z) findInterval(z, c(-0.6, 0.1, 0.6)))
  cont <- V[, (m + 1):(2 * m)]
  ord[matrix(stats::runif(n * m) < 0.3, n)] <- NA
  cont[matrix(stats::runif(n * m) < 0.3, n)] <- NA
  for (i in which(rowSums(!is.na(ord)) == 0)) ord[i, 1] <- findInterval(V[i, 1], c(-0.6, 0.1, 0.6))
  for (i in which(rowSums(!is.na(cont)) == 0)) cont[i, 1] <- V[i, m + 1]
  list(dat = list(ord = ord, cont = cont), type = c("ord", "cont"), argvals = tt, m = m)
}

test_that("select_df_bic returns a dimension from the grid with finite BIC", {
  d <- make_bic_data()
  sel <- suppressWarnings(select_df_bic(d$dat, d$type, d$argvals, df_grid = c(4, 5, 6), min_no_pairs = 10))
  expect_length(sel$df, 2)
  expect_true(all(sel$df %in% c(4, 5, 6)))
  expect_equal(nrow(sel$table), 6)
  for (j in 1:2) {
    tj <- sel$table[sel$table$var == j, ]
    expect_equal(sel$df[j], tj$df[which.min(tj$bic)])
  }
})

test_that("BIC data terms equal the Gaussian -2 log-likelihood for complete continuous data", {
  set.seed(4)
  m <- 8; tt <- seq(0, 1, length.out = m); n <- 250
  C0 <- exp(-abs(outer(tt, tt, "-")) / 0.3)
  X <- matrix(stats::rnorm(n * m), n) %*% chol(C0)
  fit <- suppressWarnings(SGCTools::fpca.sgc.lat(X, type = "cont", argvals = tt, df = 5, min_no_pairs = 10, weights = TRUE))
  b <- M2FPCA:::.bic_margin(X, "cont", m, fit$cov, 5, 0, 1e-3)
  C <- as.matrix(Matrix::nearPD(stats::cov2cor(as.matrix(fit$cov)), corr = TRUE, posd.tol = 1e-3)$mat)
  V <- SGCTools::getLatentPreds(X, type = rep("cont", m), lat.cov.est = C, impute.missing = TRUE)
  m2ll <- -2 * sum(mvtnorm::dmvnorm(V, sigma = C, log = TRUE)) - n * m * log(2 * pi)
  expect_equal(unname(b["logdet"] + b["quad"]), m2ll, tolerance = 1e-8)
  expect_equal(unname(b["bic"]), m2ll + 5 * 6 / 2 * log(n), tolerance = 1e-8)
})

test_that("subjects without observations are left out, and bad inputs are rejected", {
  d <- make_bic_data(seed = 6)
  d$dat$cont[1, ] <- NA
  set.seed(6)
  sel <- suppressWarnings(select_df_bic(d$dat["cont"], "cont", d$argvals, df_grid = c(4, 5), min_no_pairs = 10))
  expect_true(all(is.finite(sel$table$bic)))
  expect_true(all(sel$table$n_used == nrow(d$dat$cont) - 1))
  expect_error(select_df_bic(d$dat, d$type, d$argvals, df_grid = 3), "at least 4")
  expect_warning(try(select_df_bic(d$dat["cont"], "cont", d$argvals, df_grid = c(4, 9), min_no_pairs = 10), silent = TRUE),
                 "Skipping df = 9")
  expect_error(mpfpca.dir(d$dat, d$type, d$argvals, df = 4, df_marginal = "BIC"), "df_marginal")
  expect_error(mpfpca.dir(d$dat, d$type, d$argvals, df = 4, df_marginal = 3), "df_marginal")
})

test_that("mpfpca.dir: df_marginal = NULL reproduces the fixed-df fit, and 'bic' uses the selected dimensions", {
  d <- make_bic_data(seed = 5)
  set.seed(11); f0 <- suppressWarnings(mpfpca.dir(d$dat, d$type, d$argvals, df = 4, min_no_pairs = 10))
  set.seed(11); f1 <- suppressWarnings(mpfpca.dir(d$dat, d$type, d$argvals, df = 4, min_no_pairs = 10, df_marginal = 4))
  expect_identical(f0$raw_est, f1$raw_est)
  expect_equal(f0$df_marginal, c(4L, 4L))
  expect_null(f0$bic)
  set.seed(12); sel <- suppressWarnings(select_df_bic(d$dat, d$type, d$argvals, df_grid = c(4, 6), min_no_pairs = 10))
  set.seed(12); fb <- suppressWarnings(mpfpca.dir(d$dat, d$type, d$argvals, df = 4, min_no_pairs = 10,
                                                  df_marginal = "bic", df_grid = c(4, 6)))
  expect_equal(fb$df_marginal, sel$df)
  expect_false(is.null(fb$bic))
  m <- d$m
  for (j in 1:2) {
    ix <- ((j - 1) * m + 1):(j * m)
    expect_identical(fb$raw_est[ix, ix], as.matrix(sel$fits[[j]]$cov))   # fits are reused
  }
  ## the cross surface keeps df = 4 (fits start from random values, hence the tolerance)
  expect_equal(fb$raw_est[1:m, m + 1:m], f0$raw_est[1:m, m + 1:m], tolerance = 1e-3)
  expect_true(all(is.finite(fb$raw_est)))
})

test_that("mpfpca.dir accepts a different dimension for each margin", {
  d <- make_bic_data(seed = 7)
  set.seed(13)
  f <- suppressWarnings(mpfpca.dir(d$dat, d$type, d$argvals, df = 4, min_no_pairs = 10, df_marginal = c(4, 6)))
  expect_equal(f$df_marginal, c(4L, 6L))
  expect_true(all(is.finite(f$raw_est)))
  expect_true(all(is.finite(f$raw_est_cov)))
})
