# =============================================================================
# 02_run_m2fpca.R
#
# Fit M2FPCA (direct estimator) to the synthetic mHealth data.
#
# For a real analysis this is a single call:
#
#     mpf <- mpfpca.dir(dat_list, type, argvals, df = 5, weights = TRUE)
#
# Here we run the SAME computation but checkpoint each marginal covariance and
# each pairwise cross-covariance to disk, so the (deliberately small but still
# nontrivial) fit can resume if interrupted. The final object is byte-for-byte
# the structure returned by mpfpca.dir(). Set SDA_ONECALL=1 to instead run the
# plain one-call version.
# =============================================================================

suppressMessages(library(M2FPCA))
outdir <- Sys.getenv("SDA_OUT", unset = "output")
st  <- readRDS(file.path(outdir, "synthetic_data.rds"))
dat_list <- st$dat_list; type <- st$type; argvals <- st$argvals
m <- st$m; p <- st$p; n <- st$n; df <- 5; min_no_pairs <- 15

## ---- one-call path (used when it fits in the available time) --------------
if (Sys.getenv("SDA_ONECALL", "") == "1") {
  mpf <- mpfpca.dir(dat_list, type = type, argvals = argvals, df = df,
                    weights = TRUE, min_no_pairs = min_no_pairs)
  saveRDS(mpf, file.path(outdir, "mpf.rds"))
  cat("mpfpca.dir one-call done. L =", mpf$L, "\n"); quit(save = "no")
}

## ---- checkpointed path ----------------------------------------------------
# internal helpers reused so the assembly matches mpfpca.dir() exactly
ea        <- getFromNamespace("eigen_analysis", "M2FPCA")
trapz     <- getFromNamespace("trapz", "M2FPCA")
ginv_der  <- getFromNamespace("ginv_deriv", "M2FPCA")
up_elim   <- getFromNamespace("up.elim.diag", "M2FPCA")
call_marg <- getFromNamespace(".call_fpca_sgc_lat", "M2FPCA")

knotsT <- seq(min(argvals), max(argvals), l = df - 2)
bbasisT <- fda::create.bspline.basis(c(min(argvals), max(argvals)),
                                     length(knotsT) + 2, 4, knotsT)
bs.argvals <- fda::eval.basis(argvals, bbasisT)

t0 <- Sys.time()
budget <- as.numeric(Sys.getenv("SDA_BUDGET", "600"))
time_left <- function() budget - as.numeric(Sys.time() - t0, units = "secs")
save_atomic <- function(obj, path) { tmp <- paste0(path, ".tmp"); saveRDS(obj, tmp); file.rename(tmp, path) }

## (1) marginal covariance surfaces --------------------------------------
marg_dir <- file.path(outdir, "marginals"); dir.create(marg_dir, showWarnings = FALSE)
for (i in 1:p) {
  f <- file.path(marg_dir, sprintf("marg_%d.rds", i))
  if (file.exists(f)) next
  if (time_left() < 20) { cat("PAUSE before marginal", i, "\n"); quit(status = 3) }
  ti <- Sys.time()
  cv <- call_marg(dat_list[[i]], type = type[i], weights = TRUE,
                  argvals = argvals, df = df, min_no_pairs = min_no_pairs)
  save_atomic(cv, f)
  cat(sprintf("  marginal %d (%s) done in %.1fs\n", i, type[i],
              as.numeric(Sys.time() - ti, units = "secs")))
}

## (2) pairwise cross-covariance surfaces --------------------------------
pair_dir <- file.path(outdir, "pairs"); dir.create(pair_dir, showWarnings = FALSE)
pairs <- t(utils::combn(p, 2))
for (r in seq_len(nrow(pairs))) {
  j <- pairs[r, 1]; k <- pairs[r, 2]
  f <- file.path(pair_dir, sprintf("pair_%d_%d.rds", j, k))
  if (file.exists(f)) next
  if (time_left() < 30) { cat("PAUSE before pair", j, k, "\n"); quit(status = 3) }
  ti <- Sys.time()
  pf <- pfpca_crosscov(dat_list = list(dat_list[[j]], dat_list[[k]]),
                       type = c(type[j], type[k]), argvals = argvals,
                       df = df, weights = TRUE, min_no_pairs = min_no_pairs)
  save_atomic(pf, f)
  cat(sprintf("  cross-cov %d-%d (%s,%s) done in %.1fs\n", j, k, type[j], type[k],
              as.numeric(Sys.time() - ti, units = "secs")))
}

## (3) assemble the mpfpca.dir() object ----------------------------------
cov_list <- lapply(1:p, function(i) readRDS(file.path(marg_dir, sprintf("marg_%d.rds", i))))

Shat <- matrix(NA_real_, m * p, m * p)
MV   <- matrix(1:(m * p)^2, ncol = m * p, nrow = m * p)
pfpca_list <- vector("list", p)
for (i in 1:p) {
  idx <- ((i - 1) * m + 1):(i * m)
  Shat[idx, idx] <- as.matrix(cov_list[[i]]$cov)
  pfpca_list[[i]] <- ea(cov_list[[i]]$cov, argvals, FVEthreshold = 0.99, maxK = n * p)
}
L <- max(sapply(pfpca_list, function(x) length(x$lambda)))
for (i in 1:p) pfpca_list[[i]] <- ea(cov_list[[i]]$cov, argvals, FVEthreshold = 1, maxK = L)

# partial-separability shared eigenbasis (marginal = "equal")
cov_avg <- Reduce("+", lapply(cov_list, function(x) as.matrix(x$cov))) / p
cov_avg <- as.matrix(Matrix::nearPD(cov_avg)$mat)
pfpca_avg <- ea(cov_avg, argvals, FVEthreshold = 1, maxK = L)
for (i in 1:p) pfpca_list[[i]] <- pfpca_avg
L <- length(pfpca_avg$lambda)

S_list <- lapply(1:L, function(i) matrix(NA_real_, p, p))
for (r in seq_len(nrow(pairs))) {
  j <- pairs[r, 1]; k <- pairs[r, 2]
  pf <- readRDS(file.path(pair_dir, sprintf("pair_%d_%d.rds", j, k)))
  Chat2 <- as.matrix(pf$cov)
  jdx <- ((j - 1) * m + 1):(j * m); kdx <- ((k - 1) * m + 1):(k * m)
  Shat[jdx, kdx] <- Chat2; Shat[kdx, jdx] <- t(Chat2)
  Cphi <- sapply(1:m, function(y) sapply(1:L, function(x) trapz(argvals, Chat2[y, ] * pfpca_list[[j]]$phi[, x])))
  sjk  <- sapply(1:L, function(x) trapz(argvals, Cphi[x, ] * pfpca_list[[k]]$phi[, x]))
  for (i in 1:L) {
    S_list[[i]][j, j] <- pfpca_list[[j]]$lambda[i]
    S_list[[i]][j, k] <- S_list[[i]][k, j] <- sjk[i]
    S_list[[i]][k, k] <- pfpca_list[[k]]$lambda[i]
  }
}

mpf <- list(pfpca_results = pfpca_list, S = S_list, L = L,
            raw_est = Shat,
            cov_marginal = lapply(cov_list, function(x) as.matrix(x$cov)),
            argvals = argvals)
saveRDS(mpf, file.path(outdir, "mpf.rds"))
cat(sprintf("M2FPCA assembled. L = %d | shared-eigenbasis (partial separability)\n", L))
cat("FIT_DONE\n")
