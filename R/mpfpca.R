# Forward `weights` / `min_no_pairs` to fpca.sgc.lat only if the installed
# SGCTools version accepts them (robust across SGCTools builds).
#' @noRd
.call_fpca_sgc_lat <- function(X, type, argvals, df, weights = NULL, min_no_pairs = NULL, ...) {
  fmls <- names(formals(fpca.sgc.lat))
  a <- list(X, type = type, argvals = argvals, df = df, ...)
  if (!is.null(weights)      && "weights"      %in% fmls) a$weights      <- weights
  if (!is.null(min_no_pairs) && "min_no_pairs" %in% fmls) a$min_no_pairs <- min_no_pairs
  do.call(fpca.sgc.lat, a)
}

#' Multivariate mixed-type FPCA (direct / partially separable estimator)
#'
#' M2FPCA: jointly estimates the marginal and cross-covariance surfaces of a set
#' of mixed-type functional variables and assembles a partially separable
#' multivariate Karhunen-Loeve representation. Marginal surfaces come from
#' [SGCTools::fpca.sgc.lat()]; cross surfaces from [pfpca_crosscov()]. Under
#' partial separability the variables share a common temporal eigenbasis; the
#' cross-domain dependence is captured by the per-component p-by-p score
#' covariance matrices in `S`.
#'
#' @param dat_list List of `p` numeric matrices (each n by m): the functional
#'   variables on a common grid. `NA`s allowed (irregular sampling).
#' @param type Character vector of length `p`, each `"cont"`/`"trunc"`/`"ord"`/`"bin"`.
#' @param argvals Numeric vector of length m of argument (time) values.
#' @param df Integer; B-spline degrees of freedom per margin.
#' @param marginal Either `"equal"` (default; shared eigenbasis from the averaged
#'   marginal covariance, i.e. the partial-separability assumption) or any other
#'   value to use each variable's own eigenbasis.
#' @param weights Logical; weight the NLS objective by the number of co-observed
#'   pairs (recommended, especially for irregular data).
#' @param min_no_pairs Minimum co-observed pairs per (s, t) cell.
#' @param ... Passed through to the marginal `fpca.sgc.lat` calls.
#'
#' @return A list with
#'   \describe{
#'     \item{pfpca_results}{list of per-variable eigenanalyses (`lambda`, `phi`, `cumFVE`)}
#'     \item{S}{list of length `L`; `S[[l]]` is the p-by-p covariance of the l-th
#'       component scores across variables}
#'     \item{L}{number of components}
#'     \item{raw_est}{the raw (mp)-by-(mp) block covariance estimate}
#'     \item{raw_est_cov}{delta-method variance of `raw_est` (vectorised)}
#'   }
#' @references Dey D., Ghosal R., Merikangas K., Zipunnikov V. (2024) \doi{10.1002/sim.10240}
#' @export
mpfpca.dir <- function(dat_list, type, argvals, df = 5, marginal = "equal",
                       weights = TRUE, min_no_pairs = 30, ...) {
  p <- length(dat_list)
  m <- ncol(dat_list[[1]])
  n <- nrow(dat_list[[1]])

  knotsT <- seq(min(argvals), max(argvals), l = df - 2)
  norder <- 4
  nbasisT <- length(knotsT) + norder - 2
  dayrngT <- c(min(argvals), max(argvals))
  bbasisT <- create.bspline.basis(dayrngT, nbasisT, norder, knotsT)
  bs.argvals <- eval.basis(argvals, bbasisT)

  cov_list <- pfpca_list <- list()
  Shat <- matrix(NA_real_, ncol = m * p, nrow = m * p)
  Shat_vcov <- matrix(0, ncol = (m * p)^2, nrow = (m * p)^2)
  MV <- matrix(1:(m * p)^2, ncol = m * p, nrow = m * p)

  # marginal covariance surfaces + their delta-method variances
  for (i in 1:p) {
    tp <- type[i]
    cov_list[[i]] <- .call_fpca_sgc_lat(dat_list[[i]], type = tp, weights = weights,
                                        argvals = argvals, df = df, min_no_pairs = min_no_pairs)
    Chat2 <- as.matrix(cov_list[[i]]$cov)
    D1 <- diag(as.numeric(ginv_deriv(Chat2))) %*%
      (kronecker(bs.argvals, diag(m)) %*% kronecker(diag(df), bs.argvals)) %*% up.elim.diag(df)
    V_c <- D1 %*% as.matrix(cov_list[[i]]$vcov) %*% t(D1)
    idx <- ((i - 1) * m + 1):(i * m)
    idx2 <- as.numeric(MV[idx, idx])
    Shat[idx, idx] <- Chat2
    Shat_vcov[idx2, idx2] <- as.matrix(V_c)
    pfpca_list[[i]] <- eigen_analysis(cov_list[[i]]$cov, argvals, FVEthreshold = 0.99, maxK = n * p)
  }

  # number of components (max across variables), then refit keeping all of them
  L <- max(sapply(pfpca_list, function(x) length(x$lambda)))
  for (i in 1:p) {
    pfpca_list[[i]] <- eigen_analysis(cov_list[[i]]$cov, argvals, FVEthreshold = 1, maxK = L)
  }

  if (marginal == "equal") {
    cov_avg <- Reduce("+", lapply(cov_list, function(x) as.matrix(x$cov))) / p
    cov_avg <- as.matrix(nearPD(cov_avg)$mat)
    pfpca_avg <- eigen_analysis(cov_avg, argvals, FVEthreshold = 1, maxK = L)
    for (i in 1:p) pfpca_list[[i]] <- pfpca_avg
    L <- length(pfpca_avg$lambda)
  }

  # per-component p x p score covariance matrices
  S_list <- vector("list", L)
  for (i in 1:L) S_list[[i]] <- matrix(NA_real_, ncol = p, nrow = p)

  for (j in 1:(p - 1)) {
    for (k in (j + 1):p) {
      pf_temp <- pfpca_crosscov(dat_list = list(dat_list[[j]], dat_list[[k]]),
                                type = c(type[j], type[k]), argvals = argvals,
                                df = df, weights = weights, min_no_pairs = min_no_pairs)
      Chat2 <- as.matrix(pf_temp$cov)
      D1 <- diag(as.numeric(ginv_deriv(Chat2))) %*%
        (kronecker(bs.argvals, diag(m)) %*% kronecker(diag(df), bs.argvals))
      V_c <- D1 %*% pf_temp$vcov %*% t(D1)
      jdx <- ((j - 1) * m + 1):(j * m)
      kdx <- ((k - 1) * m + 1):(k * m)
      jkdx1 <- as.numeric(MV[jdx, kdx])
      jkdx2 <- as.numeric(MV[kdx, jdx])
      Shat[jdx, kdx] <- Chat2
      Shat[kdx, jdx] <- t(Chat2)
      Shat_vcov[jkdx1, jkdx1] <- as.matrix(V_c)
      Shat_vcov[jkdx2, jkdx2] <- t(V_c)
      Cphi <- sapply(1:m, function(y) sapply(1:L, function(x) trapz(argvals, Chat2[y, ] * pfpca_list[[j]]$phi[, x])))
      sjk <- sapply(1:L, function(x) trapz(argvals, Cphi[x, ] * pfpca_list[[k]]$phi[, x]))
      for (i in 1:L) {
        S_list[[i]][j, j] <- pfpca_list[[j]]$lambda[i]
        S_list[[i]][j, k] <- S_list[[i]][k, j] <- sjk[i]
        S_list[[i]][k, k] <- pfpca_list[[k]]$lambda[i]
      }
    }
  }

  list(pfpca_results = pfpca_list, S = S_list, L = L,
       raw_est = Shat, raw_est_cov = Shat_vcov,
       cov_marginal = lapply(cov_list, function(x) as.matrix(x$cov)),
       argvals = argvals)
}

#' Reconstruct the multivariate covariance from an mpfpca.dir fit
#'
#' Rebuilds the full (mp)-by-(mp) covariance surface from the partially separable
#' representation returned by [mpfpca.dir()] (shared eigenfunctions `phi` and the
#' per-component cross-variable score covariances `S`).
#'
#' @param mpf Output of [mpfpca.dir()].
#' @return An (mp)-by-(mp) covariance matrix.
#' @export
cov.from.mpfpca.dir <- function(mpf) {
  p <- nrow(mpf$S[[1]])
  m <- nrow(mpf$pfpca_results[[1]]$phi)
  L <- mpf$L
  M <- matrix(NA_real_, ncol = m * p, nrow = m * p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      idx <- ((i - 1) * m + 1):(i * m)
      jdx <- ((j - 1) * m + 1):(j * m)
      M[idx, idx] <- Reduce("+", lapply(1:L, function(x) {
        pfi <- mpf$pfpca_results[[i]]; pfi$lambda[x] * (pfi$phi[, x] %*% t(pfi$phi[, x]))
      }))
      M[jdx, jdx] <- Reduce("+", lapply(1:L, function(x) {
        pfj <- mpf$pfpca_results[[j]]; pfj$lambda[x] * (pfj$phi[, x] %*% t(pfj$phi[, x]))
      }))
      M[idx, jdx] <- Reduce("+", lapply(1:L, function(x) {
        pfi <- mpf$pfpca_results[[i]]; pfj <- mpf$pfpca_results[[j]]
        sij <- mpf$S[[x]][i, j]; sij * pfi$phi[, x] %*% t(pfj$phi[, x])
      }))
      M[jdx, idx] <- t(M[idx, jdx])
    }
  }
  M
}
