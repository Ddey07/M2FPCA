# Forward `ridge` to getLatentPreds only if the installed SGCTools accepts it.
#' @noRd
.call_getLatentPreds <- function(X, type, lat.cov.est, ridge = NULL, impute.missing = FALSE) {
  fmls <- names(formals(getLatentPreds))
  a <- list(X = X, type = type, lat.cov.est = lat.cov.est, impute.missing = impute.missing)
  if (!is.null(ridge) && "ridge" %in% fmls) a$ridge <- ridge
  do.call(getLatentPreds, a)
}

#' Multivariate principal-component scores from an mpfpca.dir fit
#'
#' Computes subject-level multivariate FPC scores by (1) predicting each
#' variable's latent continuous trajectory with [SGCTools::getLatentPreds()] and
#' (2) projecting those trajectories onto the (shared) temporal eigenfunctions
#' from an [mpfpca.dir()] fit:
#' \eqn{\theta_{i v \ell} = \int L_{v i}(t)\,\phi_\ell(t)\,dt}.
#'
#' @param dat_list List of `p` numeric matrices (n by m), as passed to [mpfpca.dir()].
#' @param type Character vector of length `p` of variable types.
#' @param mpf Output of [mpfpca.dir()].
#' @param argvals Numeric vector of length m of argument values (must match the fit).
#' @param npc Number of components to score; defaults to `mpf$L`.
#' @param ridge Optional ridge passed to `getLatentPreds` if the installed
#'   SGCTools supports it (stabilises the per-subject BLUP solve for poorly
#'   estimated margins). Ignored otherwise.
#' @param zcap Winsorising cap on the magnitude of latent predictions
#'   (`Inf` disables); guards against rare blow-ups before projection.
#' @param impute.missing Logical, passed to `getLatentPreds`.
#'
#' @return A list with
#'   \describe{
#'     \item{scores}{list of length `npc`; `scores[[l]]` is an n-by-p matrix of
#'       the l-th component scores for each variable}
#'     \item{scores_wide}{an n-by-(npc*p) matrix, components grouped}
#'     \item{latent}{list of `p` latent-trajectory matrices (n by m)}
#'   }
#' @export
mpfpca_scores <- function(dat_list, type, mpf, argvals, npc = NULL,
                          ridge = 0.05, zcap = 8, impute.missing = FALSE) {
  p <- length(dat_list)
  m <- ncol(dat_list[[1]])
  n <- nrow(dat_list[[1]])
  if (is.null(npc)) npc <- mpf$L
  npc <- min(npc, mpf$L)

  latent <- vector("list", p)
  for (v in 1:p) {
    # use the stored marginal latent correlation when available, else rebuild it
    lat_cov <- if (!is.null(mpf$cov_marginal)) {
      mpf$cov_marginal[[v]]
    } else {
      pf <- mpf$pfpca_results[[v]]
      pf$phi %*% diag(pf$lambda, length(pf$lambda)) %*% t(pf$phi)
    }
    z <- .call_getLatentPreds(X = dat_list[[v]], type = rep(type[v], m),
                              lat.cov.est = lat_cov,
                              ridge = ridge, impute.missing = impute.missing)
    if (is.finite(zcap)) z <- pmin(pmax(z, -zcap), zcap)
    latent[[v]] <- z
  }

  # center each variable by its cross-sectional mean (as fgm::pfpca does), then
  # project onto the temporal eigenfunctions: theta_{ivl} = int (L_vi - mu_v) phi_l
  mu_list <- lapply(latent, colMeans)
  scores <- vector("list", npc)
  for (l in 1:npc) {
    Sl <- matrix(NA, nrow = n, ncol = p)
    for (v in 1:p) {
      phi_l <- mpf$pfpca_results[[v]]$phi[, l]
      mu_v <- mu_list[[v]]
      Sl[, v] <- apply(latent[[v]], 1, function(row) trapz(argvals, (row - mu_v) * phi_l))
    }
    colnames(Sl) <- paste0("V", 1:p, "_PC", l)
    scores[[l]] <- Sl
  }
  scores_wide <- do.call(cbind, scores)

  list(scores = scores, scores_wide = scores_wide, latent = latent)
}

#' Covariance implied by a partially separable FPCA object
#'
#' Assembles the multivariate covariance \eqn{\sum_\ell \mathrm{cov}(\theta_\ell)
#' \otimes \phi_\ell \phi_\ell^\top} from a partially separable FPCA fit `pf`
#' (e.g. the output of `fgm::pfpca`), as used by the ps-M2FPCA variant.
#'
#' @param pf A partially separable FPCA object with elements `L` (number of
#'   components), `theta` (list of component scores) and `phi` (eigenfunctions).
#' @return The implied covariance matrix.
#' @export
sigma_from_fpca <- function(pf) {
  L <- pf$L
  Reduce("+", lapply(1:L, function(x) {
    sig <- cov(t(pf$theta[[x]]))
    phi <- pf$phi[x, ]
    kronecker(sig, phi %*% t(phi))
  }))
}

#' Partially separable M2FPCA (ps-M2FPCA)
#'
#' Two-stage variant: estimate each variable's marginal covariance and latent
#' trajectories, then run a partially separable FPCA (`fgm::pfpca`) on the
#' concatenated latent predictions. Requires the suggested package \pkg{fgm}.
#'
#' @inheritParams mpfpca_scores
#' @param df Integer; B-spline degrees of freedom per margin.
#' @param min_no_pairs Minimum number of co-observed pairs required per
#'   (s, t) cell in the marginal covariance fits.
#' @param ... Passed to the marginal `fpca.sgc.lat` calls.
#' @return A list with `pf` (the `fgm::pfpca` object) and `L` (list of latent
#'   trajectory matrices).
#' @export
mpfpca <- function(dat_list, type, argvals, df, min_no_pairs = 30,
                   ridge = 0.05, zcap = 8, ...) {
  if (!requireNamespace("fgm", quietly = TRUE)) {
    stop("ps-M2FPCA requires the 'fgm' package. Install it with install.packages('fgm').")
  }
  p <- length(dat_list)
  m <- ncol(dat_list[[1]])
  cov_list <- z_list <- vector("list", p)
  for (i in 1:p) {
    tp <- type[i]
    cov_list[[i]] <- .call_fpca_sgc_lat(dat_list[[i]], type = tp, argvals = argvals,
                                        df = df, min_no_pairs = min_no_pairs, ...)
    z <- .call_getLatentPreds(X = dat_list[[i]], type = rep(tp, m),
                              lat.cov.est = cov_list[[i]]$cov, ridge = ridge)
    if (is.finite(zcap)) z <- pmin(pmax(z, -zcap), zcap)
    z_list[[i]] <- z
  }
  partial_fsgc <- fgm::pfpca(z_list)
  list(pf = partial_fsgc, L = z_list)
}
