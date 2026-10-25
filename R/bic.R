#' Select the spline dimension of each marginal surface by BIC
#'
#' For each functional variable, fits the marginal latent correlation surface
#' with [SGCTools::fpca.sgc.lat()] at every spline dimension in `df_grid` and
#' returns the dimension that minimises
#' \deqn{BIC(K) = n \log|\widehat C_K| + \sum_{i=1}^n \widehat V_i^\top
#'   \widehat C_K^{-1} \widehat V_i + \frac{K(K+1)}{2}\log n,}
#' the criterion of Dey, Ghosal, Merikangas and Zipunnikov (2026, Web Appendix C.4).
#'
#' Details of the implementation:
#' * \eqn{\widehat C_K} is the fitted latent correlation of the variable,
#'   projected to the nearest correlation matrix whose eigenvalues are at least
#'   `proj_tol` times the largest (as used for prediction). For large grids the
#'   value of `proj_tol` affects the size of the log-determinant term.
#' * \eqn{\widehat V_i} is the latent BLUP of subject \eqn{i} from this variable
#'   alone, from [SGCTools::getLatentPreds()] with `impute.missing = TRUE` and the
#'   given `ridge`. For binary, ordinal and censored truncated values it is a
#'   conditional mean, so the criterion is an approximation to the likelihood.
#' * Subjects with no observation of the variable are left out, and \eqn{n} is
#'   the number of subjects used.
#' * Candidates with more spline coefficients, \eqn{K(K+1)/2}, than there are
#'   off-diagonal cells, \eqn{m(m-1)/2}, are skipped. A candidate whose fit or
#'   BLUP fails gets `NA` and a warning; ties go to the smaller dimension.
#' * The fits start from random values (in `fpca.sgc.lat`), so call
#'   [set.seed()] first for reproducible results.
#'
#' @param dat_list List of `p` numeric matrices (each n by m), `NA` for unobserved cells.
#' @param type Character vector of length `p`, each `"cont"`/`"trunc"`/`"ord"`/`"bin"`.
#' @param argvals Numeric vector of length m of argument (time) values.
#' @param df_grid Integer vector of candidate spline dimensions (each at least 4).
#' @param weights,min_no_pairs Passed to `fpca.sgc.lat`, as in [mpfpca.dir()].
#' @param ridge Ridge passed to `getLatentPreds` for the BLUP.
#' @param proj_tol Relative eigenvalue floor of the nearest-correlation projection.
#' @param verbose Logical; print progress.
#' @param ... Further arguments passed to `fpca.sgc.lat`.
#'
#' @return A list with
#'   \describe{
#'     \item{df}{integer vector of length `p`: the selected dimension for each
#'       variable (`NA` if no candidate could be evaluated)}
#'     \item{table}{data frame with one row per variable and candidate: `var`,
#'       `type`, `df`, `bic`, its two data terms `logdet` (\eqn{n\log|\widehat C_K|})
#'       and `quad` (the quadratic form), `n_used`, and `seconds`}
#'     \item{fits}{list of length `p` with the `fpca.sgc.lat` fit at the
#'       selected dimension}
#'   }
#' @seealso [mpfpca.dir()], whose `df_marginal = "bic"` calls this function.
#' @references Dey, D., Ghosal, R., Merikangas, K., & Zipunnikov, V. (2026)
#'   "Multivariate Functional Principal Component Analysis for Mixed-Type mHealth
#'   Data: An Application to Mood Disorders." \emph{arXiv:2603.11385}.
#'   \url{https://arxiv.org/abs/2603.11385}
#' @examples
#' \dontrun{
#' set.seed(1)
#' sel <- select_df_bic(dat_list, type, argvals, df_grid = c(4, 5, 6, 8, 10), min_no_pairs = 30)
#' sel$df
#' fit <- mpfpca.dir(dat_list, type, argvals, df = 5, df_marginal = sel$df)
#' }
#' @export
select_df_bic <- function(dat_list, type, argvals, df_grid = c(4, 5, 6, 8, 10),
                          weights = TRUE, min_no_pairs = 30, ridge = 0.05,
                          proj_tol = 1e-3, verbose = FALSE, ...) {
  p <- length(dat_list)
  if (length(type) != p) stop("'type' must have one entry per element of 'dat_list'.")
  if (!is.numeric(df_grid) || anyNA(df_grid) || any(df_grid != round(df_grid)) || any(df_grid < 4))
    stop("'df_grid' must contain whole numbers of at least 4.")
  df_grid <- sort(unique(as.integer(df_grid)))
  m <- ncol(dat_list[[1]])
  ok_grid <- df_grid[df_grid * (df_grid + 1) / 2 <= m * (m - 1) / 2]
  if (length(ok_grid) < length(df_grid))
    warning("Skipping df = ", paste(setdiff(df_grid, ok_grid), collapse = ", "),
            ": more spline coefficients than off-diagonal cells (m = ", m, ").", call. = FALSE)
  if (!length(ok_grid)) stop("No candidate in 'df_grid' can be fitted with m = ", m, " grid points.")
  rows <- list(); best <- rep(NA_integer_, p); fits <- vector("list", p)
  for (j in seq_len(p)) {
    best_bic <- Inf
    for (K in ok_grid) {
      t0 <- unname(proc.time()[3])
      fit <- tryCatch(.call_fpca_sgc_lat(dat_list[[j]], type = type[j], argvals = argvals, df = K,
                                         weights = weights, min_no_pairs = min_no_pairs, ...),
                      error = function(e) {
                        warning(sprintf("select_df_bic: variable %d (%s), df = %d: fit failed: %s",
                                        j, type[j], K, conditionMessage(e)), call. = FALSE)
                        NULL })
      b <- if (is.null(fit)) c(bic = NA_real_, logdet = NA_real_, quad = NA_real_, n = NA_real_) else
        .bic_margin(dat_list[[j]], type[j], m, fit$cov, K, ridge, proj_tol, j)
      rows[[length(rows) + 1]] <- data.frame(var = j, type = type[j], df = K, bic = b[["bic"]],
                                             logdet = b[["logdet"]], quad = b[["quad"]], n_used = b[["n"]],
                                             seconds = round(unname(proc.time()[3]) - t0, 2))
      if (verbose) message(sprintf("variable %d (%s), df = %d: BIC = %s", j, type[j], K, format(b[["bic"]])))
      if (is.finite(b[["bic"]]) && b[["bic"]] < best_bic) { best_bic <- b[["bic"]]; best[j] <- K; fits[[j]] <- fit }
    }
    if (is.na(best[j])) warning(sprintf("select_df_bic: no candidate could be evaluated for variable %d (%s).",
                                        j, type[j]), call. = FALSE)
  }
  tab <- do.call(rbind, rows); rownames(tab) <- NULL
  list(df = best, table = tab, fits = fits)
}

## BIC of one marginal fit and its parts; NA (with a warning) when the projection or the BLUP fails.
#' @noRd
.bic_margin <- function(X, tp, m, cov, K, ridge, proj_tol, j = NA) {
  na <- c(bic = NA_real_, logdet = NA_real_, quad = NA_real_, n = NA_real_)
  fail <- function(what, e) {
    warning(sprintf("select_df_bic: variable %s (%s), df = %d: %s failed: %s", j, tp, K, what,
                    conditionMessage(e)), call. = FALSE)
    NULL }
  C <- tryCatch(as.matrix(nearPD(stats::cov2cor(as.matrix(cov)), corr = TRUE, posd.tol = proj_tol)$mat),
                error = function(e) fail("projection", e))
  if (is.null(C)) return(na)
  X <- X[rowSums(!is.na(X)) > 0, , drop = FALSE]          # subjects with at least one observation
  V <- tryCatch(.call_getLatentPreds(X = X, type = rep(tp, m), lat.cov.est = C, ridge = ridge,
                                     impute.missing = TRUE), error = function(e) fail("latent prediction", e))
  if (is.null(V)) return(na)
  if (!all(is.finite(V))) {
    warning(sprintf("select_df_bic: variable %s (%s), df = %d: non-finite latent predictions.", j, tp, K), call. = FALSE)
    return(na) }
  n <- nrow(V)
  R <- chol(C)
  logdet <- n * 2 * sum(log(diag(R)))
  quad <- sum(backsolve(R, t(V), transpose = TRUE)^2)
  c(bic = logdet + quad + K * (K + 1) / 2 * log(n), logdet = logdet, quad = quad, n = n)
}
