# Internal helpers for M2FPCA.
#
# The bridging functions and small utilities below are vendored (copied) from
# the SGCTools package so that M2FPCA is self-contained and does not rely on
# SGCTools' unexported (:::) internals. They are intentionally NOT exported.
# The four mixed-type bridges that live only inside mixedCCA (bridgeF_tt/bt/bc/ct)
# are wrapped via mixedCCA::: at the bottom of this file.
#
# Keep these in sync with SGCTools if its bridging functions ever change.

# ---- function factory (vendored from SGCTools) ----------------------------
#' @noRd
make_function <- function(args, body, env = parent.frame()) {
  args <- as.pairlist(args)
  eval(call("function", args, body), env)
}

# ---- correlation-scale link (vendored from SGCTools) -----------------------
# Inverse of g(x) = log((1+x)/(1-x)); maps spline predictor onto (-1, 1).
#' @noRd
ginv <- function(x) {
  (exp(x) - 1) / (exp(x) + 1)
}

# Derivative of ginv() w.r.t. its argument (used for delta-method variances).
#' @noRd
ginv_deriv <- function(x) {
  2 * exp(x) / (exp(x) + 1)^2
}

# ---- tensor-product spline covariance evaluation (vendored from SGCTools) ---
#' @noRd
fb.sp <- function(u, Tj, Tl) {
  t(Tj) %*% u %*% Tl
}

#' @noRd
Chat <- function(u, s, t, bs) {
  fb.sp(u, as.numeric(eval.basis(s, bs)), as.numeric(eval.basis(t, bs)))
}
Chat <- Vectorize(Chat, vectorize.args = c("s", "t"))

# ---- within-type bridging functions (vendored from SGCTools) ----------------
#' @noRd
bridgeF_cc <- function(r) {
  2 / pi * asin(r)
}

#' @noRd
bridgeF_bb <- function(r, zratio1, zratio2) {
  de1 <- qnorm(zratio1)
  de2 <- qnorm(zratio2)
  if (r[1] == 1)  r <- r - 1e-6
  if (r[1] == -1) r <- r + 1e-6
  as.numeric(2 * (fMultivar::pnorm2d(de1, de2, rho = r) - zratio1 * zratio2))
}

# continuous - ordinal bridging function
#' @noRd
bridge_co <- function(t, delta, ...) {
  M <- 12
  delta <- c(delta, M)
  S <- cbind(c(1, 0, t / sqrt(2)), c(0, 1, -t / sqrt(2)), c(t / sqrt(2), -t / sqrt(2), 1))
  l <- length(delta)
  term1 <- sum(sapply(1:(l - 1), function(x) {
    4 * pmvnorm(upper = c(delta[x], delta[x + 1], 0), sigma = S, algorithm = Miwa()) -
      2 * pnorm(delta[x]) * pnorm(delta[x + 1])
  }))
  term1
}

# ordinal - ordinal bridging function
#' @noRd
bridge_oo <- function(t, delta1, delta2) {
  M <- 12
  l1 <- length(delta1)
  l2 <- length(delta2)
  if (t[1] == 1)  t <- t - 1e-6
  if (t[1] == -1) t <- t + 1e-6
  setting1 <- expand.grid(c(1:l1), c(1:l2))
  delta_t1 <- c(-M, delta1, M)
  delta_t2 <- c(-M, delta2, M)
  term1 <- sum(sapply(1:nrow(setting1), function(i) {
    x <- as.numeric(setting1[i, ])
    pnorm2d(delta_t1[x[1] + 1], delta_t2[x[2] + 1], rho = t) *
      (pnorm2d(delta_t1[x[1] + 2], delta_t2[x[2] + 2], rho = t) -
         pnorm2d(delta_t1[x[1] + 2], delta_t2[x[2]], rho = t))
  }))
  term2 <- sum(sapply(1:l1, function(i) {
    x <- as.numeric(setting1[i, ])
    pnorm(delta_t1[x[1] + 1]) * (pnorm2d(delta_t1[x[1] + 2], delta_t2[l2 + 1], rho = t))
  }))
  c(2 * (term1 - term2))
}

# binary - ordinal bridging function
#' @noRd
bridge_bo <- function(t, delta1, delta2) {
  M <- 12
  if (t[1] == 1)  t <- t - 1e-6
  if (t[1] == -1) t <- t + 1e-6
  l1 <- length(delta1)
  delta_t1 <- c(delta2, M)
  Splus  <- cbind(c(1, t), c(t, 1))
  Sminus <- cbind(c(1, -t), c(-t, 1))
  term1 <- sum(sapply(1:(length(delta_t1) - 1), function(x) {
    pmvnorm(lower = c(delta_t1[x], -M), upper = c(delta_t1[x + 1], -delta1), sigma = Sminus, algorithm = Miwa()) *
      pmvnorm(upper = c(delta_t1[x], delta1), sigma = Splus, algorithm = Miwa())
  }))
  term3 <- sum(sapply(1:(length(delta_t1) - 1), function(x) {
    pmvnorm(lower = c(delta_t1[x], -M), upper = c(delta_t1[x + 1], delta1), sigma = Splus, algorithm = Miwa()) *
      pmvnorm(upper = c(delta_t1[x], -delta1), sigma = Sminus, algorithm = Miwa())
  }))
  c(2 * (term1 - term3))
}

# ordinal - truncated bridging function
#' @noRd
bridge_ot <- function(t, delta1, delta2) {
  M <- 12
  delta_t1 <- c(delta1, M)
  delta_t2 <- c(-M, delta1)
  Splus  <- cbind(c(1, t), c(t, 1))
  Sminus <- cbind(c(1, -t), c(-t, 1))
  S4plus  <- cbind(c(1, 0, 0, -t / sqrt(2)), c(0, 1, -t, t / sqrt(2)), c(0, -t, 1, -1 / sqrt(2)), c(-t / sqrt(2), t / sqrt(2), -1 / sqrt(2), 1))
  S4minus <- cbind(c(1, 0, 0, -t / sqrt(2)), c(0, 1, t, -t / sqrt(2)), c(0, t, 1, -1 / sqrt(2)), c(-t / sqrt(2), -t / sqrt(2), -1 / sqrt(2), 1))
  term1 <- sum(sapply(1:(length(delta_t1) - 1), function(x) {
    pmvnorm(lower = c(delta_t1[x], -M), upper = c(delta_t1[x + 1], -delta2), sigma = Sminus, algorithm = Miwa()) *
      pmvnorm(upper = c(delta_t1[x], delta2), sigma = Splus, algorithm = Miwa()) +
      pmvnorm(lower = c(delta_t1[x], -M, -M, -M), upper = c(delta_t1[x + 1], delta_t1[x], -delta2, 0), sigma = S4plus, algorithm = Miwa())
  }))
  term2 <- sum(sapply(1:(length(delta_t2) - 1), function(x) {
    pmvnorm(lower = c(delta_t2[x], -M), upper = c(delta_t2[x + 1], -delta2), sigma = Sminus, algorithm = Miwa()) *
      pmvnorm(upper = c(-delta_t2[x + 1], delta2), sigma = Sminus, algorithm = Miwa()) +
      pmvnorm(lower = c(delta_t2[x], -M, -M, -M), upper = c(delta_t2[x + 1], -delta_t2[x + 1], -delta2, 0), sigma = S4minus, algorithm = Miwa())
  }))
  2 * (term1 - term2)
}

# ---- mixed-type bridges that live only in mixedCCA -------------------------
# These four are not in SGCTools; we wrap mixedCCA's unexported versions.
# (binary-truncated, binary-continuous, continuous-truncated, truncated-truncated)
#' @noRd
bridgeF_bt <- function(r, zratio1, zratio2) mixedCCA:::bridgeF_bt(r, zratio1, zratio2)
#' @noRd
bridgeF_bc <- function(r, zratio1)          mixedCCA:::bridgeF_bc(r, zratio1)
#' @noRd
bridgeF_ct <- function(r, zratio2)          mixedCCA:::bridgeF_ct(r = r, zratio1 = NULL, zratio2 = zratio2)
#' @noRd
bridgeF_tt <- function(r, zratio1, zratio2) mixedCCA:::bridgeF_tt(r, zratio1, zratio2)

# ---- asymptotic-variance helper matrices (vendored from mpfpca_functions) ---
#' @noRd
up.elim.diag <- function(p) {
  E_u <- matrix(0, nrow = p * (p + 1) / 2, ncol = p^2)
  E_i <- matrix(1:p^2, ncol = p)
  e <- E_i[upper.tri(E_i, diag = TRUE)]
  for (i in 1:nrow(E_u)) E_u[i, ][e[i]] <- 1
  t(E_u)
}

#' @noRd
beta.elim2 <- function(p, v) {
  p1 <- length(v)
  E_b <- matrix(0, nrow = (p - p1) * p, ncol = p^2)
  E_i <- matrix(1:p^2, ncol = p)
  e <- c(E_i[-v, v], as.numeric(E_i[-v, -v]))
  for (i in 1:nrow(E_b)) E_b[i, ][e[i]] <- 1
  E_b
}

# Renamed from deriv.b2 -> deriv_b2 so roxygen2 does not mistake it for an
# S3 method of the stats::deriv generic.
#' @noRd
deriv_b2 <- function(A, dep.v) {
  p <- ncol(A)
  Rinv <- as.matrix(solve(A[-dep.v, -dep.v]))
  p1 <- length(dep.v)
  deriv.beta1 <- kronecker(diag(p1), Rinv)
  deriv.beta2 <- -(kronecker(Rinv, Rinv)) %*% (kronecker(diag(p - p1), A[-dep.v, dep.v]))
  deriv.beta <- rbind(deriv.beta1, deriv.beta2)
  t(deriv.beta)
}

# ---- numerical integration + eigenanalysis (replace fdapace / fgm) ----------
# Trapezoidal rule (replaces fdapace::trapzRcpp).
#' @noRd
trapz <- function(x, y) {
  n <- length(x)
  sum((y[-1] + y[-n]) / 2 * diff(x))
}

# Eigenanalysis of a covariance surface on grid `regGrid`.
#
# Faithful reimplementation of fgm:::.GetEigenAnalysisResults
# (https://github.com/javzapata/fgm, R/pfpca.R), which itself follows fdapace's
# GetEigenAnalysisResults. Matches fgm exactly: gridSize = regGrid[2]-regGrid[1];
# keep eigenvalues >= 0; truncate to maxK BEFORE the FVE threshold; normalise
# eigenfunctions by the trapezoidal integral x/sqrt(int x^2); sign convention
# sum(x * muWork) >= 0 with muWork = 1:m; lambda = gridSize * eigenvalues.
# `cumFVE` is returned on the 0-100 (percent) scale, as in fgm.
#' @noRd
eigen_analysis <- function(smoothCov, regGrid, FVEthreshold = 0.99, maxK = NULL,
                           muWork = NULL) {
  smoothCov <- as.matrix(smoothCov)
  gridSize <- regGrid[2] - regGrid[1]

  eig <- eigen(smoothCov, symmetric = TRUE)
  positiveInd <- eig$values >= 0
  if (sum(positiveInd) == 0) {
    stop("All eigenvalues are negative. The covariance estimate is incorrect.")
  }
  d <- eig$values[positiveInd]
  eigenV <- eig$vectors[, positiveInd, drop = FALSE]

  if (is.null(maxK)) maxK <- length(d)
  if (maxK < length(d)) {
    d <- d[1:maxK]
    eigenV <- eigenV[, 1:maxK, drop = FALSE]
  }

  # cumulative FVE (percent), threshold applied after maxK truncation
  FVE <- cumsum(d) / sum(d) * 100
  no_opt <- min(which(FVE >= FVEthreshold * 100))

  if (is.null(muWork)) muWork <- 1:nrow(eigenV)
  phi <- apply(eigenV, 2, function(x) {
    x <- x / sqrt(trapz(regGrid, x^2))
    if (0 <= sum(x * muWork)) x else -x
  })
  lambda <- gridSize * d

  list(lambda = lambda[1:no_opt],
       phi    = phi[, 1:no_opt, drop = FALSE],
       cumFVE = FVE,
       kChoosen = no_opt)
}
