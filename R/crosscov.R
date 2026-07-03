#' Cross-covariance surface between two mixed-type functional variables
#'
#' Estimates the latent cross-covariance surface \eqn{C_{jk}(s,t)} between two
#' functional variables of possibly different types (continuous, truncated,
#' ordinal, binary), using the appropriate cross-type bridging function of
#' Kendall's tau under the Semiparametric Gaussian Copula model. A tensor-product
#' B-spline is fit to the bridged Kendall's tau surface by nonlinear least
#' squares. Handles regular and irregular sampling through the `min_no_pairs`
#' reliability cutoff (pairs with fewer co-observations are dropped).
#'
#' @param dat_list A list of two numeric matrices (each n by m); the two
#'   functional variables observed on a common grid. `NA` entries are allowed
#'   (irregular / sparse sampling).
#' @param type Character vector of length 2 giving the type of each variable,
#'   each one of `"cont"`, `"trunc"`, `"ord"`, `"bin"`.
#' @param argvals Numeric vector of length m of argument (time) values. Defaults
#'   to an equally spaced grid on \eqn{[0,1]}.
#' @param df Integer; degrees of freedom (number of B-spline basis functions per
#'   margin) for the tensor-product spline.
#' @param min_no_pairs Minimum number of co-observed pairs required to use a
#'   given (s, t) cell; cells with fewer pairs are dropped before fitting.
#' @param weights Logical; if `TRUE` the NLS objective is weighted by the number
#'   of co-observed pairs in each cell (recommended for irregular data).
#'
#' @return A list with elements
#'   \describe{
#'     \item{cov}{the estimated m by m cross-covariance surface}
#'     \item{vcov}{variance-covariance matrix of the fitted spline coefficients}
#'     \item{par}{the fitted `nlsLM` object}
#'   }
#' @references Dey D., Ghosal R., Merikangas K., Zipunnikov V. (2024)
#'   \doi{10.1002/sim.10240}
#' @export
pfpca_crosscov <- function(dat_list, type, argvals = NULL,
                           df = 5, min_no_pairs = 30, weights = TRUE) {

  if (!is.matrix(dat_list[[1]]) || !is.matrix(dat_list[[2]])) {
    stop("The elements of 'dat_list' must be n * m matrices.")
  }
  n <- nrow(dat_list[[1]])
  m <- ncol(dat_list[[1]])
  if (is.null(argvals)) argvals <- seq(0, 1, length = m)
  if (!is.numeric(df) || df <= 0 || df != round(df)) {
    stop("The 'df' parameter must be a positive integer.")
  }

  # sort the pair by type for a canonical bridging branch
  z1 <- dat_list[[order(type)[1]]]
  z2 <- dat_list[[order(type)[2]]]
  type_orig <- type
  type <- sort(type)

  # tensor-product B-spline basis
  cmb <- combn(m, 2)
  knotsT <- seq(min(argvals), max(argvals), l = df - 2)
  norder <- 4
  nbasisT <- length(knotsT) + norder - 2
  dayrngT <- c(min(argvals), max(argvals))
  bbasisT <- create.bspline.basis(dayrngT, nbasisT, norder, knotsT)
  bs.argvals <- eval.basis(argvals, bbasisT)
  formula.spline <- function(j, l) paste(paste0("u", j, l), "*", paste0("Tj", j), "*", paste0("Tl", l))
  formula.spline <- Vectorize(formula.spline, vectorize.args = c("j", "l"))
  spl.f <- paste(as.character(outer(1:df, 1:df, formula.spline)), collapse = " + ")

  if (type[1] == "bin" & type[2] == "ord") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      hatdelta_j <- qnorm(1 - mean(data1[, j], na.rm = TRUE))
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      cats_l <- sort(unique(na.omit(data2[, l])))
      hatdelta_l <- unlist(lapply(cats_l, function(x) qnorm(1 - mean(data2[, l] >= x, na.rm = TRUE))))[-1]
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, hatdelta_j, hatdelta_l, njl = njl)
      } else c(k = NA, Tj, Tl, hatdelta_j, hatdelta_l, njl = njl)
    }
    ncutoff <- length(sort(unique(na.omit(z2[, 1])))) - 1
    l <- vector(mode = "list", length = ncutoff + 2)
    names(l) <- c("t", "dj", paste0("dl", 1:ncutoff))
    bridge_bo_wrapper <- make_function(args = l, body = quote({
      args_fa <- unlist(as.list(environment()))
      delta1 <- args_fa[2]; delta2 <- args_fa[3:(ncutoff + 2)]
      bridge_bo(as.numeric(args_fa[1]), delta1, delta2)
    }))
    bridgeF_bo_v <- Vectorize(bridge_bo_wrapper, vectorize.args = c("t", "dj", paste0("dl", 1:ncutoff)))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_bo_v(ginv(", spl.f, "),",
                               paste(paste("dj", "=", "dj", collapse = ","), ", ",
                                     paste(paste0("dl", 1:ncutoff), "=", paste0("dl", 1:ncutoff), collapse = ",")), ")"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), "dj", paste0("dl", 1:ncutoff))

  } else if (type[1] == "bin" & type[2] == "bin") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      zratioj <- mean(data1[, j] == 0, na.rm = TRUE)
      zratiol <- mean(data2[, l] == 0, na.rm = TRUE)
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, zratioj, zratiol, njl = njl)
      } else c(k = NA, Tj, Tl, zratioj, zratiol, njl = njl)
    }
    bridgeF_bb_v <- Vectorize(bridgeF_bb, vectorize.args = c("r", "zratio1", "zratio2"))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_bb_v(ginv(", spl.f, "), dj, dl)"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), "dj", "dl")

  } else if (type[1] == "bin" & type[2] == "trunc") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      zratioj <- mean(data1[, j] == 0, na.rm = TRUE)
      zratiol <- mean(data2[, l] == 0, na.rm = TRUE)
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, zratioj, zratiol, njl = njl)
      } else c(k = NA, Tj, Tl, zratioj, zratiol, njl = njl)
    }
    bridgeF_bt_v <- Vectorize(bridgeF_bt, vectorize.args = c("r", "zratio1", "zratio2"))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_bt_v(ginv(", spl.f, "), dj, dl)"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), "dj", "dl")

  } else if (type[1] == "bin" & type[2] == "cont") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      zratioj <- mean(data1[, j] == 0, na.rm = TRUE)
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, zratioj, njl = njl)
      } else c(k = NA, Tj, Tl, zratioj, njl = njl)
    }
    bridgeF_bc_v <- Vectorize(bridgeF_bc, vectorize.args = c("r", "zratio1"))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_bc_v(ginv(", spl.f, "), dj)"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), "dj")

  } else if (type[1] == "ord" & type[2] == "ord") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      cats_j <- sort(unique(na.omit(data1[, j])))
      hatdelta_j <- unlist(lapply(cats_j, function(x) qnorm(1 - mean(data1[, j] >= x, na.rm = TRUE))))[-1]
      cats_l <- sort(unique(na.omit(data2[, l])))
      hatdelta_l <- unlist(lapply(cats_l, function(x) qnorm(1 - mean(data2[, l] >= x, na.rm = TRUE))))[-1]
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, hatdelta_j, hatdelta_l, njl = njl)
      } else c(k = NA, Tj, Tl, hatdelta_j, hatdelta_l, njl = njl)
    }
    ncutoff1 <- length(sort(unique(na.omit(z1[, 1])))) - 1
    ncutoff2 <- length(sort(unique(na.omit(z2[, 1])))) - 1
    l <- vector(mode = "list", length = ncutoff1 + ncutoff2 + 1)
    names(l) <- c("t", paste0("dj", 1:ncutoff1), paste0("dl", 1:ncutoff2))
    bridge_oo_fast_wrapper <- make_function(args = l, body = quote({
      args_fa <- unlist(as.list(environment()))
      delta1 <- args_fa[2:(ncutoff1 + 1)]
      delta2 <- args_fa[(ncutoff1 + 2):(ncutoff1 + ncutoff2 + 1)]
      bridge_oo(args_fa[1], delta1, delta2)
    }))
    bridgeF_oo_v <- Vectorize(bridge_oo_fast_wrapper, vectorize.args = c("t", paste0("dj", 1:ncutoff1), paste0("dl", 1:ncutoff2)))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_oo_v(ginv(", spl.f, "),",
                               paste(paste(paste0("dj", 1:ncutoff1), "=", paste0("dj", 1:ncutoff1), collapse = ","),
                                     ", ", paste(paste0("dl", 1:ncutoff2), "=", paste0("dl", 1:ncutoff2), collapse = ",")), ")"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)),
                         paste0("dj", 1:ncutoff1), paste0("dl", 1:ncutoff2))

  } else if (type[1] == "ord" & type[2] == "trunc") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      cats_j <- sort(unique(na.omit(data1[, j])))
      hatdelta_j <- unlist(lapply(cats_j, function(x) qnorm(1 - mean(data1[, j] >= x, na.rm = TRUE))))[-1]
      hatdelta_l <- qnorm(1 - mean(data2[, l], na.rm = TRUE))
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, hatdelta_j, hatdelta_l, njl = njl)
      } else c(k = NA, Tj, Tl, hatdelta_j, hatdelta_l, njl = njl)
    }
    ncutoff <- length(sort(unique(na.omit(z1[, 1])))) - 1
    l <- vector(mode = "list", length = ncutoff + 2)
    names(l) <- c("t", paste0("dj", 1:ncutoff), "dl")
    bridge_ot_wrapper <- make_function(args = l, body = quote({
      args_fa <- unlist(as.list(environment()))
      delta1 <- args_fa[2:(ncutoff + 1)]; delta2 <- args_fa[(ncutoff + 2)]
      bridge_ot(as.numeric(args_fa[1]), delta1, delta2)
    }))
    bridgeF_ot_v <- Vectorize(bridge_ot_wrapper, vectorize.args = c("t", paste0("dj", 1:ncutoff), "dl"))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_ot_v(ginv(", spl.f, "),",
                               paste(paste(paste0("dj", 1:ncutoff), "=", paste0("dj", 1:ncutoff), collapse = ","),
                                     ", ", paste("dl", "=", "dl", collapse = ",")), ")"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), paste0("dj", 1:ncutoff), "dl")

  } else if (type[1] == "cont" & type[2] == "ord") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      cats_l <- sort(unique(na.omit(data2[, l])))
      hatdelta_l <- unlist(lapply(cats_l, function(x) qnorm(1 - mean(data2[, l] >= x, na.rm = TRUE))))[-1]
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, hatdelta_l, njl = njl)
      } else c(k = NA, Tj, Tl, hatdelta_l, njl = njl)
    }
    ncutoff <- length(sort(unique(na.omit(z2[, 1])))) - 1
    l <- vector(mode = "list", length = ncutoff + 1)
    names(l) <- c("t", paste0("dl", 1:ncutoff))
    bridge_co_wrapper <- make_function(args = l, body = quote({
      args_fa <- unlist(as.list(environment()))
      delta <- args_fa[2:(ncutoff + 1)]
      bridge_co(as.numeric(args_fa[1]), delta)
    }))
    bridgeF_co_v <- Vectorize(bridge_co_wrapper, vectorize.args = c("t", paste0("dl", 1:ncutoff)))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_co_v(ginv(", spl.f, "),",
                               paste(paste(paste0("dl", 1:ncutoff), "=", paste0("dl", 1:ncutoff), collapse = ",")), ")"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), paste0("dl", 1:ncutoff))

  } else if (type[1] == "cont" & type[2] == "trunc") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      zratiol <- mean(data2[, l] == 0, na.rm = TRUE)
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, zratiol, njl = njl)
      } else c(k = NA, Tj, Tl, zratiol, njl = njl)
    }
    bridgeF_ct_v <- Vectorize(bridgeF_ct, vectorize.args = c("r", "zratio2"))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_ct_v(ginv(", spl.f, "), dl)"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), "dl")

  } else if (type[1] == "trunc" & type[2] == "trunc") {
    fjl.df <- function(data1, data2, j, l) {
      Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
      zratioj <- mean(data1[, j] == 0, na.rm = TRUE)
      zratiol <- mean(data2[, l] == 0, na.rm = TRUE)
      njl <- sum(!is.na(data1[, j] * data2[, l]))
      if (njl > min_no_pairs) {
        tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
        c(k = tjl, Tj, Tl, zratioj, zratiol, njl = njl)
      } else c(k = NA, Tj, Tl, zratioj, zratiol, njl = njl)
    }
    bridgeF_tt_v <- Vectorize(bridgeF_tt, vectorize.args = c("r", "zratio1", "zratio2"))
    eunsc <- as.formula(paste0("k ~ ", "bridgeF_tt_v(ginv(", spl.f, "), dj, dl)"))
    obj_df_colnames <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), "dj", "dl")
  }

  # continuous-continuous warm start (always run)
  fjl.df.cont <- function(data1, data2, j, l) {
    Tj <- bs.argvals[j, ]; Tl <- bs.argvals[l, ]
    njl <- sum(!is.na(data1[, j] * data2[, l]))
    if (njl > min_no_pairs) {
      tjl <- Kendall_mixed(cbind(data1[, j], data2[, l]))[1, 2]
      c(k = tjl, Tj, Tl, njl = njl)
    } else c(k = NA, Tj, Tl, njl = njl)
  }
  bridgeF_cc_v <- Vectorize(bridgeF_cc, vectorize.args = c("r"))
  eunsc_cont <- as.formula(paste0("k ~ ", "bridgeF_cc_v(ginv(", spl.f, "))"))
  obj_df_colnames_cont <- c("k", paste0("Tj", 1:ncol(bs.argvals)), paste0("Tl", 1:ncol(bs.argvals)), "njl")
  obj_df_colnames <- c(obj_df_colnames, "njl")

  # include diagonal (s == t) and both orderings of each off-diagonal pair
  c2 <- cbind(matrix(rep(1:m, each = 2), ncol = m), cmb, apply(combn(m, 2), 2, function(x) c(max(x), min(x))))
  obj_df_cont <- data.frame(t(sapply(1:ncol(c2), function(x) {
    res <- fjl.df.cont(data1 = z1, data2 = z2, c2[1, x], c2[2, x])
    if (is.na(res[1])) res <- rep(NA, length(obj_df_colnames_cont))
    res
  })))
  names(obj_df_cont) <- obj_df_colnames_cont

  init <- runif(df^2, -5, 5)
  formula.start <- function(j, l) paste(paste0("u", j, l))
  formula.start <- Vectorize(formula.start, vectorize.args = c("j", "l"))
  names(init) <- outer(1:df, 1:df, formula.start)
  obj_df_cont <- obj_df_cont[complete.cases(obj_df_cont), ]

  if (weights) {
    ns0 <- nlsLM(eunsc_cont, data = obj_df_cont, start = init,
                 control = nls.lm.control(ftol = 0.001), weights = obj_df_cont$njl)
  } else {
    ns0 <- nlsLM(eunsc_cont, data = obj_df_cont, start = init,
                 control = nls.lm.control(ftol = 0.001))
  }
  ns1 <- ns0

  if (!(type[1] == "cont" & type[2] == "cont")) {
    obj_df <- data.frame(t(sapply(1:ncol(c2), function(x) {
      res <- fjl.df(data1 = z1, data2 = z2, c2[1, x], c2[2, x])
      if (is.na(res[1])) res <- rep(NA, length(obj_df_colnames))
      res
    })))
    names(obj_df) <- obj_df_colnames
    init <- coef(summary(ns0))[, 1]
    # drop (s,t) cells with too few pairs (irregular-data fix: assignment kept)
    obj_df <- obj_df[complete.cases(obj_df), ]
    if (weights) {
      ns1 <- nlsLM(eunsc, data = obj_df, start = init,
                   control = nls.lm.control(ftol = 0.001), weights = obj_df$njl)
    } else {
      ns1 <- nlsLM(eunsc, data = obj_df, start = init,
                   control = nls.lm.control(ftol = 0.001))
    }
  }

  uhat.nls <- matrix(coef(summary(ns1))[, 1], ncol = df, nrow = df)
  Chat1 <- outer(argvals, argvals, Chat, u = uhat.nls, bs = bbasisT)
  Chat2 <- ginv(Chat1)
  if (type_orig[1] != type[1]) Chat2 <- t(Chat2)
  list(cov = Chat2, vcov = vcov(ns1), par = ns1)
}
