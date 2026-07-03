#' M2FPCA: Multivariate FPCA for Mixed-Type Functional Data
#'
#' Joint covariance estimation and multivariate functional principal component
#' analysis for multivariate mixed-type functional data (continuous, truncated,
#' ordinal, binary) via the latent Semiparametric Gaussian Copula process.
#'
#' The package estimates the marginal covariance surface of each functional
#' variable (through [SGCTools::fpca.sgc.lat()]) and the cross-covariance
#' surface between every pair of variables ([pfpca_crosscov()]), assembles them
#' into a partially separable multivariate Karhunen-Loeve representation
#' ([mpfpca.dir()]), reconstructs the full block covariance
#' ([cov.from.mpfpca.dir()]), and computes multivariate principal-component
#' scores ([mpfpca_scores()]). Regular and irregular (sparse / asynchronous)
#' sampling are both supported via the `min_no_pairs` reliability cutoff.
#'
#' @references
#' Dey D., Ghosal R., Merikangas K., Zipunnikov V. (2024) "Functional Principal
#' Component Analysis for Continuous Non-Gaussian, Truncated, and Discrete
#' Functional Data." Statistics in Medicine, 43, 5431-5445.
#' \doi{10.1002/sim.10240}
#'
#' @keywords internal
#' @importFrom SGCTools fpca.sgc.lat getLatentPreds Kendall_mixed
#' @importFrom fda create.bspline.basis eval.basis
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom Matrix nearPD
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom fMultivar pnorm2d
#' @importFrom stats as.formula coef complete.cases cov optimize qnorm pnorm runif vcov na.omit
#' @importFrom utils combn
"_PACKAGE"
