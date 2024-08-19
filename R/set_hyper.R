#' Gather model hyperparameters.
#'
#' This function is used to construct an object with hyperparameter values for
#' the FPCA or mFPCA model in a format that can then be supplied to
#' \code{\link{bayesSYNC}}.
#'
#' The \code{\link{bayesSYNC}} function can
#' also be used with default hyperparameter values (i.e., without using
#' \code{\link{set_hyper}}) by setting their argument \code{list_hyper} to
#' \code{NULL}.
#'
#' @param sigma_beta Vector of size 1 or 2 for the standard deviation of the
#'                   spline coefficients for the mean function.
#' @param A Positive real number for the top-level hyperparameter of the
#'          iterated inverse chi-square distributions.
#' @param c_0 to add
#' @param d_0 to add
#'
#' @return An object containing the hyperparameter settings to be supplied to
#'         \code{\link{bayesSYNC}}.
#'
#' @export
#'
set_hyper <- function(sigma_beta = 1e5, A = 1e5, c_0 = 1, d_0 = 1) {

  check_structure(sigma_beta, "vector", "numeric", c(1, 2))
  check_positive(sigma_beta)

  if (length(sigma_beta) == 1) sigma_beta <- rep(sigma_beta, 2)
  Sigma_beta <- sigma_beta^2*diag(2)

  check_structure(A, "vector", "numeric", 1)
  check_positive(A)

  check_structure(c_0, "vector", "numeric", 1)
  check_structure(d_0,"vector", "numeric", 1)

  sigma_zeta <- 1
  mu_beta <- rep(0, 2)

  list_hyper <- create_named_list(sigma_zeta, mu_beta, Sigma_beta, A, c_0, d_0)

  class(list_hyper) <- "hyper"

  list_hyper

}
