#' Construct spline eigenfunction bases for the simulated factors
#'
#' Builds one orthonormal spline system per latent factor using the
#' \pkg{Splinets} package. Each factor uses a different spline order, number of
#' knots and periodicity so that the simulated factors have distinct temporal
#' shapes.
#'
#' @param Q Number of latent factors. Must be at most 4.
#' @param L Number of FPCA components per factor. Must be at most 4.
#'
#' @return A named list with \code{vec_K} (number of knots per factor) and
#'         \code{splines} (list of \pkg{Splinets} objects, one per factor).
#'
#' @export
#'
construct_splines <- function(Q = 4, L = 3) {

  if (!requireNamespace("Splinets", quietly = TRUE)) {
    stop("Package 'Splinets' is required to generate simulated data. ",
         "Install it with install.packages('Splinets').")
  }

  splines <- vector("list", length = Q)
  vec_K <- NULL
  stopifnot(L <= 4)
  stopifnot(Q <= 4)

  n_int_knots <- L + 2
  for (q in 1:Q) {
    if (q == 1) {
      K <- n_int_knots + 3
      k <- 3
      bool_periodic <- FALSE
    } else if (q == 2) {
      K <- n_int_knots + 2
      k <- 2
      bool_periodic <- FALSE
    } else if (q == 3) {
      K <- n_int_knots + 4
      k <- 3
      bool_periodic <- TRUE
    } else if (q == 4) {
      K <- n_int_knots + 3
      k <- 2
      bool_periodic <- TRUE
    } else {
      stop("Not implemented for Q > 4.")
    }
    knots <- seq(0, 1, length.out = K)
    # 'degree' (not 'order') since Splinets renamed the argument, see
    # https://github.com/ranibasna/R-Splinets/commit/a8f79565257da5d292acbd4fc03af62ecd3f899e
    so <- Splinets::splinet(knots, degree = k, periodic = bool_periodic, norm = TRUE)
    splines[[q]] <- so
    vec_K <- c(vec_K, K)
  }
  create_named_list(vec_K, splines)
}


#' Simulate multivariate longitudinal data from the bayesSYNC model
#'
#' Generates high-dimensional curves following the bayesSYNC generative model:
#' a shared set of latent factors, each with subject-specific FPCA scores and
#' spline eigenfunctions, combined through variable-specific sparse loadings and
#' added to variable-specific mean functions, with Gaussian measurement noise.
#' The output is in the exact format expected by \code{\link{bayesSYNC}} and the
#' simulated ground truth can be aligned to the estimates with
#' \code{\link{match_factor_and_sign}}.
#'
#' @param N Number of subjects.
#' @param n Integer vector of length \code{N} giving the number of observation
#'          times per subject.
#' @param p Number of variables (e.g., biomarkers).
#' @param Q Number of latent factors. Must be at most 4 (see
#'          \code{\link{construct_splines}}).
#' @param L Number of FPCA components per factor. Must be at most 4.
#' @param mu_func Mean function with signature \code{mu_func(t, j, p)} returning
#'                the mean of variable \code{j} at times \code{t}.
#' @param time_obs Optional list of length \code{N} of observation-time vectors.
#'                 If \code{NULL}, times are drawn uniformly on [0, 1].
#' @param n_g Size of the dense grid on which the true mean and eigenfunctions
#'            are evaluated.
#' @param a_om,b_om Shape parameters of the Beta prior on the factor inclusion
#'                  probabilities \code{omega}. The default \code{b_om = p} yields
#'                  sparse loadings; lower \code{b_om} to make more variables load
#'                  on each factor.
#' @param vec_sd_zeta Vector of length \code{L} of score standard deviations per
#'                    component. Defaults to \code{1/(1:L)}.
#' @param vec_sd_eps Measurement-noise standard deviation, a scalar or a vector
#'                   of length \code{p}. Defaults to 1 for every variable.
#' @param seed Optional seed for reproducibility.
#'
#' @return A named list with, among others: \code{time_obs}, \code{Y} (the data
#'         in \code{\link{bayesSYNC}} input format), the true loadings \code{B},
#'         inclusion indicators \code{gamma}, scores \code{Zeta} (list of
#'         \code{Q} matrices of size \code{N x L}), eigenfunctions \code{Phi} and
#'         their dense-grid version \code{Phi_g} (list of \code{Q} matrices of
#'         size \code{n_g x L}), and the mean functions \code{mu} and \code{mu_g}.
#'
#' @seealso \code{\link{bayesSYNC}}, \code{\link{match_factor_and_sign}},
#'          \code{\link{construct_splines}}
#'
#' @importFrom stats rbeta rbinom rnorm runif
#' @export
#'
generate_bayesSYNC_data <- function(N, n, p, Q, L, mu_func,
                                    time_obs = NULL,
                                    n_g = 1000,
                                    a_om = 1, b_om = p,
                                    vec_sd_zeta = NULL,
                                    vec_sd_eps = NULL,
                                    seed = NULL) {

  if (!requireNamespace("Splinets", quietly = TRUE)) {
    stop("Package 'Splinets' is required to generate simulated data. ",
         "Install it with install.packages('Splinets').")
  }

  check_structure(seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  check_structure(n, "vector", "numeric", N)

  if (is.null(time_obs)) {
    time_obs <- lapply(n, runif)
    time_obs <- lapply(time_obs, sort)
  }

  if (is.null(vec_sd_zeta)) {
    vec_sd_zeta <- 1 / (1:L)
  } else {
    stopifnot(length(vec_sd_zeta) == L)
  }

  if (is.null(vec_sd_eps)) {
    vec_sd_eps <- rep(1, p)
  } else if (length(vec_sd_eps) == 1) {
    vec_sd_eps <- rep(vec_sd_eps, p)
  } else {
    stopifnot(length(vec_sd_eps) == p)
  }

  omega <- rbeta(Q, shape1 = a_om, shape2 = b_om)
  gamma <- matrix(rbinom(p * Q, size = 1, prob = omega), nrow = p, ncol = Q)

  # whether all factors should be "active",
  # i.e., represent the contribution of at least one variable j = 1, ..., p
  bool_active <- TRUE
  if (bool_active) {
    if (any(colSums(gamma) == 0)) {
      for (q in which(colSums(gamma) == 0)) {
        gamma[sample(1:p, 1), q] <- 1
      }
    }
  }

  B_normal <- matrix(rnorm(p * Q, mean = 0, sd = 1), nrow = p, ncol = Q)

  B <- B_normal * gamma

  Zeta <- lapply(1:Q, function(q) { # list of length Q, containing matrices of size N x L
    matrix(rnorm(N * L, 0, rep(vec_sd_zeta, each = N)), nrow = N, ncol = L)
  })

  mu <- lapply(1:N, function(i) {
    lapply(1:p, function(j) {
      mu_func(time_obs[[i]], j = j, p = p)
    })
  })

  res_splines <- construct_splines(Q = Q, L = L)
  splines <- res_splines$splines
  vec_K <- res_splines$vec_K

  list_ind_spline_functions <- lapply(1:Q, function(q) {
    sample(1:(ncol(Splinets::evspline(splines[[1]]$os)) - 1), L) # remove 1 as the first column is the vector of timepoints
  })

  Phi <- lapply(1:Q, function(q) {
    lapply(1:N, function(i) {
      eval_qi <- Splinets::evspline(splines[[q]]$os, x = time_obs[[i]])[, -1, drop = FALSE]
      eval_qi[, list_ind_spline_functions[[q]], drop = FALSE]
    })
  })

  Y <- lapply(1:N, function(i) {
    Y_i <- lapply(1:p, function(j) {

      tmp <- sapply(1:Q, function(q) B[j, q] * Phi[[q]][[i]] %*% Zeta[[q]][i, ])
      if (is.vector(tmp)) { # to deal with cases with a single observation
        tmp <- t(as.matrix(tmp))
      }
      mu[[i]][[j]] +
        rowSums(tmp) +
        rnorm(length(time_obs[[i]]), mean = 0, sd = vec_sd_eps[j])
    })
    names(Y_i) <- paste0("variable_", 1:p)
    Y_i
  })

  names(Y) <- paste0("subject_", 1:N)

  time_g <- seq(0, 1, length.out = n_g)

  mu_g <- lapply(1:p, function(j) mu_func(time_g, j = j, p = p))

  Phi_g <- lapply(1:Q, function(q) {
    eval_q <- Splinets::evspline(splines[[q]]$os, x = time_g)[, -1, drop = FALSE]
    eval_q[, list_ind_spline_functions[[q]], drop = FALSE]
  })

  create_named_list(time_obs, Y, B, omega, gamma, mu, Phi, Zeta,
                    vec_sd_zeta, vec_sd_eps,
                    time_g, mu_g, Phi_g, vec_K)
}
