# Builds the inputs of compute_factor_pve() from synthetic variational
# parameters passed through the post-processing (orthonormalise), without
# fitting a model. c_scale rescales the loadings by c_scale and the scores by
# 1 / c_scale before the post-processing.
make_synthetic_inputs <- function(c_scale = 1, sigsq_eps = NULL,
                                  N = 30, p = 8, Q = 3, L = 2, K = 8,
                                  n_g = 200, seed = 1) {
  set.seed(seed)
  time_obs <- lapply(1:N, function(i) sort(runif(10)))
  grid_obj <- get_grid_objects(time_obs, K, n_g = n_g, format_univ = TRUE)

  mu_q_nu_mu <- lapply(1:p, function(j) rnorm(K + 2))
  mu_q_nu_phi <- lapply(1:Q, function(q) matrix(rnorm((K + 2) * L), K + 2, L))
  mu_q_zeta <- lapply(1:Q, function(q) matrix(rnorm(N * L), N, L))
  Sigma_q_zeta <- lapply(1:Q, function(q) lapply(1:N, function(i) diag(0.1, L)))
  mu_q_normal_b <- matrix(rnorm(p * Q), p, Q)
  Sigma_q_normal_b <- matrix(runif(p * Q, 0.01, 0.1), p, Q)
  mu_q_gamma <- matrix(runif(p * Q), p, Q)
  if (is.null(sigsq_eps)) sigsq_eps <- runif(p, 0.1, 1)

  # rescaling, after all random draws so that they do not depend on c_scale
  mu_q_normal_b <- c_scale * mu_q_normal_b
  Sigma_q_normal_b <- c_scale^2 * Sigma_q_normal_b
  mu_q_zeta <- lapply(mu_q_zeta, function(M) M / c_scale)
  Sigma_q_zeta <- lapply(Sigma_q_zeta, function(ll) lapply(ll, function(S) S / c_scale^2))
  mu_q_b <- mu_q_gamma * mu_q_normal_b

  res_orth <- orthonormalise(N, p, Q, L, grid_obj$time_g, grid_obj$C_g,
                             mu_q_nu_mu, mu_q_b, mu_q_zeta, mu_q_nu_phi,
                             Sigma_q_zeta, sigsq_eps)

  ppi <- mu_q_gamma
  colnames(ppi) <- paste0("Factor_", 1:Q)

  list(B_slab_mean = mu_q_normal_b,
       B_slab_var = Sigma_q_normal_b,
       ppi = ppi,
       factor_ppi = 1 - exp(colSums(log1p(-ppi))),
       list_eigenvalues = lapply(res_orth$list_Zeta_hat, function(Z) apply(Z, 2, var)),
       sigsq_eps = sigsq_eps,
       time_g = grid_obj$time_g,
       list_Phi_hat = res_orth$list_list_Phi_hat)
}

pve_from_inputs <- function(inputs, ...) {
  compute_factor_pve(inputs$B_slab_mean, inputs$B_slab_var, inputs$ppi,
                     inputs$factor_ppi, inputs$list_eigenvalues,
                     inputs$sigsq_eps, inputs$time_g, ...)
}

# Largest absolute deviation of the trapezoidal Gram matrix of the
# eigenfunctions from the identity.
gram_deviation <- function(Phi, time_g) {
  L <- ncol(Phi)
  G <- matrix(NA, L, L)
  for (l in 1:L) {
    for (l_tilde in 1:L) {
      G[l, l_tilde] <- trapint(time_g, Phi[, l] * Phi[, l_tilde])
    }
  }
  max(abs(G - diag(L)))
}

inputs <- make_synthetic_inputs()
option_grid <- expand.grid(factors = c("all", "mpm"),
                           loadings = c("posterior", "mpm"),
                           stringsAsFactors = FALSE)

test_that("structured PVE sums to one over the included factors", {
  for (k in seq_len(nrow(option_grid))) {
    res <- pve_from_inputs(inputs, factors = option_grid$factors[k],
                           loadings = option_grid$loadings[k])
    expect_equal(sum(res$structured, na.rm = TRUE), 1)
  }
})

test_that("total PVE is below structured PVE, with equality only without noise", {
  for (k in seq_len(nrow(option_grid))) {
    res <- pve_from_inputs(inputs, factors = option_grid$factors[k],
                           loadings = option_grid$loadings[k])
    inc <- res$factors_included
    expect_true(all(res$total[inc] <= res$structured[inc]))
    # strict for factors with non-zero structured variance (with
    # loadings = "mpm", a factor with no selected loading has both equal to 0)
    inc_nonzero <- inc[res$structured[inc] > 0]
    expect_true(all(res$total[inc_nonzero] < res$structured[inc_nonzero]))
  }

  inputs_no_noise <- make_synthetic_inputs(sigsq_eps = rep(0, 8))
  res <- pve_from_inputs(inputs_no_noise, factors = "all")
  expect_equal(res$total, res$structured)
})

test_that("renormalising total PVE gives structured PVE", {
  for (k in seq_len(nrow(option_grid))) {
    res <- pve_from_inputs(inputs, factors = option_grid$factors[k],
                           loadings = option_grid$loadings[k])
    expect_equal(res$total / sum(res$total, na.rm = TRUE), res$structured)
  }
})

test_that("PVE is invariant to rescaling loadings by c and scores by 1 / c", {
  inputs_scaled <- make_synthetic_inputs(c_scale = 3.7)
  for (k in seq_len(nrow(option_grid))) {
    res <- pve_from_inputs(inputs, factors = option_grid$factors[k],
                           loadings = option_grid$loadings[k])
    res_scaled <- pve_from_inputs(inputs_scaled, factors = option_grid$factors[k],
                                  loadings = option_grid$loadings[k])
    expect_equal(res_scaled$total, res$total)
    expect_equal(res_scaled$structured, res$structured)
  }
})

test_that("eigenfunctions are orthonormal under the trapezoidal rule", {
  # The post-processing makes the eigenfunctions orthogonal in the Euclidean
  # inner product on the grid, so on an evenly spaced grid the trapezoidal
  # cross-products are only zero up to a term of order max(diff(time_g)).
  for (Phi in inputs$list_Phi_hat) {
    expect_lt(gram_deviation(Phi, inputs$time_g), 10 * max(diff(inputs$time_g)))
  }
})

test_that("excluded factors get NA and options are recorded", {
  res <- pve_from_inputs(inputs, factors = c(3, 1), loadings = "mpm")
  expect_equal(unname(res$factors_included), c(1L, 3L))
  expect_true(is.na(res$total[2]) && is.na(res$structured[2]))
  expect_equal(length(res$total), ncol(inputs$ppi))
  expect_setequal(res$order, c("Factor_1", "Factor_3"))
  expect_equal(res$options, list(factors = c(3, 1), loadings = "mpm"))
})

test_that("PVE of irrelevant factors is close to zero under an overspecified Q", {
  skip_on_cran()
  skip_if_not_installed("Splinets")

  Q_true <- 2
  Q_fit <- 3
  N <- 60
  sim <- generate_bayesSYNC_data(N = N, n = rep(20, N), p = 10, Q = Q_true, L = 2,
                                 mu_func = function(t, j, p) rep(0, length(t)),
                                 vec_sd_eps = 0.2, a_om = 2, b_om = 2, seed = 2)
  fit <- suppressWarnings(bayesSYNC(sim$time_obs, sim$Y, L = 2, Q = Q_fit,
                                    n_g = 200, verbose = FALSE))

  expect_equal(get_factor_pve(fit), fit$factor_pve)

  res <- get_factor_pve(fit, factors = "all")
  total_sorted <- sort(res$total)
  expect_true(all(total_sorted[1:(Q_fit - Q_true)] < 1e-3))
  expect_true(all(total_sorted[(Q_fit - Q_true + 1):Q_fit] > 0.05))

  for (Phi in fit$list_list_Phi_hat) {
    expect_lt(gram_deviation(Phi, fit$time_g), 10 * max(diff(fit$time_g)))
  }
})
