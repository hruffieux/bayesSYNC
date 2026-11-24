# Regression tests for a single FPCA component per factor (L = 1), which
# previously failed in bayesSYNC_core() (empty within-factor sum and diag() of a
# scalar), in the post-processing (NULL ellipse elements deleted from the list)
# and in generate_bayesSYNC_data() (too few knots for the cubic spline systems).

test_that("generate_bayesSYNC_data() works with L = 1", {
  skip_if_not_installed("Splinets")

  for (Q in 1:4) {
    dat <- generate_bayesSYNC_data(N = 10, n = rep(8, 10), p = 4, Q = Q, L = 1,
                                   mu_func = function(t, j, p) 0, n_g = 50,
                                   seed = Q)
    expect_length(dat$Zeta, Q)
    expect_equal(dim(dat$Zeta[[1]]), c(10, 1))
    expect_equal(dim(dat$Phi_g[[1]]), c(50, 1))
  }
})

test_that("bayesSYNC() runs with L = 1", {
  set.seed(1)
  N <- 15
  p <- 4
  Q <- 2
  time_obs <- lapply(1:N, function(i) sort(runif(10)))
  Y <- lapply(1:N, function(i) {
    zeta_i <- rnorm(Q)
    lapply(1:p, function(j) {
      j * zeta_i[1] * sin(2 * pi * time_obs[[i]]) +
        zeta_i[2] * cos(2 * pi * time_obs[[i]]) * (j <= 2) +
        rnorm(10, sd = 0.3)
    })
  })

  fit <- suppressWarnings(
    bayesSYNC(time_obs, Y, L = 1, Q = Q, K = 8, anneal = c(1, 1.5, 6),
              n_g = 50, maxit = 20, verbose = FALSE, seed = 1)
  )

  expect_equal(fit$L, 1)
  expect_length(fit$list_Zeta_hat, Q)
  expect_equal(dim(fit$list_Zeta_hat[[1]]), c(N, 1))
  expect_equal(dim(fit$list_list_Phi_hat[[1]]), c(50, 1))
  expect_equal(dim(fit$list_Cov_zeta_hat[[1]][[1]]), c(1, 1))
  expect_length(fit$list_list_zeta_ellipse, Q)
  expect_null(fit$list_list_zeta_ellipse[[1]])
  expect_true(is.finite(fit$ELBO_iter))
})
