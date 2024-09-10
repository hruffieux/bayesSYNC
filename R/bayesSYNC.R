#' BAYESian model for high-dimensional functional factor analysis of Shared latent dYNamiCs (bayesSYNC)
#'
#' This function performs functional analysis to uncover shared latent dynamics
#' using Bayesian methods.
#'
#' @param time_obs List of vectors or list of lists of vectors containing time
#'                 of observations for univariate or multivariate curves, respectively.
#' @param Y List of lists of vectors
#'          (for multivariate curves) with the functional data.
#' @param L Number of latent dimensions.
#' @param Q Number of latent factors.
#' @param K Number of spline functions (if applicable). Default is NULL.
#' @param list_hyper Hyperparameter settings constructed using the function
#'         \code{\link{set_hyper}}. If \code{NULL}, default hyperparameters will be used.
#' @param n_g Desired size for dense grid.
#' @param time_g Dense grid provided as a vector of size \code{n_g}. If provided,
#'               then \code{n_g} must be \code{NULL} as it will be taken to be
#'               \code{length(time_g)}.
#' @param tol_abs Tolerance on the absolute changes in the ELBO.
#' @param tol_rel Tolerance on the relative changes in the ELBO.
#' @param maxit Number of iterations for the mean-field variational Bayes algorithm. Default maxit = 500.
#' @param n_cpus User-sepcieid number of CPU to be used for parallel execution. Default n_cpus = 1 for serial execution.
#' @param verbose Boolean indicating whether messages should be printed during
#'                the run. Default is TRUE.
#' @param seed User-specified seed for reproducibility.
#'
#' @return An object containing the resulting estimates.
#'
#' @export
#'
# bayesSYNC <- function(time_obs, Y, L, K = NULL, list_hyper = NULL, maxit = 500,
#                   n_g = 1000, time_g = NULL, Phi_g = NULL, verbose = TRUE, seed = NULL) {
#   # Function implementation here
# }

bayesSYNC <- function(time_obs, Y, L, Q, K = NULL,
                      list_hyper = NULL,
                      n_g = 1000, time_g = NULL,
                      tol_abs = 1e-3,
                      tol_rel = 1e-5, maxit = 500, n_cpus = 1,
                      verbose = TRUE, seed = NULL) {

  check_structure(seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(seed)) {
    cat(paste0("== Seed set to ", seed, " ==\n\n"))
    set.seed(seed)
  }

  p <- length(Y[[1]])

  stopifnot(is.list(time_obs))
  stopifnot(is.list(Y))
  stopifnot(isTRUE(all.equal(length(time_obs), length(Y))))

  check_structure(L, "vector", "numeric", 1)
  check_natural(L)

  if (is.null(list_hyper)) {
    list_hyper <- set_hyper(d_0 = p)
  } else if (!inherits(list_hyper, "hyper")) {
    stop(paste0("The provided list_hyper must be an object of class ",
                "``hyper''. \n *** you must either use the ",
                "function set_hyper to set your own hyperparameters or ",
                "list_hyper to NULL for automatic choice. ***"))
  }

  check_structure(n_g, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(n_g)) check_natural(n_g)

  check_structure(time_g, "vector", "numeric", null_ok = TRUE)

  check_structure(maxit, "vector", "numeric", 1)
  check_natural(maxit)

  check_structure(verbose, "vector", "logical", 1)

  N <- length(time_obs)

  if (is.null(names(Y))) {
    stopifnot(is.null(names(time_obs)))
    subj_names <- paste0("subj_", 1:N)
    names(Y) <- names(time_obs) <- subj_names
  } else {
    subj_names <- names(Y)
    stopifnot(isTRUE(all.equal(names(time_obs), subj_names)) | is.null(names(time_obs)))
  }
  names(time_obs) <- subj_names

  p <-length(Y[[1]])

  if (is.null(K)) {  # Ruppert (2002) sets a simple default value for K as min(nobs/4,40), where nobs is the number of observations.
    # here since nobs differs for each i, we take the nobs / 4 = round(median(obs_i)/4), and do this for each variable j = 1, ..., p
    # and we enforce that K>=7
    K <- max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i))/4), 40)), 7)

    # if supplied K is such that length(K) = 1, then will be set to K <- rep(K, p)
  } else {
    check_structure(K, "vector", "numeric", 1)
  }
  check_natural(K)

  grid_obj <- get_grid_objects(time_obs, K, n_g = n_g, time_g = time_g,
                               format_univ = TRUE)

  C <- grid_obj$C
  n_g <- grid_obj$n_g
  time_g <- grid_obj$time_g
  C_g <- grid_obj$C_g


  tmp_names <- lapply(Y, function(Y_i) names(Y_i))
  tmp_names_time <-  lapply(time_obs, function(time_obs_i) names(time_obs_i))
  if (is.null(unlist(unique(tmp_names)))) {
    stopifnot(is.null(unlist(unique(tmp_names_time))))
    var_names <- paste0("variable_", 1:p)
    Y <- lapply(Y, function(Y_i) { names(Y_i) <- var_names; Y_i})
  } else if (!all_same(tmp_names)) {
    stop("Variable names in Y must be the same for all individuals.")
  } else {
    var_names <- names(Y[[1]])
  }

  if (!is.null(unlist(unique(tmp_names_time))) && (!all_same(tmp_names_time) | !isTRUE(all.equal(tmp_names_time[[1]], var_names)))) {
    stop("Variable names for each individual in time_obs must be the same as in Y.")
  }

  debug <- TRUE # whether to throw an error when the ELBO is not increasing monotonically

  bayesSYNC_core(N = N, p=p, L=L, Q=Q, K=K, C = C, Y = Y, list_hyper,
                 time_obs = time_obs, n_g =n_g, time_g = time_g, C_g = C_g,
                 tol_abs = tol_abs, tol_rel = tol_rel, maxit = maxit,
                 n_cpus = n_cpus, debug = debug, verbose = verbose)

}


bayesSYNC_core <- function(N, p, L,Q, K, C, Y, list_hyper, time_obs, n_g, time_g,
                           C_g, tol_abs, tol_rel, maxit, n_cpus, debug, verbose) {

  eps <- .Machine$double.eps^0.5
  eps_elbo <- 1e-11 # could probably be taken as eps, simply

  sigma_zeta <- list_hyper$sigma_zeta
  mu_beta <- list_hyper$mu_beta
  Sigma_beta <- list_hyper$Sigma_beta
  A <- list_hyper$A
  c_0 <- list_hyper$c_0
  d_0 <- list_hyper$d_0

  T_vec <- sapply(Y, length)
  sum_obs <- sum(sapply(time_obs, function(x) length(x)))

  inv_Sigma_zeta <- 1/sigma_zeta^2*diag(L)
  inv_Sigma_beta <- solve(Sigma_beta)

  mu_q_zeta <- lapply(1:Q, function(q) matrix(0.5, nrow = N, ncol = L))
  Sigma_q_zeta <- lapply(1:Q, function(q) lapply(1:N, function(i) diag(L)))

  mu_q_recip_sigsq_mu  <- mu_q_recip_sigsq_eps <- rep(1, p)
  mu_q_recip_sigsq_phi <- matrix(1, nrow = Q, ncol = L)

  kappa_q_sigsq_mu <- kappa_q_sigsq_phi <- K/2 + 0.5
  kappa_q_a_mu <- kappa_q_a_eps <- kappa_q_a_phi <- 1
  kappa_q_sigsq_eps <- sum_obs/2 + 0.5

  mu_q_normal_b <- matrix(rnorm(p*Q, 0, 1), nrow = p, ncol = Q)
  mu_q_b <- mu_q_normal_b
  Sigma_q_normal_b <- mu_q_gamma <- matrix(1, nrow= p, ncol= Q)

  mu_q_nu_phi <- lapply(1:Q, function(q) matrix(0.5, nrow = K+2, ncol = L))
  mu_q_recip_a_mu <- mu_q_recip_a_eps <- rep(1, p)
  mu_q_recip_a_phi <- matrix(1, nrow = Q, ncol = L)

  list_cp_C <- parallel::mclapply(1:N, function(i) crossprod(C[[i]]), mc.cores = n_cpus)
  sum_list_cp_C <- Reduce("+", list_cp_C)

  list_cp_C_Y <- parallel::mclapply(1:N, function(i) sapply(1:p, function(j) crossprod(C[[i]], Y[[i]][[j]])), mc.cores = n_cpus)
  list_cp_Y <- simplify2array(parallel::mclapply(1:p, function(j) sapply(1:N, function(i) crossprod(Y[[i]][[j]])), mc.cores = n_cpus))

  term_b <- (Sigma_q_normal_b + mu_q_normal_b^2)*mu_q_gamma

  inv_Sigma_q_nu_mu <- parallel::mclapply(1:p, function(j) blkdiag(inv_Sigma_beta,
                                                                   mu_q_recip_sigsq_mu[j]*diag(K)), mc.cores = n_cpus)

  inv_Sigma_q_nu_phi <- parallel::mclapply(1:Q, function(q) lapply(1:L, function(l) blkdiag(inv_Sigma_beta,
                                                                                              mu_q_recip_sigsq_phi[q,l]*diag(K))), mc.cores = n_cpus)

  Sigma_q_nu_phi <- lapply(1:Q, function(q) lapply(1:L, function(i) matrix(NA, nrow = K+2, ncol = K+2)))

  ELBO <- NULL
  for(i_iter in 1:maxit) {

    if (verbose) {
      cat(paste0("Iteration ", format(i_iter), "... \n"))
    }

    # Update of q(nu_mu):

    Sigma_q_nu_mu <- parallel::mclapply(1:p, function(j) solve(inv_Sigma_q_nu_mu[[j]] + mu_q_recip_sigsq_eps[j]* sum_list_cp_C), mc.cores = n_cpus)


    list_sweep <- parallel::mclapply(1:Q, function(q) lapply(1:N, function(i) sweep(C[[i]]%*%mu_q_nu_phi[[q]], 2, mu_q_zeta[[q]][i,], `*`)), mc.cores = n_cpus)

    mu_q_nu_mu <- parallel::mclapply(1:p, function(j) {

      sum_term_mu_j <- rowSums(sapply(1:N, function(i) {

        sum_val_i_j <- rowSums(sapply(1:Q, function(q) {
          rowSums(mu_q_b[j,q]*list_sweep[[q]][[i]])
        }))

        list_cp_C_Y[[i]][,j] - crossprod(C[[i]], sum_val_i_j)

      }))

      as.vector(Sigma_q_nu_mu[[j]] %*% (mu_q_recip_sigsq_eps[j]*sum_term_mu_j))

    }, mc.cores = n_cpus)


    # Update of q(nu_phi):

    sum_vec <- colSums(sweep(term_b, 1, mu_q_recip_sigsq_eps, `*`)) # vector of size Q

    E_q_zeta_square <- parallel::mclapply(1:N, function(i) lapply(1:Q, function(q) diag(Sigma_q_zeta[[q]][[i]]) + mu_q_zeta[[q]][i,]^2),  mc.cores = n_cpus)
    list_sum <- parallel::mclapply(1:Q, function(q) lapply(1:L, function(l) Reduce("+", lapply(1:N, function(i) E_q_zeta_square[[i]][[q]][l]*list_cp_C[[i]]))), mc.cores = n_cpus)

    mu_H_q_phi <- tr_term <- vector("list", length = Q)
    tr_qi <- matrix(NA, nrow = Q, ncol = N)

    list_cp_C_nu_mu <- parallel::mclapply(1:N, function(i) sapply(1:p, function(j) list_cp_C[[i]] %*% mu_q_nu_mu[[j]]), mc.cores = n_cpus)

    for (q in 1:Q) { # could be run in parallel

      term_q <- mu_q_recip_sigsq_eps * term_b[,q]
      sum_sigma <- sum(term_q)

      # do not move outside the above for loop as updated for q_tilde
      list_cp_C_nu_phi <- parallel::mclapply(1:N, function(i) lapply(1:Q, function(q) list_cp_C[[i]] %*% mu_q_nu_phi[[q]]), mc.cores = n_cpus)

      for (l in 1:L) {

        sum_term_mu_ql <- Reduce("+", parallel::mclapply(1:N, function(i) {

          tmp_term_qil <- tmp_term_qil_2 <- vector("list", length = Q)

          for (q_tilde in setdiff(1:Q, q)) {

            # avoid recomputing for each j
            tmp_term_qil[[q_tilde]] <- rowSums(sapply(1:L, function(l_tilde) { mu_q_zeta[[q_tilde]][i,l_tilde] * mu_q_zeta[[q]][i,l] *
                list_cp_C_nu_phi[[i]][[q_tilde]][,l_tilde, drop = F]
            }))
            tmp_term_qil[[q_tilde]] <- tmp_term_qil[[q_tilde]] * sum(mu_q_recip_sigsq_eps * mu_q_b[, q] * mu_q_b[, q_tilde])

          }

          # avoid recomputing for each j
          tmp_term_qil[[q]] <- rowSums(sapply(setdiff(1:L, l), function(l_tilde) { (mu_q_zeta[[q]][i,l_tilde]*mu_q_zeta[[q]][i,l] + Sigma_q_zeta[[q]][[i]][l,l_tilde])*
              list_cp_C[[i]] %*% mu_q_nu_phi[[q]][,l_tilde]
          }))

          # don't use the version below as leads to decreasing ELBO. The above version uses the updated mu_q_nu_phi[[q]] for l_tilde....
          # tmp_term_qil[[q]] <- rowSums(sapply(setdiff(1:L, l), function(l_tilde) { (mu_q_zeta[[q]][i,l_tilde]*mu_q_zeta[[q]][i,l] + Sigma_q_zeta[[q]][[i]][l,l_tilde])*
          #     list_cp_C_nu_phi[[i]][[q]][,l_tilde, drop = F]
          # }))

          tmp_term_qil[[q]] <- tmp_term_qil[[q]] * sum_sigma

          list_term_qil_2 <- Reduce("+", tmp_term_qil)

          list_term_qil_1 <- rowSums(sweep(mu_q_zeta[[q]][i,l]*(list_cp_C_Y[[i]] - list_cp_C_nu_mu[[i]]), 2, mu_q_b[,q]*mu_q_recip_sigsq_eps, "*"))

          as.vector(list_term_qil_1 -  list_term_qil_2)
        }, mc.cores = n_cpus))

      Sigma_q_nu_phi[[q]][[l]] <- solve(sum_vec[q] * list_sum[[q]][[l]] + inv_Sigma_q_nu_phi[[q]][[l]])
      mu_q_nu_phi[[q]][,l] <- Sigma_q_nu_phi[[q]][[l]]%*%sum_term_mu_ql
      }

      mu_H_q_phi[[q]] <- tr_term[[q]] <- vector("list", length = N)
      for (i in 1:N) {

        # # Update of q(zeta_iq):
        mod_h_i_nu <- rowSums(sapply(setdiff(1:Q, q), function(q_tilde) {
          tmp_iq_tilde <- mu_q_nu_phi[[q_tilde]] %*% mu_q_zeta[[q_tilde]][i,]
          rowSums(sapply(1:p, function(j) mu_q_b[j,q]*mu_q_recip_sigsq_eps[j]*mu_q_b[j,q_tilde]*tmp_iq_tilde))
        }))

        sum_mu <- rowSums(sapply(1:p, function(j) mu_q_b[j,q]*mu_q_recip_sigsq_eps[j]*(list_cp_C_Y[[i]][,j] - list_cp_C_nu_mu[[i]][,j]))) - list_cp_C[[i]] %*% mod_h_i_nu

        tr_term[[q]][[i]] <- diag(sapply(1:L, function(l) tr(list_cp_C[[i]] %*% Sigma_q_nu_phi[[q]][[l]])))
        mu_H_q_phi[[q]][[i]] <- crossprod(mu_q_nu_phi[[q]], list_cp_C[[i]]) %*% mu_q_nu_phi[[q]] + tr_term[[q]][[i]]
        Sigma_q_zeta[[q]][[i]]<- solve(sum_sigma*mu_H_q_phi[[q]][[i]] + inv_Sigma_zeta)
        mu_q_zeta[[q]][i,] <- as.vector(Sigma_q_zeta[[q]][[i]]%*%crossprod(mu_q_nu_phi[[q]], sum_mu))
        tr_qi[q, i] <- tr(mu_H_q_phi[[q]][[i]]%*%(Sigma_q_zeta[[q]][[i]]+tcrossprod(mu_q_zeta[[q]][i,])))

      }

    }

    list_tcp_nu_phi_zeta <- parallel::mclapply(1:Q, function(q) tcrossprod(mu_q_nu_phi[[q]], mu_q_zeta[[q]]), mc.cores = n_cpus)

    lambda_q_sigsq_eps <- mu_q_recip_a_eps + simplify2array(parallel::mclapply(1:p, function(j) {
      0.5 * sum(sapply(1:N, function(i) {
        list_cp_Y[i,j]-
          2*crossprod(mu_q_nu_mu[[j]] +
                        Reduce("+",
                               lapply(1:Q, function(q) mu_q_b[j,q]*list_tcp_nu_phi_zeta[[q]][,i])), list_cp_C_Y[[i]][,j]) +
          crossprod(mu_q_nu_mu[[j]], list_cp_C_nu_mu[[i]][,j]) + tr(crossprod(list_cp_C[[i]], Sigma_q_nu_mu[[j]])) +
          2*crossprod(Reduce("+",
                             lapply(1:Q, function(q) mu_q_b[j,q]*list_tcp_nu_phi_zeta[[q]][,i])), list_cp_C_nu_mu[[i]][,j]) +
          Reduce("+",
                 lapply(1:Q, function(q) {

                   sum_val_phi<- term_b[j, q]*tr_qi[q, i]

                   sum_val_phi_q_tilde <- Reduce("+", lapply(setdiff(1:Q, q), function(q_tilde) {
                     mu_q_b[j,q_tilde]*C[[i]]%*%list_tcp_nu_phi_zeta[[q_tilde]][,i]

                   }))

                   sum_val_phi + mu_q_b[j,q]*crossprod(C[[i]] %*% list_tcp_nu_phi_zeta[[q]][,i], sum_val_phi_q_tilde)

                 }))
      }))
    }, mc.cores = n_cpus))

    mu_q_recip_sigsq_eps <- kappa_q_sigsq_eps/lambda_q_sigsq_eps
    mu_q_log_sigsq_eps <- log(lambda_q_sigsq_eps)-digamma(kappa_q_sigsq_eps)

    lambda_q_a_eps <- mu_q_recip_sigsq_eps + 1/A^2
    mu_q_recip_a_eps <- kappa_q_a_eps/lambda_q_a_eps
    mu_q_log_a_eps <- log(lambda_q_a_eps)-digamma(kappa_q_a_eps)

    lambda_q_sigsq_mu <- unlist(parallel::mclapply(1:p, function(j) {
      mu_q_u_mu_j <- mu_q_nu_mu[[j]][-c(1:2)]
      Sigma_q_u_mu_j <- Sigma_q_nu_mu[[j]][-c(1:2), -c(1:2)]
      0.5*(crossprod(mu_q_u_mu_j)+tr(Sigma_q_u_mu_j))+mu_q_recip_a_mu[j]
    }, mc.cores = n_cpus))

    mu_q_recip_sigsq_mu <- kappa_q_sigsq_mu/lambda_q_sigsq_mu
    mu_q_log_sigsq_mu <- log(lambda_q_sigsq_mu) - digamma(kappa_q_sigsq_mu)

    lambda_q_a_mu <- mu_q_recip_sigsq_mu + 1/A^2
    mu_q_recip_a_mu <- kappa_q_a_mu/lambda_q_a_mu
    mu_q_log_a_mu <- log(lambda_q_a_mu)- digamma(kappa_q_a_mu)


    lambda_q_sigsq_phi <- mu_q_recip_a_phi + 0.5*sapply(1:L, function(l) sapply(1:Q, function(q) tr(Sigma_q_nu_phi[[q]][[l]][-c(1:2), -c(1:2)])+crossprod(mu_q_nu_phi[[q]][-c(1:2),l])))

    mu_q_recip_sigsq_phi <- kappa_q_sigsq_phi/ lambda_q_sigsq_phi
    mu_q_log_sigsq_phi <- log(lambda_q_sigsq_phi)-digamma(kappa_q_sigsq_phi)

    lambda_q_a_phi <- mu_q_recip_sigsq_phi + 1/A^2
    mu_q_recip_a_phi <- kappa_q_a_phi/lambda_q_a_phi
    mu_q_log_a_phi <- log(lambda_q_a_phi)- digamma(kappa_q_a_phi)

    cs_mu_q_gamma <- colSums(mu_q_gamma)

    c_1_omega <- c_0 + cs_mu_q_gamma
    d_1_omega <- d_0 + p - cs_mu_q_gamma

    dig <- digamma(c_0 + d_0 + p)
    mu_q_log_omega <- digamma(c_1_omega) - dig
    mu_q_log_1_omega <- digamma(d_1_omega) - dig

    # alpha was the precision of the slab on the B entries, now this precision is fixed to 1
    # so the code below is not used anymore # see commit d7d7b7cc09b33e581e7b9e9288da04e61d06d721
    # for version that includes the corresponding updates in all updates (commented)
    #
    # cs_term_b <- colSums(term_b)
    #
    # alpha_1_alpha <- alpha_0 + 0.5*cs_mu_q_gamma
    # beta_1_alpha <- beta_0 + 0.5*cs_term_b
    #
    # mu_q_alpha <- (alpha_0 + 0.5*cs_mu_q_gamma) / (beta_0 + 0.5*cs_term_b)
    # mu_q_log_alpha <- digamma(alpha_1_alpha)- log (beta_1_alpha)


    #Update of q(b_jq |gamma_jq)
    rs_tr_qi <- rowSums(tr_qi)
    Sigma_q_b <- matrix(NA, nrow=p, ncol=Q)
    # list_tcp_nu_phi_zeta <- parallel::mclapply(1:Q, function(q_tilde) tcrossprod(mu_q_nu_phi[[q_tilde]], mu_q_zeta[[q_tilde]]), mc.cores = n_cpus)

    for(q in 1:Q){

      sum_mu_q <- unlist(parallel::mclapply(1:p, function(j) sum(sapply(1:N, function(i) crossprod(list_tcp_nu_phi_zeta[[q]][,i],
                                                                                list_cp_C_Y[[i]][,j]-
                                                                                  list_cp_C_nu_mu[[i]][,j]-
                                                                                  list_cp_C[[i]]%*%rowSums(sapply(setdiff(1:Q, q), function(q_tilde) mu_q_gamma[j,q_tilde]*mu_q_normal_b[j,q_tilde]*list_tcp_nu_phi_zeta[[q_tilde]][,i]))))),
                                                    mc.cores = n_cpus))

      # don't move outside the for(q in 1:Q) loop as used in the update for sum_mu_q, leads to a less efficient scheme if outside
      Sigma_q_normal_b[,q] <- 1/(mu_q_recip_sigsq_eps * rs_tr_qi[q] + 1)
      mu_q_normal_b[,q] <- Sigma_q_normal_b[,q]*mu_q_recip_sigsq_eps*sum_mu_q
      mu_q_gamma[,q] <- 1 / (1 + sqrt(mu_q_recip_sigsq_eps*rs_tr_qi[q]+ 1) *
                           exp(mu_q_log_1_omega[q]-mu_q_log_omega[q] -
                                 0.5*(mu_q_normal_b[,q]^2)*(mu_q_recip_sigsq_eps*rs_tr_qi[q]+ 1)))

    }
    mu_q_b<- mu_q_gamma*mu_q_normal_b
    term_b <- (Sigma_q_normal_b + mu_q_normal_b^2)*mu_q_gamma
    Sigma_q_b <- term_b - mu_q_b^2

    inv_Sigma_q_nu_mu <- parallel::mclapply(1:p, function(j) blkdiag(inv_Sigma_beta,
                                                                           mu_q_recip_sigsq_mu[j]*diag(K)), mc.cores = n_cpus)

    inv_Sigma_q_nu_phi <- parallel::mclapply(1:Q, function(q) lapply(1:L, function(l) blkdiag(inv_Sigma_beta,
                                                                                              mu_q_recip_sigsq_phi[q,l]*diag(K))), mc.cores = n_cpus)

    # COMPUTE ELBO
    #
    elbo_y <- - sum(sum_obs/2*(log(2*pi) + mu_q_log_sigsq_eps) + mu_q_recip_sigsq_eps*(lambda_q_sigsq_eps - mu_q_recip_a_eps))

    vec_term_list_mu <- sapply(1:p, function(j) {
      log_det_j_obj <- determinant(Sigma_q_nu_mu[[j]], logarithm = TRUE)
      log_det_j_obj$modulus * log_det_j_obj$sign - crossprod(mu_q_nu_mu[[j]], inv_Sigma_q_nu_mu[[j]] %*% mu_q_nu_mu[[j]]) - tr(inv_Sigma_q_nu_mu[[j]] %*% Sigma_q_nu_mu[[j]])})

    # doesn't seem faster
    # vec_term_list_mu <- unlist(parallel::mclapply(1:p, function(j) {
    #   log_det_j_obj <- determinant(Sigma_q_nu_mu[[j]], logarithm = TRUE)
    #   log_det_j_obj$modulus * log_det_j_obj$sign - crossprod(mu_q_nu_mu[[j]], inv_Sigma_q_nu_mu[[j]] %*% mu_q_nu_mu[[j]]) - tr(inv_Sigma_q_nu_mu[[j]] %*% Sigma_q_nu_mu[[j]])}, mc.cores = n_cpus))

    elbo_mu <- -p*log(prod(diag(Sigma_beta)))/2 + sum(0.5*vec_term_list_mu - 0.5*K*mu_q_log_sigsq_mu + K/2 + 1 +
      K/2*mu_q_log_sigsq_mu - mu_q_recip_sigsq_mu*(mu_q_recip_a_mu-lambda_q_sigsq_mu)-
      0.5*mu_q_log_a_mu -kappa_q_sigsq_mu*log(lambda_q_sigsq_mu)-lgamma(0.5)+lgamma(kappa_q_sigsq_mu)+
      (3/2)*mu_q_log_a_mu - mu_q_recip_a_mu*(1/A^2 - lambda_q_a_mu) - log(lambda_q_a_mu) - lgamma(0.5) + lgamma(1) + 0.5*log(1/A^2)+
      (kappa_q_sigsq_eps -0.5)*mu_q_log_sigsq_eps - mu_q_recip_sigsq_eps*(mu_q_recip_a_eps - lambda_q_sigsq_eps)-
      0.5*mu_q_log_a_eps - kappa_q_sigsq_eps*log(lambda_q_sigsq_eps)-lgamma(0.5)+ lgamma(kappa_q_sigsq_eps)+
      (3/2)*mu_q_log_a_eps -mu_q_recip_a_eps*(1/A^2-lambda_q_a_eps)-
      log(lambda_q_a_eps)-lgamma(0.5)+lgamma(1)+ 0.5*log(1/A^2))

    elbo_b_g <- sum(0.5*mu_q_gamma*(log(Sigma_q_normal_b)+ 1) - 0.5*term_b +
      sweep(mu_q_gamma, 2, mu_q_log_omega, "*") +
      sweep(1 - mu_q_gamma, 2, mu_q_log_1_omega, "*") -
      mu_q_gamma*log(mu_q_gamma + eps_elbo) - (1-mu_q_gamma)*log(1- mu_q_gamma + eps_elbo))

    elbo_omega <- sum((c_0 - c_1_omega) * mu_q_log_omega + (d_0 - d_1_omega) * mu_q_log_1_omega +
      lbeta(c_1_omega, d_1_omega) - lbeta(c_0,d_0))

    elbo_sig_phi <- sum(K/2*mu_q_log_sigsq_phi - (mu_q_recip_a_phi - lambda_q_sigsq_phi)*mu_q_recip_sigsq_phi -
      0.5*mu_q_log_a_phi - kappa_q_sigsq_phi*log(lambda_q_sigsq_phi) - lgamma(0.5) + lgamma(kappa_q_sigsq_phi) +
      (3/2)*mu_q_log_a_phi - (1/A^2-lambda_q_a_phi)*mu_q_recip_a_phi +
      0.5*log(1/A^2) - log(lambda_q_a_phi) - lgamma(0.5) + lgamma(1))

    mat_term_list_phi <- sapply(1:L, function(l) sapply(1:Q, function(q) {
      log_det_q_l_obj <- determinant(Sigma_q_nu_phi[[q]][[l]], logarithm = TRUE)
      log_det_q_l_obj$modulus * log_det_q_l_obj$sign - crossprod(mu_q_nu_phi[[q]][,l], inv_Sigma_q_nu_phi[[q]][[l]] %*% mu_q_nu_phi[[q]][,l]) - tr(inv_Sigma_q_nu_phi[[q]][[l]] %*% Sigma_q_nu_phi[[q]][[l]])}))

    elbo_phi <- -Q*L*log(prod(diag(Sigma_beta)))/2 + sum(0.5*mat_term_list_phi - (K/2)*mu_q_log_sigsq_phi + K/2 + 1)

    elbo_zeta <- sum(sapply(1:Q, function(q)
      sum(sapply(1:N, function (i){
        log_det_i_q_obj <- determinant(Sigma_q_zeta[[q]][[i]], logarithm = TRUE)
        log_det_i_q <- log_det_i_q_obj$modulus * log_det_i_q_obj$sign

        0.5*log_det_i_q - L/2 * log(sigma_zeta) + L/2 - 0.5/sigma_zeta *(crossprod(mu_q_zeta[[q]][i,]) + tr(Sigma_q_zeta[[q]][[i]]))}))))

    ELBO_iter <- elbo_y + elbo_mu + elbo_b_g + elbo_omega + elbo_sig_phi + elbo_phi + elbo_zeta

    if (verbose) {
      cat(paste0("ELBO = ", format(ELBO_iter), "\n\n"))
    }

    ELBO <- c(ELBO, ELBO_iter)

    if (i_iter >1){

      ELBO_diff <- ELBO[i_iter] - ELBO[i_iter-1]

      if (ELBO_diff < -eps){

        if (debug) {
          stop(paste0("ELBO not increasing monotonically. Difference last consecutive ELBO values: ", ELBO_diff, " Exit."))
        } else {
          warning(paste0("ELBO not increasing monotonically. Difference last consecutive ELBO values: ", ELBO_diff))
        }
      }

      rel_converged <- (abs(ELBO[i_iter] / ELBO[i_iter-1] - 1) < tol_rel)
      abs_converged <- (abs(ELBO_diff) < tol_abs)

      if(rel_converged | abs_converged) {

        mess_criterion_met <- paste0(ifelse(rel_converged, "Relative ", "Absolute "), "convergence criterion met. \n")

        if (verbose) {
          cat(paste0("Convergence obtained after ", format(i_iter), " iterations. \n",
                     mess_criterion_met,
                     "Optimal marginal log-likelihood variational lower bound ",
                     "(ELBO) = ", format(ELBO[i_iter]), ". \n\n"))
        }

        break
      } else if (i_iter == maxit) {
        warning(paste0("Maximal number of iterations reached before convergence. Difference last consecutive ELBO values: ",
                       ELBO_diff, " Exit."))
      }

    }


  }

  sigsq_eps <- lambda_q_sigsq_eps/(kappa_q_sigsq_eps - 1)

  res_orth <- orthonormalise(N, p, Q, L, time_g, C_g, # see what she has used?
                             mu_q_nu_mu, # not used for the orthonormalisation
                             mu_q_b,
                             mu_q_zeta,
                             mu_q_nu_phi,
                             Sigma_q_zeta,
                             sigsq_eps)

  list_Y_hat <- res_orth$list_Y_hat
  list_Y_low <- res_orth$list_Y_low
  list_Y_upp <- res_orth$list_Y_upp
  list_gbl_hat <- res_orth$list_gbl_hat
  list_mu_hat <- res_orth$list_mu_hat
  list_list_Phi_hat <- res_orth$list_list_Phi_hat
  list_Zeta_hat <- res_orth$list_Zeta_hat
  list_Cov_zeta_hat <- res_orth$list_Cov_zeta_hat
  list_list_zeta_ellipse <- res_orth$list_list_zeta_ellipse

  if (!is.null(res_orth$mu_q_b_norm)) {
    mu_q_b  <- res_orth$mu_q_b_norm # QUESTION
  }

  B_hat <- mu_q_b
  ppi <- mu_q_gamma

  # probabilities of activity of each factor:
  factor_ppi <- 1 - colProds(1-ppi) # probability that each factor contains at least one contribution by a variable

  res <- create_named_list(K, list_Y_hat, list_Y_low, list_Y_upp,
                           list_mu_hat, list_list_Phi_hat,
                           list_Zeta_hat, list_Cov_zeta_hat, list_list_zeta_ellipse,
                           # Sigma_q_nu_mu,mu_q_nu_mu, Sigma_q_nu_phi,mu_q_nu_phi,mu_q_zeta, Sigma_q_zeta, sigsq_eps,
                           # lambda_q_sigsq_phi, lambda_q_a_phi, lambda_q_sigsq_mu, lambda_q_a_mu, lambda_q_sigsq_eps,
                           # lambda_q_a_eps, Sigma_q_normal_b,
                           B_hat,
                           ppi,
                           factor_ppi,
                           # time_g, C_g,
                           ELBO, i_iter, n_g)


}



orthonormalise <- function(N, p, Q, L, time_g, C_g, # see what she has used?
                           mu_q_nu_mu, # not used for the orthonormalisation
                           mu_q_b,
                           mu_q_zeta,
                           mu_q_nu_phi,
                           Sigma_q_zeta,
                           sigsq_eps) {
  # Orthogonalisation:

  list_mu_q_mu <- lapply(mu_q_nu_mu, function(mu_q_nu_mu_j) as.vector(C_g%*%mu_q_nu_mu_j))

  list_M_q_Phi <- lapply(mu_q_nu_phi, function(mu_q_nu_phi_q) {
    # M_q_V_phi <- Reduce(cbind, mu_q_nu_phi_q) # matrix K+2 x L
    M_q_Phi <- C_g %*% mu_q_nu_phi_q
  })

  # if(L > 1){
  #   list_M_q_Zeta <- lapply(1:Q, function(q) { t(sapply(mu_q_zeta, function(mu_q_zeta_i) mu_q_zeta_i[[q]]))})
  # } else {
  #   list_M_q_Zeta <- lapply(1:Q, function(q) { matrix(sapply(mu_q_zeta, function(mu_q_zeta_i) mu_q_zeta_i[[q]]), nrow = N)})
  # }
  if(L > 1) {
    list_M_q_Zeta <- mu_q_zeta
  } else {
    list_M_q_Zeta <- lapply(mu_q_zeta, function(ll) as.matrix(ll))
  }

  one_N <- rep(1, N)
  list_mu_mat <- lapply(list_mu_q_mu, function(mu_q_mu) tcrossprod(mu_q_mu, one_N))
  list_Y_mat <- lapply(1:p, function(j) { list_mu_mat[[j]] +
      Reduce('+', lapply(1:Q, function(q) mu_q_b[j, q] * tcrossprod(list_M_q_Phi[[q]], list_M_q_Zeta[[q]])))})

  bool_sim <- F # deal with identifiability. renormalise B and zeta. dirty - valid only for simulations, and B should be passed as an argument
  if (bool_sim) {
    # B is identifiable up to multiplicative sign on its columns
    norm_col_B <- sqrt(colSums(B^2)) # this is unknown in practice.
    norm_mu_q_b <- sqrt(colSums(mu_q_b^2))
    mu_q_b_norm <- sapply(1:Q, function(q) mu_q_b[,q] * norm_col_B[q] / norm_mu_q_b[q])
    list_M_q_Zeta <- lapply(1:Q, function(q) list_M_q_Zeta[[q]] * norm_mu_q_b[q] / norm_col_B[q])
  } else {
    mu_q_b_norm <- NULL
  }

  list_M_q_Phi_svd <- lapply(list_M_q_Phi, function(M_q_Phi) svd(M_q_Phi))
  list_U_orth <- lapply(list_M_q_Phi_svd, function(M_q_Phi_svd) M_q_Phi_svd$u)
  if (L > 1) {
    list_D_diag <- lapply(list_M_q_Phi_svd, function(M_q_Phi_svd) diag(M_q_Phi_svd$d))
  } else {
    list_D_diag <- lapply(list_M_q_Phi_svd, function(M_q_Phi_svd) matrix(M_q_Phi_svd$d))
  }
  list_V_orth <- lapply(list_M_q_Phi_svd, function(M_q_Phi_svd) M_q_Phi_svd$v)

  list_M_q_Zeta_rotn <- lapply(1:Q, function(q) list_M_q_Zeta[[q]] %*% list_V_orth[[q]] %*% list_D_diag[[q]])

  bool_center <- F
  if (bool_center) {
    list_m_zeta <- lapply(list_M_q_Zeta_rotn, function(M_q_Zeta_rotn) apply(M_q_Zeta_rotn, 2, mean))
    list_eigen_M_q_Zeta_shift <- lapply(1:Q, function(q) eigen(cov(list_M_q_Zeta_rotn[[q]] - tcrossprod(one_N, list_m_zeta[[q]]))))

    add_to_mean <- lapply(1:Q, function(q) as.vector(list_U_orth[[q]] %*% list_m_zeta[[q]]))
    for (j in 1:p) {
      list_mu_q_mu[[j]] <- list_mu_q_mu[[j]] + Reduce("+", lapply(1:Q, function(q) mu_q_b[j, q]* add_to_mean[[q]]))
    }

  } else {
    list_eigen_M_q_Zeta_shift <- lapply(list_M_q_Zeta_rotn, function(M_q_Zeta_rotn) eigen(cov(M_q_Zeta_rotn)))
  }

  list_Q <- lapply(list_eigen_M_q_Zeta_shift, function(eigen_M_q_Zeta_shift) eigen_M_q_Zeta_shift$vectors)
  if (L > 1) {
    list_Lambda <- lapply(list_eigen_M_q_Zeta_shift, function(eigen_M_q_Zeta_shift) diag(eigen_M_q_Zeta_shift$values + 1e-10))
    list_Lambda_inv <- lapply(list_eigen_M_q_Zeta_shift, function(eigen_M_q_Zeta_shift) diag(1/(eigen_M_q_Zeta_shift$values + 1e-10)))
  } else {
    list_Lambda <- lapply(list_eigen_M_q_Zeta_shift, function(eigen_M_q_Zeta_shift) matrix(eigen_M_q_Zeta_shift$values + 1e-10))
    list_Lambda_inv <- lapply(list_eigen_M_q_Zeta_shift, function(eigen_M_q_Zeta_shift) matrix(1/(eigen_M_q_Zeta_shift$values + 1e-10)))
  }

  list_S <- lapply(1:Q, function(q) list_Q[[q]]%*%sqrt(list_Lambda[[q]]))
  list_S_inv <- lapply(1:Q, function(q) tcrossprod(sqrt(list_Lambda_inv[[q]]), list_Q[[q]]))

  list_Phi_hat <- lapply(1:Q, function(q) list_U_orth[[q]]%*%list_S[[q]])

  if (bool_center) {
    list_Zeta_hat <- lapply(1:Q, function(q) tcrossprod(list_M_q_Zeta_rotn[[q]] - tcrossprod(one_N, list_m_zeta[[q]]), list_S_inv[[q]]))
  } else {
    list_Zeta_hat <- lapply(1:Q, function(q) tcrossprod(list_M_q_Zeta_rotn[[q]], list_S_inv[[q]]))

    # q <- 1
    # print(all.equal(tcrossprod(list_Phi_hat[[q]], list_Zeta_hat[[q]]),
    #                 tcrossprod(list_M_q_Phi[[q]], list_M_q_Zeta[[q]])))
    # q <- 2
    # print(all.equal(tcrossprod(list_Phi_hat[[q]], list_Zeta_hat[[q]]),
    #                 tcrossprod(list_M_q_Phi[[q]], list_M_q_Zeta[[q]])))
  }

  norm_const <- matrix(NA, nrow = Q, ncol = L)
  for (q in 1:Q) {
    for(l in 1:L) {

      norm_const[q, l] <- sqrt(trapint(time_g, (list_Phi_hat[[q]][,l])^2))
      if(norm_const[q, l]!=0) {

        list_Phi_hat[[q]][,l] <- list_Phi_hat[[q]][,l]/norm_const[q, l]
        list_Zeta_hat[[q]][,l] <- norm_const[q, l]*list_Zeta_hat[[q]][,l]

        # if(!is.null(Phi_g)) { # TODO: we also need to know which simulated pathway corresponds to which inferred pathway...
        #
        #   cprod_sign <- sign(cprod(Phi_hat[,l], Phi_g[,l]))
        #   Phi_hat[,l] <- cprod_sign*Phi_hat[,l]
        #   Zeta_hat[,l] <- cprod_sign*Zeta_hat[,l]
        #
        # }
      }
    }
  }

  # normalise columns of B
  # B_hat <- apply(mu_q_b, 2, function(mu_q_b_q) mu_q_b_q / sum(abs(mu_q_b_q)))


  list_mu_hat <- vector("list", length = p)
  for (j in 1:p) {
    list_mu_hat[[j]] <- list_mu_q_mu[[j]]
  }

  list_Cov_zeta_hat <- vector("list", length = Q)
  list_list_zeta_ellipse <-  vector("list", length=Q)
  list_list_Phi_hat <- vector("list", length = Q)
  for (q in 1:Q) {

    if (L > 1) {
      scale_mat <- diag(norm_const[q,])
    } else {
      scale_mat <- matrix(norm_const[q,])
    }
    # list_Cov_zeta_hat[[q]] <- vector("list", length = N)
    # for(i in 1:N) {
    #   mat_transform <- list_S_inv[[q]]%*%tcrossprod(list_D_diag[[q]], list_V_orth[[q]])
    #   Cov_zeta_dot_hat <- tcrossprod(mat_transform%*%Sigma_q_zeta[[q]][[i]], mat_transform)
    #   list_Cov_zeta_hat[[q]][[i]] <- tcrossprod(scale_mat %*% Cov_zeta_dot_hat, scale_mat)
    # }

    list_Cov_zeta_hat[[q]] <- lapply(1:N, function(i) {
      mat_transform <- list_S_inv[[q]] %*% tcrossprod(list_D_diag[[q]], list_V_orth[[q]])
      Cov_zeta_dot_hat <- tcrossprod(mat_transform%*%Sigma_q_zeta[[q]][[i]], mat_transform)
      tcrossprod(scale_mat %*% Cov_zeta_dot_hat, scale_mat)
    })

    if (L > 1) {
      list_list_zeta_ellipse[[q]] <-  vector("list", length=N)
      for(i in 1:N) {

        zeta_mean <- list_Zeta_hat[[q]][i,][1:2]

        zeta_ellipse <- ellipse::ellipse(
          list_Cov_zeta_hat[[q]][[i]][1:2, 1:2],
          centre = zeta_mean,
          level = 0.95
        )

        list_list_zeta_ellipse[[q]][[i]] <- zeta_ellipse

      }
    } else {
      list_list_zeta_ellipse[[q]] <- NULL
    }


    list_list_Phi_hat[[q]] <- list_Phi_hat[[q]]

  }

  list_Y_hat <- list_Y_low <- list_Y_upp <- vector("list", length = N)

  # list_var_vec <- lapply(1:Q, function(q) lapply(1:N, function(i)
  #   diag(tcrossprod(list_M_q_Phi[[q]]%*%Sigma_q_zeta[[q]][[i]], list_M_q_Phi[[q]])))) # functions assumed to be known exactly (we can use all objects pre-orthogonalisation to construct y)

  # more efficient
  list_var_vec <- lapply(1:Q, function(q) lapply(1:N, function(i)
    rowSums((list_M_q_Phi[[q]]%*%Sigma_q_zeta[[q]][[i]]) * list_M_q_Phi[[q]])))

  # print(isTRUE(all.equal(list_var_vec, list_var_vec_2)))
  # list_sd_vec <- lapply(1:Q, function(q) lapply(1:N, function(i)
  #   sqrt(diag(tcrossprod(list_list_Phi_hat[[q]]%*%list_Cov_zeta_hat[[q]][[i]], list_list_Phi_hat[[q]]))))) # functions assumed to be known exactly

  for(i in 1:N) {

    list_Y_hat[[i]] <- list_Y_low[[i]] <- list_Y_upp[[i]] <- vector("list", length = p)

    for (j in 1:p) {
      sd_i_j <- sqrt(Reduce("+", lapply(1:Q, function(q) mu_q_b[j, q]^2 * list_var_vec[[q]][[i]])) + sigsq_eps[j])# B assumed to be known exactly (variance = b^T Phi^T Cov_zeta Phi b) so sqrt(b^2) for sd


      list_Y_hat[[i]][[j]] <- list_Y_mat[[j]][,i]
      list_Y_low[[i]][[j]] <- list_Y_mat[[j]][,i] + qnorm(0.025) * sd_i_j
      list_Y_upp[[i]][[j]] <- list_Y_mat[[j]][,i] + qnorm(0.975) * sd_i_j
    }

  }

  create_named_list(list_Y_hat, list_Y_low, list_Y_upp,
                    list_mu_hat, list_list_Phi_hat,
                    list_Zeta_hat, list_Cov_zeta_hat, list_list_zeta_ellipse,
                    mu_q_b_norm)
}
