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
#' @param Psi_g Reference eigenfunctions (if available, e.g., in simulations)
#'              used to flip the sign of the resulting scores and eigenfunctions.
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
#                   n_g = 1000, time_g = NULL, Psi_g = NULL, verbose = TRUE, seed = NULL) {
#   # Function implementation here
# }

bayesSYNC <- function(time_obs, Y, L, Q, K = NULL,
                      list_hyper = NULL,
                      n_g = 1000, time_g = NULL,
                      Psi_g = NULL, tol_abs = 1e-3,
                      tol_rel = 1e-5, maxit = 500, n_cpus = 1, verbose = TRUE, seed = NULL) {

  check_structure(seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(seed)) {
    cat(paste0("== Seed set to ", seed, " ==\n\n"))
    set.seed(seed)
  }

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


  bayesSYNC_core(N= N, p=p, L=L, Q=Q, K=K, C = C, Y = Y, list_hyper,
                      time_obs = time_obs,n_g =n_g,
                      time_g =time_g,C_g =C_g, tol_abs = tol_abs,
                      tol_rel = tol_rel, maxit = maxit, n_cpus = n_cpus, verbose = verbose)

}


bayesSYNC_core <- function(N, p, L,Q, K, C, Y, list_hyper, time_obs, n_g,
                           time_g, C_g, tol_abs, tol_rel, maxit, n_cpus, verbose) {

  eps <- .Machine$double.eps^0.5

  sigma_zeta <- list_hyper$sigma_zeta
  mu_beta <- list_hyper$mu_beta
  Sigma_beta <- list_hyper$Sigma_beta
  A <- list_hyper$A
  c_0 <- list_hyper$c_0
  d_0 <- list_hyper$d_0
  alpha_0 <- list_hyper$alpha_0
  beta_0 <- list_hyper$beta_0

  T_vec <- sapply(Y, length)

  inv_Sigma_zeta <- 1/sigma_zeta^2*diag(L)
  inv_Sigma_beta <- solve(Sigma_beta)

  # mu_q_zeta <- vector("list", length = N)
  # Sigma_q_zeta <- vector("list", length = N)

  # mu_q_zeta <- lapply(1:N, function(i) lapply(1:Q, function(q) return(rep(0.5, L))))
  mu_q_zeta <- lapply(1:Q, function(q) matrix(0.5, nrow = N, ncol = L))
  Sigma_q_zeta <- lapply(1:N, function(i) lapply(1:Q, function(q) diag(L)))

  mu_q_recip_sigsq_mu  <- mu_q_recip_sigsq_eps <- rep(1, p)
  mu_q_recip_sigsq_phi <- matrix(1, nrow = Q, ncol = L)

  kappa_q_sigsq_mu <- kappa_q_sigsq_phi <- K/2 + 0.5
  kappa_q_a_mu <- kappa_q_a_eps <- kappa_q_a_phi <- 1
  kappa_q_sigsq_eps <- sum(sapply(time_obs, function(x) length (x)))/2 + 0.5

  mu_q_normal_b <- matrix(rnorm(p*Q, 0, 1), nrow = p, ncol = Q)
  mu_q_b <- mu_q_normal_b
  Sigma_q_normal_b <- mu_q_gamma <- matrix(1, nrow= p, ncol= Q)

  mu_q_alpha <- rep(1, Q)

  # mu_q_nu_phi <- lapply(1:Q, function(q) lapply(1:L, function(l) rep(0.5, K+2)))
  mu_q_nu_phi <- lapply(1:Q, function(q) matrix(0.5, nrow = K+2, ncol = L))
  mu_q_recip_a_mu <- mu_q_recip_a_eps <- rep(1, p)
  mu_q_recip_a_phi <- matrix(1, nrow = Q, ncol = L)

  list_cp_C <- lapply(1:N, function(i) crossprod(C[[i]]))
  list_cp_C_Y <- lapply(1:N, function(i) sapply(1:p, function(j) crossprod(C[[i]], Y[[i]][[j]])))
  list_cp_Y <- sapply(1:p, function(j) sapply(1:N, function(i) crossprod(Y[[i]][[j]])))


  ELBO <- NULL
  for(i_iter in 1:maxit) {

    cat("Iteration", i_iter, "\n")

    # Update of q(nu_mu):
    Sigma_q_nu_mu <- vector("list", length = p)
    mu_q_nu_mu <- vector("list", length = p)

    sum_list_cp_C <- Reduce("+", list_cp_C)

    Sigma_q_nu_mu <- lapply(1:p, function(j) solve(blkdiag(inv_Sigma_beta,
                                                           mu_q_recip_sigsq_mu[j]*diag(K)) + mu_q_recip_sigsq_eps[j]* sum_list_cp_C))

    mu_q_nu_mu <- parallel::mclapply(1:p, function(j) {

      sum_term_mu_j <- rowSums(sapply(1:N, function(i) {

        sum_val_i_j <- rowSums(sapply(1:Q, function(q) {
          rowSums(mu_q_b[j,q]*sweep(C[[i]]%*%mu_q_nu_phi[[q]], 2, mu_q_zeta[[q]][i,], `*`))
        }))

        list_cp_C_Y[[i]][,j] - crossprod(C[[i]], sum_val_i_j)

      }))

      as.vector(Sigma_q_nu_mu[[j]] %*% (mu_q_recip_sigsq_eps[j]*sum_term_mu_j))

    }, mc.cores = n_cpus)



    # Update of q(nu_phi):

    sum_vec <- colSums(sweep((Sigma_q_normal_b + mu_q_normal_b^2)*mu_q_gamma, 1, mu_q_recip_sigsq_eps, `*`)) # vector of size Q

    E_q_inv_Sigma_nu_phi <- lapply(1:Q, function(q) lapply(1:L, function(l) blkdiag(inv_Sigma_beta,
                                                                                    mu_q_recip_sigsq_phi[q,l]*diag(K))))

    E_q_zeta_square <- lapply(1:N, function(i) lapply(1:Q, function(q) diag(Sigma_q_zeta[[i]][[q]]) + mu_q_zeta[[q]][i,]^2))
    list_sum <- lapply(1:Q, function(q) lapply(1:L, function(l) Reduce("+", lapply(1:N, function(i) E_q_zeta_square[[i]][[q]][l]*list_cp_C[[i]]))))


    for (q in 1:Q) {

      term_q <- mu_q_recip_sigsq_eps * (mu_q_normal_b[, q]^2 + Sigma_q_normal_b[, q]) * mu_q_gamma[,q]
      sum_sigma <- sum(term_q)

      for (l in 1:L) {

        sum_term_mu_ql <- Reduce("+", parallel::mclapply(1:N, function(i) {

          tmp_term_qil <- vector("list", length = Q)

          for (q_tilde in setdiff(1:Q, q)) {

            # avoid recomputing for each j
            tmp_term_qil[[q_tilde]] <- rowSums(sapply(1:L, function(l_tilde) { mu_q_zeta[[q_tilde]][i,l_tilde] * mu_q_zeta[[q]][i,l] *
                list_cp_C[[i]] %*% mu_q_nu_phi[[q_tilde]][,l_tilde]
            }))

            tmp_term_qil[[q_tilde]] <- tmp_term_qil[[q_tilde]] * sum(mu_q_recip_sigsq_eps * mu_q_b[, q] * mu_q_b[, q_tilde])

          }

          # avoid recomputing for each j
          tmp_term_qil[[q]] <- rowSums(sapply(setdiff(1:L, l), function(l_tilde) { (mu_q_zeta[[q]][i,l_tilde]*mu_q_zeta[[q]][i,l] + Sigma_q_zeta[[i]][[q]][l,l_tilde])*
              list_cp_C[[i]] %*% mu_q_nu_phi[[q]][,l_tilde]
          }))

          # tmp_term_qil[[q]] <- tmp_term_qil[[q]] * colSums(crossprod(mu_q_recip_sigsq_eps * (mu_q_normal_b[, q]^2 + Sigma_q_normal_b[, q]), mu_q_gamma[, q]))
          tmp_term_qil[[q]] <- tmp_term_qil[[q]] * sum_sigma

          list_term_qil_2 <- Reduce("+", tmp_term_qil)

          list_term_qil_1 <- rowSums(sapply(1:p, function(j) mu_q_recip_sigsq_eps[j]*(mu_q_b[j,q]*mu_q_zeta[[q]][i,l]*(list_cp_C_Y[[i]][,j] - list_cp_C[[i]]%*% mu_q_nu_mu[[j]]))))

          as.vector(list_term_qil_1 -  list_term_qil_2)
        }, mc.cores = n_cpus))

      Sigma_q_nu_phi[[q]][[l]] <- solve(sum_vec[q] * list_sum[[q]][[l]] + E_q_inv_Sigma_nu_phi[[q]][[l]])
      mu_q_nu_phi[[q]][,l] <- Sigma_q_nu_phi[[q]][[l]]%*%sum_term_mu_ql
      }

      new_version <- T
      # if (new_version) {
      #
      #   for (i in 1:N) {
      #
      #     mod_h_i_nu <- rowSums(sapply(setdiff(1:Q, q), function(q_tilde) {
      #       tmp_iq_tilde <- mu_q_nu_phi[[q_tilde]] %*% mu_q_zeta[[q_tilde]][i,]
      #       rowSums(sapply(1:p, function(j) mu_q_b[j,q]*mu_q_recip_sigsq_eps[j]*mu_q_b[j,q_tilde]*tmp_iq_tilde))
      #     }))
      #
      #     sum_mu <- rowSums(sapply(1:p, function(j) mu_q_b[j,q]*mu_q_recip_sigsq_eps[j]*(Y[[i]][[j]] - C[[i]] %*% mu_q_nu_mu[[j]]))) - C[[i]]%*%mod_h_i_nu
      #
      #     tr_term <- diag(sapply(1:L, function(l) tr(list_cp_C[[i]] %*% Sigma_q_nu_phi[[q]][[l]])))
      #     mu_H_q_phi <- crossprod(mu_q_nu_phi[[q]], list_cp_C[[i]]) %*% mu_q_nu_phi[[q]] + tr_term
      #     Sigma_q_zeta[[i]][[q]]<- solve(sum_sigma*mu_H_q_phi + inv_Sigma_zeta)
      #     mu_q_zeta[[q]][i,] <- as.vector(Sigma_q_zeta[[i]][[q]]%*%crossprod(C[[i]] %*% mu_q_nu_phi[[q]], sum_mu))
      #   }
      # }


      if (new_version) {

        for (i in 1:N) {

          mod_h_i_nu <- rowSums(sapply(setdiff(1:Q, q), function(q_tilde) {
            tmp_iq_tilde <- mu_q_nu_phi[[q_tilde]] %*% mu_q_zeta[[q_tilde]][i,]
            rowSums(sapply(1:p, function(j) mu_q_b[j,q]*mu_q_recip_sigsq_eps[j]*mu_q_b[j,q_tilde]*tmp_iq_tilde))
          }))

          sum_mu <- rowSums(sapply(1:p, function(j) mu_q_b[j,q]*mu_q_recip_sigsq_eps[j]*(list_cp_C_Y[[i]][,j] - list_cp_C[[i]] %*% mu_q_nu_mu[[j]]))) - list_cp_C[[i]]%*%mod_h_i_nu

          tr_term <- diag(sapply(1:L, function(l) tr(list_cp_C[[i]] %*% Sigma_q_nu_phi[[q]][[l]])))
          mu_H_q_phi <- crossprod(mu_q_nu_phi[[q]], list_cp_C[[i]]) %*% mu_q_nu_phi[[q]] + tr_term
          Sigma_q_zeta[[i]][[q]]<- solve(sum_sigma*mu_H_q_phi + inv_Sigma_zeta)
          mu_q_zeta[[q]][i,] <- as.vector(Sigma_q_zeta[[i]][[q]]%*%crossprod(mu_q_nu_phi[[q]], sum_mu))
        }
      }


   }




   # # Update of q(zeta_iq):

    if (!new_version) {
      for (i in 1:N){
        for (q in 1:Q){
          sum_sigma <-0
          sum_mu <- rep(0, length(time_obs[[i]]))
          for (j in 1:p){
            sum_sigma<- sum_sigma+(mu_q_normal_b[j,q]^2+Sigma_q_normal_b[j,q])*mu_q_gamma[j,q]*mu_q_recip_sigsq_eps[j]
            h_i_nu <- rep(0, K+2)
            for ( q_tilde in 1:Q){
              if (q_tilde != q){
                mu_V_q_tilde_phi <- matrix(NA, nrow = K+2, ncol = L)
                for (l in 1:L){
                  mu_V_q_tilde_phi[,l]<- mu_q_nu_phi[[q_tilde]][,l]
                }
                # term_i_nu <- mu_q_b[j, q_tilde] * mu_V_q_tilde_phi%*%mu_q_zeta[[i]][[q_tilde]]
                term_i_nu <- mu_q_b[j, q_tilde] * mu_V_q_tilde_phi%*%mu_q_zeta[[q_tilde]][i,]
                h_i_nu <- h_i_nu + term_i_nu
              }
            }
            sum_mu <- sum_mu + mu_q_b[j,q]*mu_q_recip_sigsq_eps[j]*(Y[[i]][[j]]-C[[i]]%*% mu_q_nu_mu[[j]]- C[[i]]%*%h_i_nu)
          }

          tr_term <- diag(L)
          mu_V_q_phi <- matrix(NA,nrow = K+2, ncol = L)
          for (l in 1:L){
            mu_V_q_phi[,l]<- mu_q_nu_phi[[q]][,l]
            tr_term[l,l] <- tr(list_cp_C[[i]]%*% Sigma_q_nu_phi[[q]][[l]])
          }
          mu_H_q_phi <- t(mu_V_q_phi)%*%list_cp_C[[i]]%*%mu_V_q_phi + tr_term
          Sigma_q_zeta[[i]][[q]]<- solve(sum_sigma*mu_H_q_phi + inv_Sigma_zeta)
          # mu_q_zeta[[i]][[q]] <- as.vector(Sigma_q_zeta[[i]][[q]]%*%(t(mu_V_q_phi)%*%t(C[[i]])%*%sum_mu))
          mu_q_zeta[[q]][i,] <- as.vector(Sigma_q_zeta[[i]][[q]]%*%(t(mu_V_q_phi)%*%t(C[[i]])%*%sum_mu))
        }
      }
    }

    # Update q(sigsa_eps_j) & Update of q(a_eps)
    #Update of q(sigma_mu) & Update of q(a_mu)

    mu_q_recip_a_eps <- rep(1, p)
    lambda_q_sigsq_eps <- rep(NA,p)
    lambda_q_a_eps <- rep(NA,p)
    mu_q_log_sigsq_eps <- rep(NA,p)
    mu_q_log_a_eps <- rep(NA,p)

    mu_q_recip_a_mu <- rep(1, p)
    lambda_q_sigsq_mu <- rep(NA,p)
    lambda_q_a_mu <- rep(NA,p)
    mu_q_log_sigsq_mu <- rep(NA,p)
    mu_q_log_a_mu <- rep(NA,p)

    for (j in 1:p){
      sum_val_y <-0
      sum_val_mu <- 0
      sum_val_phi <-0
      for (i in 1:N){
        sum_val_y <- sum_val_y+ list_cp_Y[i,j]-2*t(mu_q_nu_mu[[j]])%*%list_cp_C_Y[[i]][,j]
        sum_val_mu <- sum_val_mu+ t(mu_q_nu_mu[[j]])%*%list_cp_C[[i]]%*%mu_q_nu_mu[[j]]+ tr(list_cp_C[[i]]%*% Sigma_q_nu_mu[[j]])
        for (q in 1:Q){
          mu_V_q_phi <- matrix(NA,nrow = K+2, ncol = L)
          tr_term <- diag(L)
          for (l in 1:L){
            mu_V_q_phi[,l]<- mu_q_nu_phi[[q]][,l]
            tr_term[l,l] <- tr(list_cp_C[[i]]%*% Sigma_q_nu_phi[[q]][[l]])
          }
          mu_H_q_phi <- t(mu_V_q_phi)%*%list_cp_C[[i]]%*%mu_V_q_phi + tr_term
          # term_val_y <- 2*mu_q_b[j,q]*t(mu_q_zeta[[i]][[q]])%*%t(mu_V_q_phi)%*%list_cp_C_Y[[i]][,j]
          term_val_y <- 2*mu_q_b[j,q]*t(mu_q_zeta[[q]][i,])%*%t(mu_V_q_phi)%*%list_cp_C_Y[[i]][,j]
          sum_val_y <- sum_val_y - term_val_y

          sum_val_mu <- sum_val_mu+2*mu_q_b[j,q]*t(mu_q_zeta[[q]][i,])%*%t(mu_V_q_phi)%*%list_cp_C[[i]]%*%mu_q_nu_mu[[j]]
          sum_val_phi<- sum_val_phi + (mu_q_normal_b[j,q]^2+Sigma_q_normal_b[j,q])*mu_q_gamma[j,q]*
            tr(mu_H_q_phi%*%(Sigma_q_zeta[[i]][[q]]+tcrossprod(mu_q_zeta[[q]][i,])))
          sum_val_phi_q_tilde <- rep(0, length(time_obs[[i]]))
          for (q_tilde in 1:Q){
            if (q_tilde !=q){
              mu_V_q_tilde_phi<- matrix(NA,nrow = K+2, ncol = L)
              for (l in 1:L){
                mu_V_q_tilde_phi[,l]<- mu_q_nu_phi[[q_tilde]][,l]
              }
              sum_val_phi_q_tilde<- sum_val_phi_q_tilde+mu_q_b[j,q_tilde]*C[[i]]%*%mu_V_q_tilde_phi%*%mu_q_zeta[[q_tilde]][i,]
            }
          }
          sum_val_phi <- sum_val_phi +mu_q_b[j,q]*t(mu_q_zeta[[q]][i,])%*%t(mu_V_q_phi)%*%t(C[[i]])%*%sum_val_phi_q_tilde
        }
      }
      lambda_q_sigsq_eps[j]<- mu_q_recip_a_eps[j]+ 0.5*(unlist(sum_val_y+sum_val_mu+sum_val_phi))
      mu_q_recip_sigsq_eps[j]<- kappa_q_sigsq_eps/lambda_q_sigsq_eps[j]
      #mu_q_recip_sigsq_eps[j]<- 1
      mu_q_log_sigsq_eps[j]<- log(lambda_q_sigsq_eps[j])-digamma(kappa_q_sigsq_eps)

      lambda_q_a_eps[j]<- mu_q_recip_sigsq_eps[j] + 1/A^2
      mu_q_recip_a_eps[j] <- kappa_q_a_eps/lambda_q_a_eps[j]
      mu_q_log_a_eps[j]<- log(lambda_q_a_eps[j])-digamma(kappa_q_a_eps)


      mu_q_u_mu <- mu_q_nu_mu[[j]][-c(1:2)]
      Sigma_q_u_mu <- Sigma_q_nu_mu[[j]][-c(1:2), -c(1:2)]

      lambda_q_sigsq_mu[j] <-  0.5*(t(mu_q_u_mu)%*%mu_q_u_mu+tr(Sigma_q_u_mu))+mu_q_recip_a_mu[j]
      mu_q_recip_sigsq_mu[j]<- kappa_q_sigsq_mu/lambda_q_sigsq_mu[j]
      mu_q_log_sigsq_mu[j]<- log(lambda_q_sigsq_mu[j])-digamma(kappa_q_sigsq_mu)

      lambda_q_a_mu[j]<- mu_q_recip_sigsq_mu[j]+ 1/A^2
      mu_q_recip_a_mu[j]<- kappa_q_a_mu/lambda_q_a_mu[j]
      mu_q_log_a_mu[j]<- log(lambda_q_a_mu[j])- digamma(kappa_q_a_mu)

    }
    #Update of q(sigma_phi) & Update of q(a_phi)

    lambda_q_sigsq_phi <- matrix(NA, nrow = Q, ncol = L)
    mu_q_recip_sigsq_phi  <- matrix(NA, nrow = Q, ncol = L)
    mu_q_log_sigsq_phi<- matrix(NA, nrow = Q, ncol = L)

    lambda_q_a_phi <-  matrix(NA, nrow = Q, ncol = L)
    mu_q_log_a_phi <-matrix(NA, nrow = Q, ncol = L)

    for (q in 1:Q){
      for( l in 1:L){
        mu_q_u_phi <- mu_q_nu_phi[[q]][-c(1:2),l]
        Sigma_q_u_phi <- Sigma_q_nu_phi[[q]][[l]][-c(1:2), -c(1:2)]

        lambda_q_sigsq_phi[q,l] <- mu_q_recip_a_phi[q,l] + 0.5*( tr(Sigma_q_u_phi)+t(mu_q_u_phi)%*%mu_q_u_phi)
        mu_q_recip_sigsq_phi[q,l] <- kappa_q_sigsq_phi/ lambda_q_sigsq_phi[q,l]
        mu_q_log_sigsq_phi[q,l]<- log(lambda_q_sigsq_phi[q,l])-digamma(kappa_q_sigsq_phi)

        lambda_q_a_phi[q,l]<- mu_q_recip_sigsq_phi[q,l]+ 1/A^2
        mu_q_recip_a_phi[q,l]<- kappa_q_a_phi/lambda_q_a_phi[q,l]
        mu_q_log_a_phi[q,l]<- log(lambda_q_a_phi[q,l])- digamma(kappa_q_a_phi)
      }
    }

    #Update of q(omega_q) and Update q(alpha_q)
    c_1_omega <- rep(NA,Q)
    d_1_omega <- rep(NA,Q)
    alpha_1_alpha<- rep(NA,Q)
    beta_1_alpha <- rep(NA,Q)

    mu_q_log_omega <- rep(NA,Q)
    mu_q_log_1_omega <- rep(NA,Q)
    mu_q_log_alpha <- rep(0,Q)
    #mu_q_alpha <- rep(NA,Q)
    for (q in 1:Q){
      sum_term_gamma <-0
      sum_term_alpha_beta <- 0
      for (j in 1:p){
        sum_term_gamma <- sum_term_gamma + mu_q_gamma[j,q]
        sum_term_alpha_beta <- sum_term_alpha_beta+ mu_q_gamma[j,q]*(mu_q_normal_b[j,q]^2 + Sigma_q_normal_b[j,q])
      }
      c_1_omega[q]<- c_0 + sum_term_gamma
      d_1_omega[q]<- d_0+ p -sum_term_gamma

      #alpha_1_alpha[q]<- alpha_0 + 0.5* sum_term_gamma
      #beta_1_alpha[q]<- beta_0 + 0.5*sum_term_alpha_beta

      mu_q_log_omega[q]<- digamma(c_1_omega[q])- digamma(c_0+ d_0 + p)
      mu_q_log_1_omega[q]<- digamma(d_1_omega[q])-digamma(c_0+ d_0 + p)

      #mu_q_alpha[q]<- (alpha_0 + 0.5* sum_term_gamma)/(beta_0 + 0.5*sum_term_alpha_beta)
      #mu_q_log_alpha[q]<- digamma(alpha_1_alpha[q])- log ( beta_1_alpha[q] )
    }
    #Update of q(b_jq |gamma_jq)
    Sigma_q_b <- matrix(NA, nrow=p, ncol=Q)
    for (j in 1:p){
      for(q in 1:Q){
        sum_sigma <- 0
        sum_mu <- 0
        for( i in 1:N){
          mu_V_q_phi <- matrix(NA, nrow = K+2, ncol = L)
          trace_term <- diag(L)
          for (l in 1:L){
            mu_V_q_phi[,l]<- mu_q_nu_phi[[q]][,l]
            trace_term[l,l]<- tr(list_cp_C[[i]]%*% Sigma_q_nu_phi[[q]][[l]])
          }
          mu_H_q_phi <- t(mu_V_q_phi)%*%list_cp_C[[i]]%*%mu_V_q_phi+ trace_term
          sum_sigma<-sum_sigma+tr(mu_H_q_phi%*%(Sigma_q_zeta[[i]][[q]]+tcrossprod(mu_q_zeta[[q]][i,])))
          sum_mu_tilde <- rep(0, K+2 )
          for (q_tilde in 1:Q){
            if(q_tilde != q){
              mu_V_q_tilde_phi <- matrix(NA, nrow = K+2, ncol = L)
              for (l in 1:L){
                mu_V_q_tilde_phi[,l]<- mu_q_nu_phi[[q_tilde]][,l]
              }
              sum_mu_tilde<- sum_mu_tilde+ mu_q_b[j,q_tilde]*mu_V_q_tilde_phi%*%mu_q_zeta[[q_tilde]][i,]
            }
          }
          sum_mu <- sum_mu +t(mu_q_zeta[[q]][i,])%*%
            t(mu_V_q_phi)%*%(list_cp_C_Y[[i]][,j]-list_cp_C[[i]]%*%mu_q_nu_mu[[j]]-list_cp_C[[i]]%*%sum_mu_tilde)
        }
        Sigma_q_normal_b[j,q]<- 1/(mu_q_recip_sigsq_eps[j]*sum_sigma+ mu_q_alpha[q])
        mu_q_normal_b[j,q]<- Sigma_q_normal_b[j,q]*mu_q_recip_sigsq_eps[j]*sum_mu
        mu_q_gamma[j,q]<- (1 + (mu_q_recip_sigsq_eps[j]*sum_sigma+ mu_q_alpha[q])^(1/2)*
                             exp(mu_q_log_1_omega[q]-mu_q_log_omega[q]-
                                   0.5*(mu_q_normal_b[j,q]^2)*(mu_q_recip_sigsq_eps[j]*sum_sigma+ mu_q_alpha[q])-
                                   0.5*mu_q_log_alpha[q]) )^(-1)

        mu_q_b[j,q]<- mu_q_gamma[j,q]*mu_q_normal_b[j,q]
        Sigma_q_b[j,q]<- mu_q_gamma[j,q]*(Sigma_q_normal_b[j,q]+mu_q_normal_b[j,q]^2)-mu_q_b[j,q]^2
      }
    }

    #COMPUTE ELBO
    ELBO_iter <- sum(sapply(1:p, function (j) {-(sum(sapply(time_obs, function(x) length (x)))/2)*(log(2*pi)+ mu_q_log_sigsq_eps[j])-
        mu_q_recip_sigsq_eps[j]*( lambda_q_sigsq_eps[j]- mu_q_recip_a_eps[j])} ))#term of Y
    ELBO_iter <- ELBO_iter+ sum(sapply(1:p, function(j){
      E_q_inv_Sigma_nu_mu_j <- blkdiag(
        inv_Sigma_beta,
        mu_q_recip_sigsq_mu[j]*diag(K))
      log_det_j_obj <- determinant(Sigma_q_nu_mu[[j]], logarithm = TRUE)
      log_det_j <- log_det_j_obj$modulus * log_det_j_obj$sign
      -1*log(sigma_zeta)+ 0.5*log_det_j-(K/2)*(mu_q_log_sigsq_mu[j])-0.5*(unlist(t(mu_q_nu_mu[[j]])%*%E_q_inv_Sigma_nu_mu_j%*%mu_q_nu_mu[[j]])+ tr(E_q_inv_Sigma_nu_mu_j %*% Sigma_q_nu_mu[[j]])) +K/2 + 1+
        K/2*mu_q_log_sigsq_mu[j] - mu_q_recip_sigsq_mu[j]*(mu_q_recip_a_mu[j]-lambda_q_sigsq_mu[j])-
        0.5*mu_q_log_a_mu[j] -kappa_q_sigsq_mu*log(lambda_q_sigsq_mu[j])-lgamma(0.5)+lgamma(kappa_q_sigsq_mu)+
        (3/2)*mu_q_log_a_mu[j]-mu_q_recip_a_mu[j]*(1/A^2-lambda_q_a_mu[j] )-log(lambda_q_a_mu[j])-lgamma(0.5)+lgamma(1)+ 0.5*log(1/A^2)+
        (kappa_q_sigsq_eps -0.5)*mu_q_log_sigsq_eps[j]-mu_q_recip_sigsq_eps[j]*(mu_q_recip_a_eps[j]-lambda_q_sigsq_eps[j])-
        0.5*mu_q_log_a_eps[j]-kappa_q_sigsq_eps*log(lambda_q_sigsq_eps[j])-lgamma(0.5)+ lgamma(kappa_q_sigsq_eps)+
        (3/2)*mu_q_log_a_eps[j] -mu_q_recip_a_eps[j]*(1/A^2-lambda_q_a_eps[j])-
        log(lambda_q_a_eps[j])-lgamma(0.5)+lgamma(1)+ 0.5*log(1/A^2)}))
    ELBO_iter<-ELBO_iter+ sum( sapply(1: p, function(j){
      sum( sapply(1:Q, function(q){
        0.5*mu_q_gamma[j,q]*(log(Sigma_q_normal_b[j,q])+ mu_q_log_alpha[q]+1)-
          0.5*mu_q_gamma[j,q]*mu_q_alpha[q]*(Sigma_q_normal_b[j,q]+(mu_q_normal_b[j,q]^2))+mu_q_gamma[j,q]*mu_q_log_omega[q]+
          (1-mu_q_gamma[j,q])*mu_q_log_1_omega[q]-mu_q_gamma[j,q]*log(mu_q_gamma[j,q]+1e-11)- (1-mu_q_gamma[j,q])*log(1-mu_q_gamma[j,q]+1e-11)
      }))
    }))
    ELBO_iter <- ELBO_iter+ sum(sapply(1:Q, function(q){
      (c_0 -c_1_omega[q])*mu_q_log_omega[q]+ (d_0 -d_1_omega[q])*mu_q_log_1_omega[q]+ log( beta(c_1_omega[q], d_1_omega[q]))-log(beta(c_0,d_0))
    })) #term omega
    ELBO_iter<- ELBO_iter+ sum(sapply(1:Q, function(q){
      sum(sapply(1:L, function(l){
        (K/2)*mu_q_log_sigsq_phi[q,l]- (mu_q_recip_a_phi[q,l]-lambda_q_sigsq_phi[q,l])*mu_q_recip_sigsq_phi[q,l]-
          0.5* mu_q_log_a_phi[q,l]-kappa_q_sigsq_phi*log(lambda_q_sigsq_phi[q,l]) -lgamma(0.5)+lgamma(kappa_q_sigsq_phi)+
          (3/2)*(mu_q_log_a_phi[q,l])-(1/A^2-lambda_q_a_phi[q,l])*mu_q_recip_a_phi[q,l]+
          0.5*log(1/A^2)- log(lambda_q_a_phi[q,l])-lgamma(0.5)+lgamma(1)
      }))
    }))# q_phi, sigsq_phi
    ELBO_iter<- ELBO_iter+ sum(sapply(1:Q, function(q){
      sum(sapply(1:L, function(l){
        E_q_inv_Sigma_nu_phi_q_l <- blkdiag(
          inv_Sigma_beta,
          mu_q_recip_sigsq_phi[q,l]*diag(K)
        )
        log_det_q_l_obj <- determinant(Sigma_q_nu_phi[[q]][[l]], logarithm = TRUE)
        log_det_q_l <- log_det_q_l_obj$modulus * log_det_q_l_obj$sign
        -1*log(sigma_zeta)+ 0.5*log_det_q_l-(K/2)*(mu_q_log_sigsq_phi[q,l])-
          0.5*(unlist((t(mu_q_nu_phi[[q]][,l])%*%E_q_inv_Sigma_nu_phi_q_l%*%mu_q_nu_phi[[q]][,l])+ tr(E_q_inv_Sigma_nu_phi_q_l%*% Sigma_q_nu_phi[[q]][[l]])))+
          K/2 + 1 #mu_phi
      }))
    }))
    ELBO_iter<- ELBO_iter + sum( sapply(1:Q, function(q){
      sum(sapply(1:N, function (i){
        log_det_i_q_obj <- determinant(Sigma_q_zeta[[i]][[q]], logarithm = TRUE)
        log_det_i_q <- log_det_i_q_obj$modulus * log_det_i_q_obj$sign
        unlist(0.5*log_det_i_q-(L/2)*log(sigma_zeta)+L/2-(0.5/sigma_zeta)*(crossprod(mu_q_zeta[[q]][i,]) + tr(Sigma_q_zeta[[i]][[q]]))) #term of zeta

      }))
    }))
    ELBO <- c(ELBO, ELBO_iter)

    if (i_iter >1 && (ELBO[i_iter] - ELBO[i_iter-1]) < -eps){
      stop(paste0("Error : THE ELBO IS DECREASING, difference: ", ELBO[i_iter] - ELBO[i_iter-1]))
    }
    if (i_iter >1){
      message(i_iter, ":", ELBO[i_iter] - ELBO[i_iter-1])
    }
    if (i_iter >1){
      rel_converged <- (abs(ELBO[i_iter]/ELBO[i_iter-1] - 1) < tol_rel)
      abs_converged <- ((ELBO[i_iter] - ELBO[i_iter-1])/N< tol_abs)
    }
    if(i_iter >1 && (rel_converged | abs_converged)) {
      print(abs(ELBO[i_iter]/ELBO[i_iter-1] - 1))
      print(abs_converged)
      break
    }


  }
  sigsq_eps <- lambda_q_sigsq_eps/(kappa_q_sigsq_eps-1)

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

  res <- create_named_list(K, list_Y_hat, list_Y_low, list_Y_upp,
                           list_mu_hat, list_list_Phi_hat,
                           list_Zeta_hat, list_Cov_zeta_hat, list_list_zeta_ellipse,
                           Sigma_q_nu_mu,mu_q_nu_mu, Sigma_q_nu_phi,mu_q_nu_phi,mu_q_zeta, Sigma_q_zeta, sigsq_eps,
                           lambda_q_sigsq_phi, lambda_q_a_phi, lambda_q_sigsq_mu, lambda_q_a_mu, lambda_q_sigsq_eps,
                           lambda_q_a_eps, mu_q_b,Sigma_q_normal_b,mu_q_gamma,mu_q_alpha, ELBO, n_g, time_g, C_g)


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
    list_Cov_zeta_hat[[q]] <- vector("list", length = N)
    for(i in 1:N) {

      mat_transform <- list_S_inv[[q]]%*%tcrossprod(list_D_diag[[q]], list_V_orth[[q]])
      Cov_zeta_dot_hat <- tcrossprod(mat_transform%*%Sigma_q_zeta[[i]][[q]], mat_transform)
      list_Cov_zeta_hat[[q]][[i]] <- tcrossprod(scale_mat %*% Cov_zeta_dot_hat, scale_mat)
    }

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

  list_Y_hat <- list_Y_low <- list_Y_upp <- vector("list", length = p)

  list_var_vec <- lapply(1:Q, function(q) lapply(1:N, function(i)
    diag(tcrossprod(list_M_q_Phi[[q]]%*%Sigma_q_zeta[[i]][[q]], list_M_q_Phi[[q]])))) # functions assumed to be known exactly (we can use all objects pre-orthogonalisation to construct y)

  # list_sd_vec <- lapply(1:Q, function(q) lapply(1:N, function(i)
  #   sqrt(diag(tcrossprod(list_list_Phi_hat[[q]]%*%list_Cov_zeta_hat[[q]][[i]], list_list_Phi_hat[[q]]))))) # functions assumed to be known exactly

  for (j in 1:p) {

    list_Y_hat[[j]] <- list_Y_low[[j]] <- list_Y_upp[[j]] <- vector("list", length = N)
    for(i in 1:N) {

      sd_i_j <- sqrt(Reduce("+", lapply(1:Q, function(q) mu_q_b[j, q]^2 * list_var_vec[[q]][[i]])) + sigsq_eps[j])# B assumed to be known exactly (variance = b^T Psi^T Cov_zeta Psi b) so sqrt(b^2) for sd


      list_Y_hat[[j]][[i]] <- list_Y_mat[[j]][,i]
      list_Y_low[[j]][[i]] <- list_Y_mat[[j]][,i] + qnorm(0.025) * sd_i_j
      list_Y_upp[[j]][[i]] <- list_Y_mat[[j]][,i] + qnorm(0.975) * sd_i_j
    }

  }

  create_named_list(list_Y_hat, list_Y_low, list_Y_upp,
                    list_mu_hat, list_list_Phi_hat,
                    list_Zeta_hat, list_Cov_zeta_hat, list_list_zeta_ellipse,
                    mu_q_b_norm)
}
