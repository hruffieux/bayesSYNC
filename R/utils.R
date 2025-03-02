create_named_list <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}

all_same <- function(x) length(unique(x)) == 1

check_natural <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps | abs(x - round(x)) > eps)) {
    stop(paste0(deparse(substitute(x)),
                " must be natural."))
  }
}

check_positive <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps)) {
    err_mess <- paste0(deparse(substitute(x)), " must be positive (larger than precision zero in R).")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_zero_one <- function(x){
  if (any(x < 0) | any(x > 1)) {
    err_mess <- paste0(deparse(substitute(x)), " must lie between 0 and 1.")
    if (length(x) > 1) err_mess <- paste0("All entries of ", err_mess)
    stop(err_mess)
  }
}

check_structure <- function(x, struct, type, size = NULL,
                            null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "a double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "an integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "a numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "a boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }

  bool_size <- TRUE # for case size = NULL (no assertion on the size/dimension)
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0) # not an empty vector
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste0(" of length ", paste0(size, collapse=" or "))
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0) # not an empty matrix
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste0(" of dimension ", size[1], " x ", size[2])
    }
  }

  correct_obj <- bool_struct & bool_type & bool_size

  bool_null <- is.null(x)

  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }

    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }

  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }

  if(!(correct_obj)) {
    stop(paste0(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
                size_mess, inf_mess, na_mess, null_mess, "."))
  }
}


#' Create a temporal-grid objects for use in FPCA algorithms.
#'
#' This function is used to create a time grid of desired density and design
#' matrices based on an observed time grid object.
#'
#' @param time_obs Vector or list of vectors containing time of observations for
#'                 univariate or multivariate curves, respectively.
#' @param K Number of O'Sulivan spline functions to be used in the FPCA
#'          algorithms. If set to \code{NULL} will be set according to the rule
#'          of Ruppert (2002), also enforcing K >=7.
#' @param n_g Desired size for dense grid.
#' @param time_g Dense grid provided as a vector of size \code{n_g}. If provided,
#'               then \code{n_g} must be \code{NULL} as will be taken to be
#'               \code{length(time_g)}.
#' @param int_knots Position of interior knots. Default is \code{NULL} for
#'                  evenly placed.
#' @param format_univ Boolean indicating whether the univariate format is used
#'                    in case of univariate curves.
#' @return C Design matrix C(t) constructed from the set of K spline functions
#'           based on the observation time grid.
#' @return n_g Number of time points in the dense grid.
#' @return time_g Dense time grid constructed.
#' @return C_g Design matrix C(t) constructed from the set of K spline functions
#'             based on the dense time grid.
#'
#' @export
#'
get_grid_objects <- function(time_obs, K, n_g = 1000, time_g = NULL,
                             int_knots = NULL,
                             format_univ = FALSE) {


  if (is.null(int_knots)) {
    if(format_univ) {

      if (is.null(K)) {  # Ruppert (2002) sets a simple default value for K as min(nobs/4,40), where nobs is the number of observations.
        # here since nobs differs for each i, we take the nobs / 4 = round(median(obs_i)/4)
        # and we enforce that K>=7
        K <- max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i))/4), 40)), 7)
      }

      unique_time_obs <- sort(unique(Reduce(c, time_obs)))
      int_knots <- quantile(unique_time_obs, seq(0, 1, length=K)[-c(1,K)])
    } else {

      p <- length(time_obs[[1]])
      if (is.null(K)) {   # Ruppert (2002) sets a simple default value for K as min(nobs/4,40), where nobs is the number of observations.
        # here since nobs differs for each i, we take the nobs / 4 = round(median(obs_i)/4)
        # and we enforce that K>=7
        K <- sapply(1:p, function(j) max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i[[j]]))/4), 40)), 7))
      } else {
        check_structure(K, "vector", "numeric", c(1, p))
        if (length(K)==1) K <- rep(K, p)
      }

      unique_time_obs <- unname(sort(unlist(time_obs)))
      int_knots <- lapply(
        K,
        function(x) quantile(unique_time_obs, seq(0, 1, length = x)[-c(1, x)])
      )
    }
  }

  N <- length(time_obs)

  if (is.null(time_g)) {
    stopifnot(!is.null(n_g))
    time_g <- seq(0, 1, length.out = n_g)
  } else {
    if(is.null(n_g)) n_g <- length(time_g)
  }

  C <- vector("list", length=N)

  if (format_univ) {
    for(i in 1:N) {
      X <- X_design(time_obs[[i]])
      Z <- ZOSull(time_obs[[i]], range.x=c(0, 1), intKnots=int_knots)
      C[[i]] <- cbind(X, Z)
    }

    X_g <- X_design(time_g)
    Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
    C_g <- cbind(X_g, Z_g)
  } else {

    p <- length(time_obs[[1]])

    for(i in 1:N) {

      C[[i]] <- vector("list", length = p)
      for(j in 1:p) {

        X <- X_design(time_obs[[i]][[j]])
        Z <- ZOSull(time_obs[[i]][[j]], range.x = c(0, 1), intKnots = int_knots[[j]])
        C[[i]][[j]] <- cbind(X, Z)
      }
    }

    X_g <- X_design(time_g)
    Z_g <- lapply(
      int_knots,
      function(x) ZOSull(time_g, range.x = c(0, 1), intKnots = x)
    )
    C_g <- lapply(Z_g, function(Z) cbind(X_g, Z))
  }

  create_named_list(C, n_g, time_g, C_g)
}


X_design <- function(x) {

  if(is.list(x)) {

    x <- do.call(cbind, x)
  }

  X <- cbind(1, x)
  return(X)
}


is_int <- function(x, tol = .Machine$double.eps^0.5) {

  abs(x - round(x)) < tol
}


tr <- function(X) {

  if(nrow(X)!=ncol(X)) stop("X must be a square matrix.")

  ans <- sum(diag(X))
  return(ans)
}


cprod <- function(x, y) {

  if(missing(y)) {

    if(!is.vector(x)) {

      stop("Use the crossprod function for matrix inner products")
    }

    y <- x
  }

  if(!is.vector(y) & !is.vector(x)) {

    stop("Use the crossprod function for matrix inner products")
  }

  ans <- as.vector(crossprod(x, y))
  return(ans)
}


normalise <- function(x) {

  ans <- x/sqrt(cprod(x))
  return(ans)
}


E_cprod <- function(mean_1, Cov_21, mean_2, A) {

  if(missing(A)) {

    A <- diag(length(mean_1))
  }

  tr_term <- tr(Cov_21 %*% A)
  cprod_term <- cprod(mean_1, A %*% mean_2)
  ans <- tr_term + cprod_term
  return(ans)
}


E_h <- function(L, mean_1, Cov_21, mean_2, A) {

  if(missing(A)) {

    A <- diag(length(mean_1))
  }

  d_1 <- length(mean_1)
  inds_1 <- matrix(1:d_1, ncol = L)

  ans <- rep(NA, L)
  for(l in 1:L) {

    mean_1_l <- mean_1[inds_1[, l]]
    Cov_21_l <- Cov_21[, inds_1[, l]]
    ans[l] <- E_cprod(mean_1_l, Cov_21_l, mean_2, A)
  }

  return(ans)
}


E_H <- function(L_1, L_2, mean_1, Cov_21, mean_2, A) {

  if(missing(A)) {

    A <- diag(length(mean_1))
  }

  d_1 <- length(mean_1)
  d_2 <- length(mean_2)
  L <- L_1 + L_2

  inds_1 <- matrix(1:d_1, ncol = L_1)
  inds_2 <- matrix(1:d_2, ncol = L_2)

  ans <- matrix(NA, L_1, L_2)
  for(l in 1:L_1) {

    mean_1_l <- mean_1[inds_1[, l]]
    for(k in 1:L_2) {

      mean_2_k <- mean_2[inds_2[, k]]

      Cov_21_kl <- Cov_21[inds_2[, k], inds_1[, l]]

      ans[l, k] <- E_cprod(mean_1_l, Cov_21_kl, mean_2_k, A)
    }
  }

  return(ans)
}


vec <- function(A) {
  return(as.vector(A))
}


vecInverse <- function(a) {
  is.wholenumber <- function(x,tol=sqrt(.Machine$double.eps))
    return(abs(x-round(x))<tol)

  a <- as.vector(a)
  if (!is.wholenumber(sqrt(length(a))))
    stop("input vector must be a perfect square in length")

  dmnVal <- round(sqrt(length(a)))

  A <- matrix(NA,dmnVal,dmnVal)
  for (j in 1:dmnVal)
    A[,j] <- a[((j-1)*dmnVal+1):(j*dmnVal)]

  return(A)
}

#' Function performing trapesoidal integration.
#'
#' @param xgrid Grid.
#' @param fgrid Function on the grid.
#' @return Integration result.
#'
#' @export
#'
trapint <- function(xgrid,fgrid) {
  ng <- length(xgrid)

  xvec <- xgrid[2:ng] - xgrid[1:(ng-1)]
  fvec <- fgrid[1:(ng-1)] + fgrid[2:ng]

  integ <- sum(xvec*fvec)/2

  return(integ)
}

wait <- function() {
  cat("Hit Enter to continue\n")
  ans <- readline()
  invisible()
}


#' Flip the sign of eigenfunctions and corresponding scores
#'
#' This function is used to alter the sign of specified eigenfunctions and
#' scores. Choosing one eigenfunction over its opposite sign has no effect on
#' the resulting fits, although one choice of sign may provide more natural
#' interpretation of the eigenfunction.
#'
#' @param vec_flip Vector of size L whose entry l is -1 if sign needs to be
#'                 flipped for eigenfunction l = 1, ..., L, or 1 if sign is keep
#'                 the same.
#' @param list_Psi_hat List or matrix of eigenfunctions as returned by
#'                     \code{\link{bayesSYNC}}.
#' @param Zeta_hat Matrix of estimated scores as returned by
#'                     \code{\link{bayesSYNC}}.
#' @param zeta_ellipse 95\% posterior credible boundaries for the scores as
#'                     return by \code{\link{bayesSYNC}}.
#'
#' @return An object containing the eigenfunctions, scores and credible
#'         boundaries with altered sign.
#'
#' @export
#'
flip_sign <- function(vec_flip, list_Psi_hat, Zeta_hat, zeta_ellipse = NULL) {

  stopifnot(all(vec_flip == 1 | vec_flip == -1))

  list_Psi_hat[,seq_along(vec_flip)] <- sapply(seq_along(vec_flip), function(ll) {
    list_Psi_hat[,ll] * vec_flip[ll]
  })

  Zeta_hat[,seq_along(vec_flip)] <- sweep(Zeta_hat[,seq_along(vec_flip)], 2, vec_flip, "*")

  if (!is.null(zeta_ellipse)) {
    zeta_ellipse <- lapply(zeta_ellipse, function(zeta_ellipse_subj) sweep(zeta_ellipse_subj, 2, vec_flip, "*"))
  }

  create_named_list(list_Psi_hat, Zeta_hat, zeta_ellipse)
}

# Function to compute the Frobenius norm between two matrices
frobenius_norm <- function(A, B) {
  sqrt(sum((A - B)^2))
}


#' @export
match_factor_and_sign <- function(B, B_hat, ppi, factor_ppi, Zeta, list_Zeta_hat,
                                  list_list_Phi_hat, list_cumulated_pve,
                                  list_Cov_zeta_hat) {

  perm_factor <- match_factors(B, B_hat, ppi, factor_ppi, list_Zeta_hat,
                               list_list_Phi_hat, list_cumulated_pve,
                               list_Cov_zeta_hat)

  best_perm <- perm_factor$best_perm
  perm_sign_factor <- perm_factor$perm_sign_factor
  perm_list_cumulated_pve <- perm_factor$perm_list_cumulated_pve
  perm_list_Cov_zeta_hat <- perm_factor$perm_list_Cov_zeta_hat

  perm_B_hat <- perm_factor$perm_B_hat
  perm_ppi <- perm_factor$perm_ppi
  perm_factor_ppi <- perm_factor$perm_factor_ppi
  perm_list_Zeta_hat <- perm_factor$perm_list_Zeta_hat
  perm_list_list_Phi_hat <- perm_factor$perm_list_list_Phi_hat

  perm_B_hat_untrimmed <- perm_factor$perm_B_hat_untrimmed
  perm_ppi_untrimmed <- perm_factor$perm_ppi_untrimmed
  perm_list_Zeta_hat_untrimmed <- perm_factor$perm_list_Zeta_hat_untrimmed
  perm_list_list_Phi_hat_untrimmed <- perm_factor$perm_list_list_Phi_hat_untrimmed

  perm_sign <- match_sign_components(Zeta, perm_list_Zeta_hat, perm_list_list_Phi_hat)

  perm_sign_fpca <- perm_factor$perm_sign_fpca

  perm_list_Zeta_hat <- perm_sign$perm_list_Zeta_hat
  perm_list_list_Phi_hat <- perm_sign$perm_list_list_Phi_hat

  bool_rescale_loadings_and_scores <- T # so they match the simulated loadings and scores, respectively
  if (bool_rescale_loadings_and_scores) {
    # B is identifiable up to multiplicative sign on its columns
    norm_col_B <- sqrt(colSums(B^2)) # this is unknown in practice.
    norm_col_B_hat <- sqrt(colSums(perm_B_hat^2))

    Q_true <- ncol(B)
    for (q in 1:Q_true) {
      perm_B_hat[,q] <- perm_B_hat_untrimmed[,q] <- perm_B_hat[,q] * norm_col_B[q] / norm_col_B_hat[q]
      perm_list_Zeta_hat[[q]] <- perm_list_Zeta_hat_untrimmed[[q]] <- perm_list_Zeta_hat[[q]] * norm_col_B_hat[q] / norm_col_B[q]
      perm_list_Cov_zeta_hat[[q]] <- lapply(perm_list_Cov_zeta_hat[[q]], function(perm_list_Cov_zeta_hat_q_i)  perm_list_Cov_zeta_hat_q_i * norm_col_B_hat[q]^2 / norm_col_B[q]^2)# check
    }

  }

  create_named_list(best_perm, perm_sign_factor, perm_sign_fpca, perm_factor_ppi,
                    perm_list_cumulated_pve, perm_list_Cov_zeta_hat,
                    perm_B_hat, perm_ppi, perm_list_Zeta_hat, perm_list_list_Phi_hat,
                    perm_B_hat_untrimmed, perm_ppi_untrimmed,
                    perm_list_Zeta_hat_untrimmed, perm_list_list_Phi_hat_untrimmed)

}

#
match_factors <- function(B, B_hat, ppi, factor_ppi, list_Zeta_hat,
                          list_list_Phi_hat, list_cumulated_pve, list_Cov_zeta_hat) {

  Q_true <- ncol(B)
  Q <- ncol(B_hat)

  if (Q_true > Q) {
    stop("The number of estimated factors, Q, must be larger than the number of true factors.")
  } else if (Q_true < Q) {
    warning("Dropping superfluous factors based on the true loadings.")
  }
  true_pat <- 1*(abs(B) > 0)

  # Generate all possible permutations of column indices
  perms <- gtools::permutations(Q, Q_true)

  # Initialize minimum distance and best permutation
  min_distance <- Inf
  best_perm <- NULL

  # Iterate over all permutations to find the best one
  for (i in 1:nrow(perms)) {
    permuted_ppi <- ppi[, perms[i, ], drop = F]
    distance <- frobenius_norm(true_pat, permuted_ppi)

    if (distance < min_distance) {
      min_distance <- distance
      best_perm <- perms[i, ]
    }
  }

  # Reorder B with the best permutation found
  perm_ppi <- ppi[, best_perm, drop = F]
  perm_B_hat <- B_hat[, best_perm, drop = F]

  perm_sign_factor <- rep(NA, Q_true)
  perm_list_Zeta_hat <- perm_list_list_Phi_hat <- vector("list", Q_true)
  perm_list_Zeta_hat_untrimmed <- perm_list_list_Phi_hat_untrimmed <- vector("list", Q)
  for (q in 1:Q) {

    if (q <= Q_true) {
      perm_sign_factor[q] <- ifelse(frobenius_norm(B[,q], perm_B_hat[,q]) > frobenius_norm(B[,q], -perm_B_hat[,q]), -1, 1)
      perm_B_hat[,q] <- perm_sign_factor[q]*perm_B_hat[,q]

      perm_list_Zeta_hat[[q]] <- perm_list_Zeta_hat_untrimmed[[q]] <- perm_sign_factor[q]*list_Zeta_hat[[best_perm[q]]]
      perm_list_list_Phi_hat[[q]] <- perm_list_list_Phi_hat_untrimmed[[q]] <- list_list_Phi_hat[[best_perm[q]]]
    } else {
      perm_list_Zeta_hat_untrimmed[[q]] <- list_Zeta_hat[[setdiff(1:Q, best_perm)[q-Q_true]]]
      perm_list_list_Phi_hat_untrimmed[[q]] <- list_list_Phi_hat[[setdiff(1:Q, best_perm)[q-Q_true]]]
    }
  }

  perm_ppi_untrimmed <- cbind(perm_ppi, ppi[, -best_perm, drop = F])
  perm_B_hat_untrimmed <- cbind(perm_B_hat, B_hat[, -best_perm, drop = F])
  perm_factor_ppi <- c(factor_ppi[best_perm], factor_ppi[-best_perm]) # posterior probabilities of inclusion of the permuted factors
  perm_list_cumulated_pve <- c(list_cumulated_pve[best_perm], list_cumulated_pve[-best_perm])
  perm_list_Cov_zeta_hat <- c(list_Cov_zeta_hat[best_perm], list_Cov_zeta_hat[-best_perm])

  create_named_list(best_perm, perm_sign_factor, perm_factor_ppi,
                    perm_B_hat_untrimmed, perm_ppi_untrimmed,
                    perm_B_hat, perm_ppi,
                    perm_list_Zeta_hat, perm_list_list_Phi_hat,
                    perm_list_Zeta_hat_untrimmed, perm_list_list_Phi_hat_untrimmed,
                    perm_list_cumulated_pve, perm_list_Cov_zeta_hat)
}


# we also potentially need to swap the signs of zeta and phi for a given component
#
match_sign_components <- function(Zeta,
                                  list_Zeta_hat, # trimmed versions of the estimate, i.e., after discarding the irrelevant factors using match_factors
                                  list_list_Phi_hat){

  Q_true <- length(Zeta)
  N <- nrow(Zeta[[1]])

  L_true <- ncol(Zeta[[1]])
  L <- ncol(list_Zeta_hat[[1]])

  if (L_true > L) {
    stop("The number of estimated components, L, must be larger than the number of true components.")
  } else if (L_true < L) {
    warning("Number of estimated components L is greater than true number of components L_true. \n Swapping sign of the first L_true components only.")
  }

  perm_list_Zeta_hat <- list_Zeta_hat
  perm_list_list_Phi_hat <- list_list_Phi_hat

  perm_matrix <- matrix(0, nrow= Q_true, ncol = L_true)
  perm_sign_fpca <- matrix(0, nrow = Q_true, ncol = L_true)

  for (q in 1:Q_true){
    for (l in 1:L_true){
      Zeta_lq <- sapply(1:N, function(i) Zeta[[q]][i, l])
      #Zeta_hat_orth <- sapply(1:nrow(perm_list_Zeta_hat[[q]]),function(i) perm_list_Zeta_hat[[q]][i,l])
      corr_zeta <- sapply(1:L_true, function(l_tilde)
        sapply (1:Q_true, function(q_tilde)
          cor(Zeta_lq, sapply(1:N, function(i) list_Zeta_hat[[q_tilde]][i,l_tilde]))))

      corr_zeta[is.na(corr_zeta)] <- 0
      if (Q_true > 1) {
        q_tilde <- which(abs(corr_zeta) == max(abs(corr_zeta)), arr.ind = TRUE)[1,1]
        l_tilde <- which(abs(corr_zeta) == max(abs(corr_zeta)), arr.ind = TRUE)[1,2]
        perm_sign_fpca[q, l] <- sign(corr_zeta[q_tilde, l_tilde])
      } else {
        q_tilde <- 1
        l_tilde <- which(abs(corr_zeta) == max(abs(corr_zeta)))
        perm_sign_fpca[q, l] <- sign(corr_zeta[l_tilde])
      }

      perm_list_Zeta_hat[[q]][, l] <- perm_sign_fpca[q,l]*list_Zeta_hat[[q_tilde]][, l_tilde]
      perm_list_list_Phi_hat[[q]][, l]<- perm_sign_fpca[q,l]*list_list_Phi_hat[[q_tilde]][, l_tilde]
    }
  }
  res <- create_named_list(perm_sign_fpca, perm_list_Zeta_hat, perm_list_list_Phi_hat)
}


log_one_plus_exp_ <- function(x) { # computes log(1 + exp(x)) avoiding
  # numerical overflow
  m <- x
  m[x < 0] <- 0

  log(exp(x - m) + exp(- m)) + m
}

log1pExp <- function(x) {
  ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}


get_annealing_ladder_ <- function(anneal, verbose) {

  # ladder set following:
  # Importance Tempering, Robert B. Gramacy & Richard J. Samworth, pp.9-10, arxiv v4

  k_m <- 1 / anneal[2]
  m <- anneal[3]

  if(anneal[1] == 1) {

    type <- "geometric"

    delta_k <- k_m^(1 / (1 - m)) - 1

    ladder <- (1 + delta_k)^(1 - m:1)

  } else if (anneal[1] == 2) { # harmonic spacing

    type <- "harmonic"

    delta_k <- ( 1 / k_m - 1) / (m - 1)

    ladder <- 1 / (1 + delta_k * (m:1 - 1))

  } else { # linear spacing

    type <- "linear"

    delta_k <- (1 - k_m) / (m - 1)

    ladder <- k_m + delta_k * (1:m - 1)
  }

  if (verbose != 0)
    cat(paste0("** Annealing with ", type," spacing ** \n\n"))

  ladder

}

# Internal function implementing sanity checks for the annealing schedule
# specification.
#
check_annealing <- function(anneal, verbose) {

  check_structure(anneal, "vector", "numeric", 3, null_ok = TRUE)

  if (!is.null(anneal)) {

    if (verbose) cat("== Checking the annealing schedule ... \n\n")

    check_natural(anneal[c(1, 3)])
    check_positive(anneal[2])

    if (!(anneal[1] %in% 1:3))
      stop(paste0("The annealing spacing scheme must be set to 1 for geometric ",
                  "2 for harmonic or 3 for linear spacing."))

    if (anneal[2] >= 2)
      stop(paste0("Initial annealing temperature must be strictly smaller than 2.\n ",
                  "Please decrease it."))

    if (anneal[3] > 1000)
      stop(paste0("Temperature grid size very large. This may be unnecessarily ",
                  "computationally demanding. Please decrease it."))

    if (verbose) cat("... done. == \n\n")


  }

}

