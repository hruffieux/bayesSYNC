#' Proportion of variance explained by each latent factor
#'
#' Computes the model-based decomposition of the proportion of variance
#' explained (PVE) by each latent factor of a fitted bayesSYNC model. This is
#' the factor-level counterpart of \code{list_cumulated_pve}, which gives the
#' PVE of each FPCA component within a factor. \code{\link{bayesSYNC}} returns
#' this decomposition with the default arguments as \code{factor_pve}; this
#' function recomputes it with other arguments without refitting.
#'
#' @details
#' Under the model, the variance of variable \eqn{j} about its mean function,
#' integrated over the time domain \eqn{T}, decomposes as
#' \deqn{\int_T Var\{y_i^{(j)}(t)\} dt = \sum_q b_{jq}^2 \Lambda_q + |T| \sigma^2_{\epsilon, j},}
#' where \eqn{\Lambda_q} is the trace of the covariance operator of the latent
#' process of factor \eqn{q}. The decomposition relies on the model assumption
#' that the FPCA scores are independent across factors. The values are implied
#' by the fitted model; they are not an \eqn{R^2} against the observed data.
#'
#' The structured variance of factor \eqn{q} is estimated by
#' \eqn{V_q = S_q \Lambda_q}. \eqn{S_q = \sum_j E(b_{jq}^2)} uses the posterior
#' second moments of the loadings rather than the squared posterior means,
#' which would understate the contribution of factors supported by many weakly
#' included variables. \eqn{\Lambda_q} is the sum of the eigenvalues of factor
#' \eqn{q} over the \code{L} fitted components, so it is truncated at \code{L};
#' the within-factor PVE in \code{list_cumulated_pve} is normalised by the same
#' truncated sum. The noise variance is
#' \eqn{V_{noise} = |T| \sum_j \sigma^2_{\epsilon, j}}, summed over all
#' variables whatever the value of \code{factors} and \code{loadings}, with
#' \eqn{|T|} the length of the domain spanned by \code{time_g} and
#' \eqn{\sigma^2_{\epsilon, j}} the posterior mean of the noise variance.
#'
#' Taking \eqn{\Lambda_q} as the sum of the eigenvalues assumes that the
#' estimated eigenfunctions are orthonormal in \eqn{L^2(T)}. The post-processing
#' makes them orthogonal in the Euclidean inner product on \code{time_g} and
#' gives them unit norm under the trapezoidal rule. On an evenly spaced grid
#' (the default) this approximates \eqn{L^2(T)} orthonormality up to an error
#' of order \code{1 / n_g}. On an unevenly spaced \code{time_g}, the eigenvalues,
#' and hence both this decomposition and \code{list_cumulated_pve}, may not be
#' on the \eqn{L^2(T)} scale.
#'
#' The model does not identify the scale of the loadings and of the latent
#' processes separately: \eqn{b_q \to c b_q} with \eqn{h^{(q)} \to h^{(q)} / c}
#' leaves the likelihood unchanged. Only the product \eqn{S_q \Lambda_q} is
#' identified, which is why neither term is returned on its own.
#'
#' If the variables were not standardised before fitting
#' (\code{bool_scale = FALSE}), \eqn{S_q} is dominated by the variables with
#' the largest variance, and the PVE reflects that.
#'
#' Only the factors in \code{factors_included} enter the decomposition, in both
#' the numerator and the denominator. \code{total} is therefore the share of the
#' structured variance of the retained factors plus the noise variance, not of
#' the whole fitted variance. Factors that are left out get \code{NA}.
#'
#' Unlike the components within a factor, the factors are not ordered by the
#' model. The factors, scores and loadings of the fitted object are not
#' reordered; the ordering by PVE is given separately in \code{order}.
#'
#' @param fit Object returned by \code{\link{bayesSYNC}}.
#' @param factors Factors to include in the decomposition. Either
#'        \code{"mpm"} (default) for the factors retained under the median
#'        probability model on the posterior factor inclusion probabilities
#'        (\code{factor_ppi > 0.5}), \code{"all"} for all factors, or a vector
#'        of factor indices.
#' @param loadings Either \code{"posterior"} (default), which uses the full
#'        posterior second moments of the loadings,
#'        \eqn{\gamma_{jq}(m_{jq}^2 + s^2_{jq})}, or \code{"mpm"}, which uses the
#'        second moments conditional on inclusion,
#'        \eqn{m_{jq}^2 + s^2_{jq}}, for the loadings with
#'        \eqn{\gamma_{jq} > 0.5} and zero for the others. The \code{"mpm"}
#'        option gives the PVE of the selected module rather than of the whole
#'        factor, and the two can differ noticeably.
#'
#' @return A list with elements:
#' \describe{
#'   \item{total}{Vector of length \code{Q}: share of the variance of the
#'         retained factors plus noise attributable to each factor. These
#'         proportions do not sum to one; the shortfall is the measurement
#'         error.}
#'   \item{structured}{Vector of length \code{Q}: share of the structured
#'         variance of the retained factors attributable to each factor. These
#'         proportions sum to one over the retained factors.}
#'   \item{order}{Names of the retained factors by decreasing PVE.}
#'   \item{factors_included}{Indices of the factors that entered the
#'         decomposition, named after the columns of \code{B_hat}.}
#'   \item{options}{List recording the values of \code{factors} and
#'         \code{loadings} used.}
#' }
#' The vectors are proportions, not percentages as in
#' \code{list_cumulated_pve}, and their entries line up with the columns of
#' \code{B_hat}.
#'
#' @seealso \code{\link{bayesSYNC}}
#'
#' @export
#'
get_factor_pve <- function(fit, factors = "mpm", loadings = "posterior") {

  required_fields <- c("B_slab_mean", "B_slab_var", "ppi", "factor_ppi",
                       "list_eigenvalues", "sigsq_eps", "time_g")
  missing_fields <- setdiff(required_fields, names(fit))
  if (length(missing_fields) > 0) {
    stop(paste0("fit does not contain ", paste(missing_fields, collapse = ", "),
                ". It may have been produced by an earlier version of ",
                "bayesSYNC; refit the model to use get_factor_pve."))
  }

  compute_factor_pve(B_slab_mean = fit$B_slab_mean,
                     B_slab_var = fit$B_slab_var,
                     ppi = fit$ppi,
                     factor_ppi = fit$factor_ppi,
                     list_eigenvalues = fit$list_eigenvalues,
                     sigsq_eps = fit$sigsq_eps,
                     time_g = fit$time_g,
                     factors = factors,
                     loadings = loadings)
}


compute_factor_pve <- function(B_slab_mean, B_slab_var, ppi, factor_ppi,
                               list_eigenvalues, sigsq_eps, time_g,
                               factors = "mpm", loadings = "posterior") {

  Q <- ncol(ppi)
  factor_names <- colnames(ppi)

  check_structure(loadings, "vector", "string", 1)
  if (!(loadings %in% c("posterior", "mpm"))) {
    stop("loadings must be either \"posterior\" or \"mpm\".")
  }

  if (is.character(factors)) {
    check_structure(factors, "vector", "string", 1)
    if (factors == "mpm") {
      factors_included <- which(factor_ppi > 0.5)
    } else if (factors == "all") {
      factors_included <- 1:Q
    } else {
      stop("factors must be \"mpm\", \"all\" or a vector of factor indices.")
    }
  } else {
    check_structure(factors, "vector", "numeric")
    check_natural(factors)
    if (any(factors > Q) || anyDuplicated(factors)) {
      stop(paste0("factors must contain distinct indices between 1 and ", Q, "."))
    }
    factors_included <- sort(as.integer(factors))
  }
  factors_included <- unname(factors_included)
  names(factors_included) <- factor_names[factors_included]

  # Posterior second moments of the loadings. Under the spike-and-slab
  # variational family, E(b_jq^2) = gamma_jq * (m_jq^2 + s2_jq). The "mpm"
  # option keeps the second moment conditional on inclusion for the entries
  # with gamma_jq > 0.5 and sets the others to zero.
  slab_second_moment <- B_slab_mean^2 + B_slab_var
  if (loadings == "posterior") {
    B2 <- ppi * slab_second_moment
  } else {
    B2 <- (ppi > 0.5) * slab_second_moment
  }
  S <- colSums(B2)

  # Trace of the covariance operator of each latent process, truncated to the
  # L fitted components. This equals the sum of the eigenvalues because the
  # eigenfunctions are (approximately) orthonormal in L2 over time_g.
  Lambda <- sapply(list_eigenvalues, sum)

  # Only the product S_q * Lambda_q is identified (b_q -> c * b_q and
  # h^(q) -> h^(q) / c leave the likelihood unchanged), so neither term is
  # returned on its own.
  V <- unname(S * Lambda)

  # |T| is the length of the domain over which the eigenfunctions are
  # normalised. It does not cancel, as it only appears in the noise term.
  V_noise <- diff(range(time_g)) * sum(sigsq_eps)

  total <- structured <- rep(NA_real_, Q)
  names(total) <- names(structured) <- factor_names

  if (length(factors_included) == 0) {
    warning(paste0("No factor was included in the factor-level PVE ",
                   "decomposition (with factors = \"mpm\", no factor has ",
                   "factor_ppi > 0.5); all values are NA."))
    pve_order <- character(0)
  } else {
    V_struct <- sum(V[factors_included])
    if (V_struct <= 0) {
      warning(paste0("The structured variance of the included factors is zero ",
                     "(e.g., with loadings = \"mpm\", no loading has ppi > 0.5 ",
                     "in these factors); the structured PVE is undefined."))
    }
    total[factors_included] <- V[factors_included] / (V_struct + V_noise)
    structured[factors_included] <- V[factors_included] / V_struct
    pve_order <- names(sort(total[factors_included], decreasing = TRUE))
  }

  options <- list(factors = factors, loadings = loadings)

  list(total = total,
       structured = structured,
       order = pve_order,
       factors_included = factors_included,
       options = options)
}
