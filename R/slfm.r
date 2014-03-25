#' Sparse Latent Factor Model
#'
#' This function is used to fit a Bayesian sparse
#' latent factor model.
#'
#' @param x matrix with the pre-processed data
#' @param ite number of iterations of the MCMC algorith
#' @param a prior shape parameter for Gamma distribution 
#' @param b prior scale parameter for Gamma distribution
#' @param gamma_a prior parameter for Beta distribution
#' @param gamma_b prior parameter for Beta distribution
#' @param omega prior variance of the slab component
#' @param omega_1 prior variance of the spike component
#' @param burnin burn-in size
#' @export
#' @importFrom coda as.mcmc
#' @importFrom Rcpp evalCpp
#' @useDynLib slfm
slfm <- function(
  x, ite, a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1,
  omega = 10, omega_1 = 0.01, burnin = 500) {
  
  # Convert the x input to numeric matrix
  x <- data.matrix(x)

  res <- gibbs(x, ite, a, b, gamma_a, gamma_b, omega, omega_1)
  coda::as.mcmc(res[burnin:ite,])
}