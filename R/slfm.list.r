#' Fit SLFM to the matrices inside a directory
#'
#' This function is used to fit a Bayesian sparse
#' latent factor model to a directory of numeric matrices.
#'
#' @param path path to read the matrices from
#' @param recursive if the function should look recursively inside folders
#' @param ite number of iterations of the MCMC algorith
#' @param a prior shape parameter for Gamma distribution 
#' @param b prior scale parameter for Gamma distribution
#' @param gamma_a prior parameter for Beta distribution
#' @param gamma_b prior parameter for Beta distribution
#' @param omega prior variance of the slab component
#' @param omega_1 prior variance of the spike component
#' @param burnin burn-in size
#' @export
slfm.list <- function(
  path = ".", recursive = TRUE, ite,
  a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1,
  omega = 10, omega_1 = 0.01, burnin = 500) {

  files_list <- list.files(path, recursive = recursive)
  results_list <- lapply(files_list, function(file_name) {
    mat <- read.table(file_name)

    res <- slfm(mat, ite, a, b, gamma_a, gamma_b, omega, omega_1, burnin)
  })

}