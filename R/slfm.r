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
  omega = 10, omega_1 = 0.01, burnin = round(0.25*ite)) {
  
  # Convert the x input to numeric matrix
  x <- data.matrix(x)

  res <- gibbs(x, ite, a, b, gamma_a, gamma_b, omega, omega_1)
  p_star_matrix <- coda::as.mcmc(res[["p_star"]][burnin:ite,])
  alpha_matrix <- coda::as.mcmc(res[["alpha"]][burnin:ite,])
  lambda_matrix <- coda::as.mcmc(res[["lambda"]][burnin:ite,])
  sigma_matrix <- coda::as.mcmc(res[["sigma"]][burnin:ite,])

  stats_alpha <- summary(alpha_matrix)$statistics
  hpds_alpha <- HPDinterval(alpha_matrix)
  table_alpha <- cbind(stats_alpha, hpds_alpha)[,-4]
  colnames(table_alpha)[4:5] = c("Upper HPD", "Lower HPD")

  stats_lambda <- summary(lambda_matrix)$statistics
  hpds_lambda <- HPDinterval(lambda_matrix)
  table_lambda <- cbind(stats_lambda, hpds_lambda)[,-4]
  colnames(table_lambda)[4:5] = c("Upper HPD", "Lower HPD")

  stats_sigma <- summary(sigma_matrix)$statistics
  hpds_sigma <- HPDinterval(sigma_matrix)
  table_sigma <- cbind(stats_sigma, hpds_sigma)[,-4]
  colnames(table_sigma)[4:5] = c("Upper HPD", "Lower HPD")

  hpds_p_star <- HPDinterval(p_star_matrix)
  alpha_clas <- format_classification(class_interval(hpds_p_star))

  obj <- list(
    x = x,
    p_star = p_star_matrix,
    alpha = table_alpha,
    lambda = table_lambda,
    sigma = table_sigma,
    classification = alpha_clas)
  class(obj) <- "slfm"
  obj
}

print.slfm <- function(x) {
  cat("SLFM object", "\n")
  cat("\n")
  cat("Dimensions","\n")
  cat("- alpha:", nrow(x$alpha),"\n")
  cat("- lambda:", nrow(x$lambda),"\n")
  cat("\n")
  cat("Classification:","\n")
  print(x$classification)
}

plot.slfm <- function(x) {
  plot_slfm(x$x, as.matrix(x$p_star), x$classification)
}