#' Fit SLFM to the matrices inside a directory
#'
#' This function is used to fit a Bayesian sparse
#' latent factor model to a directory of numeric matrices.
#'
#' @param path path to read the matrices from
#' @param recursive if the function should look recursively inside folders
#' @param ite number of iterations of the MCMC algorithm
#' @param a prior shape parameter for Gamma distribution 
#' @param b prior scale parameter for Gamma distribution
#' @param gamma_a prior parameter for Beta distribution
#' @param gamma_b prior parameter for Beta distribution
#' @param omega prior variance of the slab component
#' @param omega_1 prior variance of the spike component
#' @param burnin burn-in size
#' @importFrom coda HPDinterval
#' @importFrom tools file_path_sans_ext
#' @export
slfm_list <- function(
  path = ".", recursive = TRUE, ite,
  a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1,
  omega = 10, omega_1 = 0.01, burnin = 500) {

  files_list <- list.files(path, recursive = recursive, full.names = T)

  message("* |Press Ctrl + C to cancel...")
  pb <- txtProgressBar(0, length(files_list), style = 3, title = "BLA")

  results_list <- list()
  for(i in 1:length(files_list)) {
    file_name <- files_list[i]
    mat <- read.table(file_name)

    res <- slfm(mat, ite, a, b, gamma_a, gamma_b, omega, omega_1, burnin)
    hpd <- coda::HPDinterval(res$p_star)
    clas <- class_interval(hpd)
    final_clas <- names(which.max(table(clas)))
    results_list[[i]] <- c(name=basename(tools::file_path_sans_ext(file_name)), clas=final_clas)

    setTxtProgressBar(pb, i)
  }
  close(pb)
  ret <- as.data.frame(do.call(rbind, results_list))
  names(ret) <- c("File", "Classification")
  lvls <- levels(ret[,2])
  lvls[lvls == "P"] <- "Present"
  lvls[lvls == "M"] <- "Marginal"
  lvls[lvls == "A"] <- "Absent"
  levels(ret[,2]) <- lvls
  ret
}