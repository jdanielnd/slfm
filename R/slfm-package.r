#' slfm: the sparse latent factor model package for R.
#'
#' Set of tools to find coherent patterns in microarray data
#' using a Bayesian Sparse Latent Factor Model - SLFM (Duarte and Mayrink 2015 -
#' http://link.springer.com/chapter/10.1007%2F978-3-319-12454-4_15).
#' Considerable effort has been put into making slfm fast and memory efficient,
#' turning it an interesting alternative to simpler methods in terms
#' of execution time. It implements versions of the SLFM based on two type
#' of mixtures priors for the loadings: one relying on a degenerate component at zero and the other using
#' a small variance normal distribution for the spike part of the mixture. 
#' It also implements additional functions to allow pre-processing procedures for the data 
#' and to fit the model for a large number of probesets or genes. Includes functions to:
#' 
#'   * pre-process a set of matrices;
#'   * fit models to a set of matrices;
#'   * detailed summary of model fit.
#' 
#' @docType package
#' @name slfm-package
NULL