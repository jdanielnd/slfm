\name{slfm_list}
\alias{slfm_list}
\title{Fit SLFM to the matrices inside a directory}
\usage{
slfm_list(path = ".", recursive = TRUE, ite, a = 2.1, b = 1.1,
  gamma_a = 1, gamma_b = 1, omega = 10, omega_1 = 0.01, burnin = 500)
}
\arguments{
  \item{path}{path to read the matrices from}

  \item{recursive}{if the function should look recursively
  inside folders}

  \item{ite}{number of iterations of the MCMC algorithm}

  \item{a}{prior shape parameter for Gamma distribution}

  \item{b}{prior scale parameter for Gamma distribution}

  \item{gamma_a}{prior parameter for Beta distribution}

  \item{gamma_b}{prior parameter for Beta distribution}

  \item{omega}{prior variance of the slab component}

  \item{omega_1}{prior variance of the spike component}

  \item{burnin}{burn-in size}
}
\description{
This function is used to fit a Bayesian sparse latent
factor model to a directory of numeric matrices.
}

