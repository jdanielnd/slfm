\name{slfm}
\alias{slfm}
\title{Sparse Latent Factor Model}
\usage{
slfm(x, ite, a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1, omega = 10,
  omega_1 = 0.01, burnin = round(0.25 * ite))
}
\arguments{
  \item{x}{matrix with the pre-processed data}

  \item{ite}{number of iterations of the MCMC algorithm}

  \item{a}{prior shape parameter for Gamma distribution}

  \item{b}{prior scale parameter for Gamma distribution}

  \item{gamma_a}{prior parameter for Beta distribution}

  \item{gamma_b}{prior parameter for Beta distribution}

  \item{omega}{prior variance of the slab component}

  \item{omega_1}{prior variance of the spike component}

  \item{burnin}{burn-in size}
}
\value{
x: data matrix

p_star: matrix of MCMC chains for p_star parameter

alpha: summary table of MCMC chains for alpha parameter

lambda: summary table of MCMC chains for lambda parameter

sigma: summary table of MCMC chains for sigma parameter

classification: classification of each alpha (`present`,
`marginal`, `absent`)
}
\description{
This function is used to fit a Bayesian sparse latent
factor model.
}
\examples{
mat <- matrix(rnorm(2000), nrow = 20)
slfm(mat, ite = 1000)
}

