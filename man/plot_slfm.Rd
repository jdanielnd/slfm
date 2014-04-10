\name{plot_slfm}
\alias{plot_slfm}
\title{Plot the data matrix along with probabilities box-plots}
\usage{
plot_slfm(y, result, clas, stand = TRUE, reordRow = TRUE, reordCol = TRUE,
  highContrast = TRUE)
}
\arguments{
  \item{y}{data matrix}

  \item{result}{matrix of p_value chains from SLFM}

  \item{stand}{indicates whether the data should be
  standardized}

  \item{reordRow}{indicates whether the rows should be
  ordered}

  \item{reordCol}{indicates whether the cols should be
  ordered}

  \item{highContrast}{indicates if high contrast should be
  used}
}
\description{
This function is used to fit a Bayesian sparse latent
factor model to a directory of numeric matrices.
}

