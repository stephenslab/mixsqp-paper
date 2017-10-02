\name{normmix.data}
\docType{data}
\alias{normmix.data}
\title{5000 x 20 conditional likelihood matrix computed from simulated
  data set.}
\description{
  This is an example data set that can be to compare algorithms for
  fitting mixture models. The 5000 x 20 conditional likelihood matrix
  was computed assuming the likelihood is univariate normal with zero
  mean and standard deviation \code{0.1}, and the prior is a mixture of
  20 univariate normals with zero mean and standard deviations
  \code{c(0.01,10^(seq(-2,0,length.out = 19)))}. The maximum-likelihood
  solution is also included, which can be used for testing the fitting
  algorithms.
}
\usage{data(normmix.data)}
\format{
  \code{normmix.data} is a list with three components:
  
  \describe{
    \item{L}{The 5000 x 20 conditional likelihood matrix.}

    \item{w}{The maximum-likelihood estimate of the mixture weights
      estimated using a primal-dual interior-point optimization method.}

    \item{obj}{The estimated objective minimum (\emph{i.e.}, the maximum
      log-likelihood).}    
  }
}
\details{
  To generate the conditional likelihood matrix, 5000 samples X were
  generated as follows: (1) probability \code{w[i]}, sample U from the
  univariate normal with zero mean and standard deviation \code{s[i]}q; (2)
  sample X from the univariate normal with mean U and standard
  deviation 0.1.

  The data simulation parameters were set to \code{s = c(0,0.1,0.2,0.5)}
  and \code{w = c(0.95,0.03,0.01,0.01))}. The pseudorandom number
  generator seed was set to 1; \code{set.seed(1)}.
}
\keyword{datasets}
