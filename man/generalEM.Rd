\name{generalEM}
\alias{generalEM}
\title{
  EM Algorithm for the General Model
}
\description{
  EM algorithm to estimate your variance based on your scores, in the
  general model. 
}
\usage{
generalEM(x, Sig.init = c(10^2, 10^2, 0.1 * 10 * 10), niter = 100,
          seed = NULL)
}
\arguments{
  \item{x}{
    Scores of throws aimed at the center of the dartboard.
  }
  \item{Sig.init}{
    The initial guess for the covariance matrix, represented as a vector:
    x marginal variance, then y marginal variance, then x-y covariance.  
  }
  \item{niter}{
    The number of iterations.
  }
  \item{seed}{
    The seed for the random number generator (the E-step is done by
    importance sampling). 
  }
}
\value{
  \item{Sig.final}{
    The final estimate of the covariance matrix.
  }
  \item{Sig.init}{
    The initial estimate of the covariance matrix.
  }
  \item{Sig}{
    The estimate of the covariance at each iteration. 
  }
  \item{loglik}{
    The log likelihood at each iteration---currently not implemented
    (this is just an array of 0s).
  }
  \item{niter}{
    The number of iterations.
  }
}
\author{
  Ryan Tibshirani
}
\examples{
# Scores of 100 of my dart throws, aimed at the center of the board
x = c(12,16,19,3,17,1,25,19,17,50,18,1,3,17,2,2,13,18,16,2,25,5,5,
1,5,4,17,25,25,50,3,7,17,17,3,3,3,7,11,10,25,1,19,15,4,1,5,12,17,16,
50,20,20,20,25,50,2,17,3,20,20,20,5,1,18,15,2,3,25,12,9,3,3,19,16,20,
5,5,1,4,15,16,5,20,16,2,25,6,12,25,11,25,7,2,5,19,17,17,2,12)

# Get my variance in the general Gaussian model
a = generalEM(x,niter=100,seed=0)

# The EM estimate of my covariance matrix
Sig = a$Sig.final
}