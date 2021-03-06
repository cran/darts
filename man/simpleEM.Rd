\name{simpleEM}
\alias{simpleEM}
\title{
  EM Algorithm for the Simple Model
}
\description{
  EM algorithm to estimate your variance based on your scores, in the 
  simple model. 
}
\usage{
simpleEM(x, s.init = 100, niter = 100)
}
\arguments{
  \item{x}{
    Scores of throws aimed at the center of the dartboard.
  }
  \item{s.init}{
    The initial guess for the marginal variance. 
  }
  \item{niter}{
    The number of iterations. 
  }
}
\value{
  \item{s.final}{
    The final estimate of the variance.
  }
  \item{s.init}{
    The initial estimate of the variance.
  }
  \item{s}{
    The estimate of the variance at each iteration of the
    EM algorithm.
  }
  \item{loglik}{
    The (observed) log likelihood at each iteration of the EM algorithm.   
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

# Get my variance in the simple Gaussian model
a = simpleEM(x,niter=100)

# Check the log likelihood
plot(1:a$niter,a$loglik,type="l",xlab="Iteration",ylab="Log likelihood")

# The EM estimate of my variance
s = a$s.final
}