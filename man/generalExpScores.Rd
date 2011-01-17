\name{generalExpScores}
\alias{generalExpScores}
\title{
  Expected Scores for the General Model
}
\description{
  Computes the expected score as the aiming location varies across the
  board, using the specified covariance matrix and the general
  model for dart throws. 
}
\usage{
generalExpScores(Sig)
}
\arguments{
  \item{Sig}{
    The covariance matrix.
  }
}
\value{
  \item{e}{
    The matrix of expected scores.
  }
}
\author{
  Ryan Tibshirani
}
\examples{
# Generate the matrix of expected scores if my covariance is
# Sig=c(15^2,30^2,0), in the general model. This will take a
# couple of seconds. 
e = generalExpScores(Sig=(c(15^2,30^2,0)))

# Draw a heatmap!
drawHeatmap(e)
}