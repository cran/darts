\name{simpleExpScores}
\alias{simpleExpScores}
\title{
  Expected Scores for the Simple Model
}
\description{
  Computes the expected score as the aiming location varies across the
  board, using the specified variance and the simple model for dart
  throws. 
}
\usage{
simpleExpScores(s)
}
\arguments{
  \item{s}{
    The marginal variance.
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
# Generate the matrix of expected scores if my variance is 25^2,
# in the simple model. This will take a couple of seconds.
e = simpleExpScores(s=25^2)

# Draw a heatmap!
drawHeatmap(e)
}
