\name{drawHeatmap}
\alias{drawHeatmap}
\title{
  Heatmap of Expected Scores
}
\description{
  Draws a heatmap of the expected score as the aiming location varies
  across the dartboard, as a new plot.
}
\usage{
drawHeatmap(e, col = heat.colors(30))
}
\arguments{
  \item{e}{
    A matrix of the expected scores, created by the function
    \code{simpleExpScores} or \code{generalExpScores}. 
  }
  \item{col}{
    The colors to use for the heatmap.
  }
}
\author{
  Ryan Tibshirani
}
\examples{
# This will take a couple of seconds
e = simpleExpScores(s=0)

# This is what your heatmap would look like if you had
# perfect accuracy
drawHeatmap(e)
}