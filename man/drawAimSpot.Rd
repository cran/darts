\name{drawAimSpot}
\alias{drawAimSpot}
\title{
  Optimal Aiming Spot
}
\description{
  Draws the optimal aiming location, i.e. the spot with the highest
  expected score, on top of an existing plot. 
}
\usage{
drawAimSpot(e, col = "blue", pch = 19, ...)
}
\arguments{
  \item{e}{
    A matrix of the expected scores, created by the function
    \code{simpleExpScores} or \code{generalExpScores}. 
  }
  \item{col}{
    The color of the dot.
  }
  \item{pch}{
    The plotting "character" for the dot.
  }
  \item{\dots}{
    More arguments for drawing the dot.
  }
}
\author{
  Ryan Tibshirani
}
\examples{
# This will take a couple of seconds
e = simpleExpScores(s=0)

# This is what your heatmap would look like if you had perfect accuracy  
par(mfrow=c(1,2))
drawHeatmap(e)
drawBoard(new=TRUE)
drawAimSpot(e)
}