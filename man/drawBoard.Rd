\name{drawBoard}
\alias{drawBoard}
\title{
  Dartboard
}
\description{
  Draws the dartboard, either as a new plot or on top of an existing
  one. 
}
\usage{
drawBoard(new = FALSE, lines = TRUE, numbers = TRUE, outside = TRUE,
          col = "black", ...)
}
\arguments{
  \item{new}{
    Make a new plot?
  }
  \item{lines}{
    Draw the lines separating the regions?
  }
  \item{numbers}{
    Draw the numbers?
  }
  \item{outside}{
    Draw the numbers outside (or inside) the dartboard?
  }
  \item{col}{
    The color for the lines.
  }
  \item{\dots}{
    More arguments for drawing the numbers (text). 
  }
}
\author{
  Ryan Tibshirani
}
\examples{
# Draw a new dartboard
drawBoard(new=TRUE)
}
