\name{print.partition}
\alias{print.partition}
\alias{print.summary.partition}
\title{Print methods for partition object}
\description{
  A print method for partition objects and summary partition objects,
  including various configurable options
}
\usage{
print.partition(x, mat = getOption("matrixlike"), h =
getOption("horiz"), ...)
print.summary.partition(x, ...)
}
\arguments{
  \item{x}{Object to be printed: an object of class either
    \code{partition} or \code{summary.partition}}
  \item{mat}{Boolean, with \code{TRUE} meaning to print like a matrix,
    and any other value meaning to print without column names (which
    usually results in more compact appearance)}
  \item{h}{Boolean governing the orientation of the printed matrix, with
    \code{TRUE} meaning to print with the rows being the partitions and 
    any other value (the default) meaning to print the transpose}
  \item{\dots}{Further arguments provided for compatibility}
}
\author{Robin K. S. Hankin}
\examples{
print(parts(5))

summary(parts(7))
}
\keyword{math}
