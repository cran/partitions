\name{parts}
\alias{parts}
\alias{diffparts}
\alias{restrictedparts}
\title{Enumerate the partitions of an integer}
\description{
  Given an integer, return a matrix whose columns enumerate its
  partitions.

  Function \code{parts()} returns the unrestricted partions; function
  \code{diffparts()} returns the unequal partitions; function
  \code{restrictedparts()} returns the restricted partitions.
 }
\usage{
parts(n)
diffparts(n)
restrictedparts(n, m, include.zero=TRUE, decreasing=TRUE)
}
\arguments{
  \item{n}{Integer to be partitioned}
  \item{m}{In \code{restrictedparts()}, the order of the partition}
  \item{include.zero}{In \code{restrictedparts()}, Boolean with
    default \code{FALSE} meaning to include only partitions of \eqn{n}
    into \emph{exactly} \eqn{m} parts; and \code{TRUE} meaning to
    include partitions of \eqn{n} into \emph{at most} \eqn{m} parts
    (because zero parts are included)}
  \item{decreasing}{In \code{restrictedparts()}, Boolean with default
    \code{TRUE} meaning to return partitions whose parts are in
    decreasing order and \code{FALSE} meaning to return partitions in
    lexicographical order, as appearing in Hindenburg's
    algorithm.   Note that setting to \code{decreasing} to \code{FALSE}
    has the effect of making \code{conjugate()} return garbage.}    
}
\details{
  Function \code{parts()} uses the algorithm in Andrews.  Function
  \code{diffparts()} uses a very similar algorithm that I have not seen
  elsewhere.

  Function \code{restrictedparts()} uses the algorithm in
  Andrews, originally due to Hindenburg.  For partitions into at most
  \eqn{m} parts, the same Hindenburg's algorithm is used but with a
  start vector of \code{c(rep(0,m-1),n)}.
}
\author{Robin K. S. Hankin}
\references{
G. E. Andrews. \dQuote{The Theory of Partitions},
Cambridge University Press, 1998
}
\examples{
parts(5)
diffparts(10)
restrictedparts(9,4)
restrictedparts(9,4,FALSE)
restrictedparts(9,4,decreasing=TRUE)
}
\keyword{math}
