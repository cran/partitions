\name{P}
\alias{P}
\alias{Q}
\alias{R}
\alias{S}
\title{Number of partitions of an integer}
\description{
  Given an integer, \code{P()} returns the number of additive
  partitions, \code{Q()} returns the  number of unequal
  partitions, and \code{R()} returns the number of
  restricted partitions.  Function \code{S()} returns the number of
  block partitions.
}
\usage{
P(n, give = FALSE)
Q(n, give = FALSE)
R(m, n, include.zero = FALSE)
S(n = NULL, y, include.fewer = FALSE)
}
\arguments{
  \item{n}{Integer whose partition number is desired.  In function
    \code{S()}, the default of \code{NULL} means to return the number of
  partitions of any size}
  \item{m}{In function \code{R()}, the order of the
    decomposition}
  \item{give}{Boolean, with default \code{FALSE} meaning to return just
    \code{P(n)} or \code{Q(n)} and \code{TRUE} meaning to return
    \code{P(1:n)} or \code{Q(1:n)} (this option takes no extra
    computation)}
  \item{include.zero}{In \code{restrictedparts()}, Boolean with
    default \code{FALSE} meaning to count only partitions of \eqn{n}
    into \emph{exactly} \eqn{m} parts; and \code{TRUE} meaning to
    include partitions of \eqn{n} into \emph{at most} \eqn{m} parts
    (because parts of zero are included)}
  \item{include.fewer}{In function \code{blockparts()}, Boolean
    with default \code{FALSE} meaning to return partitions into
    \emph{exactly} \code{n} and \code{TRUE} meaning to return partitions
    into \emph{at most} \code{n}}
  \item{y}{In function \code{S()}, the stack vector}
}
\details{
  Functions \code{P()} and \code{Q()} use Euler's
  recursion formula.  Function \code{R()} enumerates the partitions 
  using Hindenburg's method (see Andrews) and counts them until the
  recursion bottoms out.
}
\references{G. E. Andrews. \dQuote{The Theory of Partitions},
  Cambridge University Press, 1998.
  }
\author{Robin K. S. Hankin}
\note{
  Functions \code{P()} and \code{Q()} use \code{unsigned long long}
  integers, a type which is system-dependent.  For me, \code{P()} works
  for \eqn{n} equal to or less than 416, and \code{Q()} works for
  \eqn{n} less than or equal to 792.  YMMV;  none of the methods test
  for overflow, so use with care!
}
\examples{
P(10,give=TRUE)
Q(10,give=TRUE)
R(10,20,include.zero=FALSE)
R(10,20,include.zero=TRUE)

S(5,1:4)

}
\keyword{math}
