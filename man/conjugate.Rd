\name{conjugate}
\alias{conjugate}
\alias{durfee}
\alias{Durfee}
\title{Conjugate partitions and Durfee squares}
\description{
Given a partition, provide its conjugate
}
\usage{
conjugate(x)
durfee(x)
}
\arguments{
  \item{x}{Either a vector describing a partition, in standard form
    (ie nonincreasing); or a matrix whose columns are partitions in
    standard form}
}
\details{
  Conjugation is described in Andrews, and (eg) Hardy and Wright.

  Essentially, \code{conjugate()} carries out
  R idiom
  \code{rev(cumsum(table(factor(a[a>0],levels=max(a):1))))}, but faster.

  The \dQuote{Durfee square} of a partition is defined on page 281 of
  Hardy and Wright.   It is the largest square of nodes contained in the
  partition's Ferrers graph.  Function \code{durfee()} returns the side
  of the Durfee square which Andrews denotes
  \eqn{d(\lambda)}{d(lambda)}.  It is equivalent to R idiom
  \code{function(a){sum(a>=1:length(a))}}, but faster.
  
}
\value{
Returns either a partition in standard form, or a matrix whose
columns are partitions in standard form.
}
\author{Robin K. S. Hankin}
\note{If argument \code{x} is not nonincreasing, all bets are off.  These
  functions will not work; they will silently return garbage.  Caveat emptor!
}
\examples{
parts(5)
conjugate(parts(5))

restrictedparts(6,4)
conjugate(restrictedparts(6,4))

durfee(10:1)

}
\keyword{math}
