\name{partitions-package}
\alias{partitions-package}
\alias{partitions}
\docType{package}
\title{
Integer partitions
}
\description{
Routines to enumerate all partitions of an integer; includes restricted
and unequal partitions.
}
\details{
This tiny package comprises six functions: \code{P()}, \code{Q()}, and
\code{R()} give the number of partitions, unequal partitions, and
restricted partitions of an integer.

Functions \code{parts()}, \code{diffparts()}, and
\code{restrictedparts()} enumerate these partitions.

Function \code{conjugate()} gives the conjugate of a partition and
function \code{durfee()} gives the size of the Durfee square.

\strong{NB} the emphasis in this package is terse, efficient C code.
This means that there is a minimum of argument checking.  For example,
function \code{conjugate()} assumes that the partition is in standard
form (ie nonincreasing); supplying a vector in nonstandard form will
result in garbage being returned silently.

}
\author{Robin K. S. Hankin}
\references{
G. E. Andrews 1998 \dQuote{The Theory of Partitions},
Cambridge University Press

M. Abramowitz and I. A. Stegun 1965.  \dQuote{Handbook of Mathematical
Functions}, New York: Dover

G. H. Hardy and E. M. Wright 1985 \dQuote{An introduction to the theory of
  numbers}, Clarendon Press: Oxford (fifth edition)

R. K. S. Hankin 2006. \dQuote{Additive Integer Partitions in R}.
Journal of Statistical Software, volume 16, codesnippet 1.

}
\keyword{package}
\examples{
 parts(5)
 diffparts(9)
 restrictedparts(15,10)
 P(10,give=TRUE)
 Q(10,give=TRUE)
 R(5,10)
}
