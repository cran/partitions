\name{parts}
\alias{parts}
\alias{diffparts}
\alias{restrictedparts}
\alias{blockparts}
\alias{allparts}
\title{Enumerate the partitions of an integer}
\description{
  Given an integer, return a matrix whose columns enumerate its
  partitions.

  Function \code{parts()} returns the unrestricted partions; function
  \code{diffparts()} returns the unequal partitions; function
  \code{restrictedparts()} returns the restricted partitions;
  function \code{blockparts()} returns the partitions subject to
  specified maxima; and function \code{allparts()} returns all
  partitions, of an integer.
 }
\usage{
parts(n)
diffparts(n)
restrictedparts(n, m, include.zero=TRUE, decreasing=TRUE)
blockparts(f, n=NULL, include.fewer=FALSE)
allparts(n)
}
\arguments{
  \item{n}{Integer to be partitioned.  In function \code{blockparts()},
    the default of \code{NULL} means to return all partitions of any size}
  \item{m}{In \code{restrictedparts()}, the order of the partition}
  \item{include.zero}{In \code{restrictedparts()}, Boolean with
    default \code{FALSE} meaning to include only partitions of \eqn{n}
    into \emph{exactly} \eqn{m} parts; and \code{TRUE} meaning to
    include partitions of \eqn{n} into \emph{at most} \eqn{m} parts
    (because zero parts are included)}
  \item{include.fewer}{In function \code{blockparts()}, Boolean with
    default \code{FALSE} meaning to return partitions into
    \emph{exactly} \code{n} and \code{TRUE} meaning to return partitions
    into \emph{at most} \code{n}}
  \item{decreasing}{In \code{restrictedparts()}, Boolean with default
    \code{TRUE} meaning to return partitions whose parts are in
    decreasing order and \code{FALSE} meaning to return partitions in
    lexicographical order, as appearing in Hindenburg's
    algorithm.   Note that setting to \code{decreasing} to \code{FALSE}
    has the effect of making \code{conjugate()} return garbage}
  \item{f}{In function \code{blockparts()}, a vector of strictly
    positive integers that gives the maximal number of blocks; see
    details} 
}
\details{
  Function \code{parts()} uses the algorithm in Andrews.  Function
  \code{diffparts()} uses a very similar algorithm that I have not seen
  elsewhere.

  Function \code{restrictedparts()} uses the algorithm in
  Andrews, originally due to Hindenburg.  For partitions into at most
  \eqn{m} parts, the same Hindenburg's algorithm is used but with a
  start vector of \code{c(rep(0,m-1),n)}.

  Function \code{blockparts()} enumerates the compositions of an integer
  subject to a maximum criterion: given vector
  \eqn{y=(y_1,\ldots,y_n)}{y=(y_1,...,y_p)} all sets of
  \eqn{a=(a_1,\ldots,a_n)}{a=(a_1,...,a_p)} satisfying
  \eqn{\sum_{i=1}^pa_i=n}{sum(a_i)=n} subect to \eqn{0<a_i\leq
    y_i}{0<a_i<y_i} for all \eqn{i} are given in lexicographical order.
  If argument \code{y} includes zero elements, these are treated
  consistently (ie a position with zero capacity).

  If \code{n} takes its default value of \code{NULL}, then
  \eqn{\sum_{i=1}^pa_i=n}{sum(a_i)=n} is removed (the numbers may sum to
  anything).  Note that these solutions are not necessarily in standard
  form, so functions \code{durfee()} and \code{conjugate()} may fail.

  Function \code{allparts()} returns all \eqn{2^n} ways of partitioning
  an integer; thus \code{4+1+1} is distinct from \code{1+4+1} or
  \code{1+1+4}.  This function is different from all the others in the
  package in that it is written in \R; it is not clear that C
  would be any faster
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

blockparts(1:4)
blockparts(1:4,3)
blockparts(1:4,3,include.fewer=TRUE)

allparts(3)

}
\keyword{math}
