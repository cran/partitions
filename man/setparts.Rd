\name{setparts}
\alias{setparts}
\title{Set partitions}
\description{
  Enumeration of set partitions
}
\usage{
setparts(x)
}
\arguments{
  \item{x}{If a vector of length 1, the size of the set to be
    partitioned.  If a vector of length greater than 1, return all
    equivalence relations with equivalence classes with sizes of the
    elements of \code{x}.  If a matrix, return all equivalence classes
    with sizes of the columns of \code{x}}
}
\details{
  A \dfn{partition} of a set
  \eqn{S=\left\{1,\ldots,n\right\}}{S={1,...,n}} is a family of
  sets \eqn{T_1,\ldots,T_k}{T1,...,Tk} satisfying
  \itemize{
    \item \eqn{i\neq j\longrightarrow T_i\cap
      T_j=\emptyset}{union(Ti,Tj) empty if i != j}
    \item \eqn{\cup_{i=1}^kT_k=S}{union(T1,T2,...,Tk)=S}
    \item \eqn{T_i\neq\emptyset}{Ti not empty} for \eqn{i=1,\ldots,
      k}{1,...,k} 
  }
  
  The induced \dfn{equivalence relation} has \eqn{i\sim j}{i~j} if
  and only if \eqn{i} and \eqn{j} belong to the same partition.
  
  There are exactly fifteen ways to partition a set of four
  elements:
  \tabular{l}{
    \eqn{(1234)}\cr
    \eqn{(123)(4), (124)(3), (134)(2), (234)(1)}\cr
    \eqn{(12)(34), (13)(24), (14)(23)}\cr
    \eqn{(12)(3)(4), (13)(2)(4), (23)(1)(4), (24)(1)(3),
      (34)(1)(2)}\cr
    \eqn{(1)(2)(3)(4)}
  }
  
  Note that \eqn{(12)(3)(4)} is the same partition as, for example,
  \eqn{(3)(4)(21)} as the equivalence relation is the same.
  
  Consider partitions of a set \eqn{S} of five elements (named
  \eqn{1,2,3,4,5}) with sizes 2,2,1.  These may be enumerated as
  follows:
  
  \preformatted{
> u <- c(2,2,1)
> setparts(u)
                                  
[1,] 1 1 1 1 1 1 1 1 1 1 1 1 3 3 3
[2,] 2 2 3 1 1 1 2 2 3 2 2 3 1 1 1
[3,] 3 2 2 3 2 2 1 1 1 3 2 2 2 1 2
[4,] 2 3 2 2 3 2 3 2 2 1 1 1 2 2 1
[5,] 1 1 1 2 2 3 2 3 2 2 3 2 1 2 2

}

  See how each column has two 1s, two 2s and one 3.  This is because the
  first and second classes have size two, and the third has size one.

  The first partition, \code{x=c(1,2,3,2,1)}, is read \dQuote{class 1
  contains elements 1 and 5 (because the first and fifth element of
  \code{x} is 1); class 2 contains elements 2 and 4 (because the second
  and fourth element of \code{x} is 2); and class 3 contains element 3
  (because the third element of \code{x} is 3)}.  Formally, class
  \code{i} has elements \code{which(x==u[i])}.
}
\value{

  Returns a matrix each of whose columns show a set partition; an object
  of class \code{c("set_partition","partition")}. Type
  \code{?print.partition} to see how to change the options for printing.
}
\references{
  \itemize{
    \item R. K. S. Hankin 2006.  \emph{Additive integer partitions in
    \R}.  Journal of Statistical Software, Code Snippets 16(1)
    \item R. K. S. Hankin 2007.  \emph{Urn sampling without replacement:
      enumerative combinatorics in \R}.  Journal of Statistical
      Software, Code Snippets 17(1)
  }
}
\author{Luke G. West (\code{C++}) and Robin K. S. Hankin (\R)}
\seealso{\code{\link{parts}}, \code{\link{print.partition}}}
\examples{
setparts(4)                # all partitions of a set of 4 elements

setparts(c(3,3,2))         # all partitions of a set of 8 elements
                           # into sets of sizes 3,3,2.


jj <- restrictedparts(5,3)
setparts(jj)               # partitions of a set of 5 elements into
                           # at most 3 sets

setparts(conjugate(jj))    # partitions of a set of 5 elements into
                           # sets not exceeding 3 elements


setparts(diffparts(5))     # partitions of a set of 5 elements into
                           # sets of different sizes


}
\keyword{math}
