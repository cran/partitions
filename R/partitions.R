"parts" <-
function(n){
  nn <- P(n)
  out <- .C("allparts",
           as.integer(n),
           as.integer(n*nn),
           ans = integer(n*nn),
           PACKAGE="partitions"
           )$ans
  dim(out) <- c(n,nn)
  return(out)
}

"diffparts" <-
function(n){
  nn <- Q(n)
  n.tri <- floor( (sqrt(1+8*n)-1)/2)
  out <- .C("alldiffparts",
            as.integer(n),
            as.integer(n.tri*nn),
            as.integer(n.tri),
            ans = integer(n.tri*nn+1),
            PACKAGE="partitions"
            )$ans[-nn*n.tri-1]
    dim(out) <- c(n.tri,nn)
  return(out)
}

"restrictedparts" <- function(n, m, include.zero=TRUE, decreasing=TRUE){
  jj.n <- R(m,n,include.zero=include.zero)
  len <- m*jj.n

  jj <- .C("allrestrictedparts",
           as.integer(m),
           as.integer(n),
           as.integer(len),
           as.integer(include.zero),
           ans=integer(len),
           PACKAGE="partitions"
           )
  out <- jj$ans
  dim(out) <- c(m,jj.n)
  if(decreasing){
    return(out[m:1,])
  } else {
    return(out)
  }
}

"P" <-
function(n, give=FALSE){
  n <- n+1
  jj <- .C("numbparts",
           as.integer(n),
           ans = double(n),
           PACKAGE = "partitions"
           )
  if(give){
    return(jj$ans[-1])
  } else {
    return(jj$ans[n])
  }
}

"R" <- function(m,n, include.zero=FALSE){
  stopifnot(m <= n)
  if(include.zero){
    start <- c(rep(0,m-1),n)
  } else {
    start <- c(rep(1,m-1),n-m+1)
  }
  jj <- .C("numbrestrictedparts_R",
           as.integer(start),
           as.integer(m),
           ans = as.integer(m),
           PACKAGE="partitions"
           )
  return(jj$ans)
}

"Q" <- function(n, give=FALSE){
  n <- n+1
  jj <- .C("numbdiffparts",
           as.integer(n),
           ans = double(n),
           PACKAGE = "partitions"
           )
  if(give){
    return(jj$ans)
  } else {
    return(jj$ans[n])
  }
}


"conjugate" <- function(x){
  x <- as.matrix(x)
  mx <- max(x)
  nc <- ncol(x)
  out <- .C("conjugate",
           as.integer(x),
           as.integer(nrow(x)),
           as.integer(nc),
           as.integer(mx),
           ans=integer(mx*nc),
           PACKAGE = "partitions"
           )$ans
  dim(out) <- c(mx,nc)
  return(drop(out))
}

"durfee" <- function(x){
  x <- as.matrix(x)
  .C("durfee",
     as.integer(x),
     as.integer(nrow(x)),
     as.integer(ncol(x)),
     ans=integer(ncol(x)),
     PACKAGE="partitions")$ans
}
