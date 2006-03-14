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

"blockparts" <- function(f,n=NULL,include.fewer=FALSE){
  s <- sum(f)
  if(is.null(n)){
    return(Recall(c(s,f),s)[-1,])
  }
  if(s<n){stop("too many blocks: n<sum(f)")}
  if(include.fewer){
    return(Recall(c(s,f),n)[-1,])
  }
  nf <- names(f)
  f <- as.vector(f)
  fnz <- f[f>0]
  nb <- S(fnz,n)
  lfnz <- length(fnz)
  lf <- length(f)
  out <- .C("allblockparts",
           ans=integer(lfnz*nb),
           as.integer(fnz),
           as.integer(nb),
           as.integer(lfnz),
           as.integer(n)
           )$ans
  dim(out) <- c(lfnz,nb)
  if(any(f<1)){
    out <- replace(matrix(0,lf,ncol(out)),f>0,out)
  }
  rownames(out) <- nf
  colnames(out) <- rep(" ",nb)
  return(out)
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

"S" <- function(f,n=NULL,include.fewer=FALSE){
  if(length(n)>1){
    stop("In function S(), n [the second argument] must be an integer (or NULL).  Check for the first and second arguments being transposed")
  }
    p <- polynomial(1)
    for(i in f){
      p <- p * polynomial(rep.int(1, i + 1))
      p <- polynomial(p[1:min(length(p),n+1)])
    }
    if(include.fewer){
      return(sum(p[1:(n+1)]))
    } else{
      return(p[n+1])
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

if(FALSE){
"U" <- function(y,naive=FALSE){
  if(naive){
    return(
           factorial(sum(y))/prod(factorial(y))
           )
  } else {
    stop("not implemented")
  }
}

"perms" <- function(y){
  n <- length(y)
  x <- rep(1:n,y)
  nn <- U(y)
  
  out <- .C("allperms",
           as.integer(n),
           as.integer(n*nn),
           ans = integer(n*nn),
           PACKAGE="partitions"
           )$ans
  dim(out) <- c(n,nn)
  return(out)
}
}
