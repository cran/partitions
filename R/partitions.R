"as.matrix.partition" <- function(x, ...){
  class(x) <- "matrix"
  NextMethod("as.matrix")
}

"summary.partition" <- function(object, ...){
  jj <- list(...)
  if(length(jj)>0){
    if(!is.null(jj$n)){
      n <- jj$n
    } else {
      n <- jj[[1]]
    }
  } else {
    n <- 10
  }
  if(ncol(object)>2*n){
    jj <- cbind(object[,seq_len(n)],object[,ncol(object)+1-seq_len(n)])
    shortened <- TRUE
  } else {
    jj <- object
    shortened <- FALSE
  }
  out <- list(shortened=shortened, n=n, out=jj)
  class(out) <- "summary.partition"
  return(out)
}

print.summary.partition <- function(x, ...){
  if(x$shortened){
    jj <- seq_len(x$n)
    x <- x$out
    x <- cbind(x[,jj],"...",x[,ncol(x)+1-jj])
  } else {
    x <- x$out
  }
  print.partition(x)
}
  
"setparts" <- function(x){
  if(length(x)==1){return(Recall(parts(x)))}
  if(is.matrix(x)){
    out <- do.call("cbind",apply(x,2,setparts))
  } else {
    x <- sort(x[x>0], decreasing=TRUE)
    num.of.parts <-
      factorial(sum(x))/(prod(c(factorial(x),factorial(table(x)))))
    out <- .C("wrap",
              as.integer(x),
              as.integer(length(x)),
              ans = integer(sum(x)*num.of.parts),
              PACKAGE="partitions"
              )$ans
    dim(out) <- c(sum(x),num.of.parts)
  }
  class(out) <- c("set_partition","partition")
  return(out)
}

"print.partition" <- function(x, mat=getOption("matrixlike"), h=getOption("horiz"), ...){
  class(x) <- "matrix"
  if(!isTRUE(mat)){
    colnames(x) <- rep(" ", ncol(x))
  }
  if(isTRUE(h)){
    x <- t(x)
  }
  return(invisible(print(noquote(x))))
}

"parts" <-
function(n){
  if(length(n)>1){
    stop("argument must be of length 1")
  }
  if(n>0){
    nn <- P(n)
    out <- .C("allparts",
              as.integer(n),
              as.integer(n*nn),
              ans = integer(n*nn),
              PACKAGE="partitions"
              )$ans
    dim(out) <- c(n,nn)
  } else {
    out <- matrix(0,0,0)
  }
  class(out) <- c("integer_partition","partition")
  return(out)
}

"diffparts" <-
function(n){
  if(length(n)>1){
    stop("argument must be of length 1")
  }
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
    class(out) <- c("diff_integer_partition","partition")
  return(out)
}

"restrictedparts" <- function(n, m, include.zero=TRUE, decreasing=TRUE){
  jj.n <- R(m,n,include.zero=include.zero)
  len <- m*jj.n

  out <- .C("allrestrictedparts",
           as.integer(m),
           as.integer(n),
           as.integer(len),
           as.integer(include.zero),
           ans=integer(len),
           PACKAGE="partitions"
           )$ans
  dim(out) <- c(m,jj.n)
  if(decreasing){
    out <- out[rev(seq_len(m)),]
  }
  class(out) <- c("restrictedpartition","partition")
  return(out)
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
  class(out) <- c("block_partition","partition")
  return(out)
}
    
"compositions" <-
function(n, m=NULL, include.zero=TRUE){
    if(!is.null(m)){
      if(include.zero){
        return(blockparts(rep(n,m),n))
      } else {
        return(1L + blockparts(rep(n-1,m),n))
      }
    }
    pad <- function(i,n){c(i,integer(n-length(i)))}
    f <- function(x){diff(c(as.integer(0),which(c(x,TRUE))))}
    jj <- apply(as.matrix(expand.grid(rep(list(c(TRUE,FALSE)),n-1))),1, f)
    out <-  do.call(cbind,lapply(jj, pad, n))
    rownames(out) <- NULL
    class(out) <- c("all_partitions","partition")
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
  if(nc>1){
    class(out) <- c("conjugate_partition","partition")
  }
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
