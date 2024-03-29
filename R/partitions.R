".fac" <- function(x){  # exact factorial; NB returns (by design) a
                        # bigz; used in setparts()
  out <- as.bigz(1)
  for(n in x){
    out <- out * prod(as.bigz(seq_len(n)))
  }
  return(out)
}

"as.matrix.partition" <- function(x, ...){
  x <- unclass(x)
  NextMethod("as.matrix")
}

"as.partition" <- function(x, ...){
  storage.mode(x) <- "integer"
  class(x) <- "partition"
  return(x)
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
  if(length(x)==1){
    if (x < 1){
      stop("if single value, x must be >= 1")
    }
    else if (x == 1){
      out <- matrix(1, 1, 1)
    } else
    return(Recall(parts(x)))
  }
  if(is.matrix(x)){
    out <- do.call(
      "cbind",
      lapply(seq_len(ncol(x)), function(i) setparts(x[, i]))
    )
  } else {
    x <- sort(x[x>0], decreasing=TRUE)
    num.of.parts <-
      as.integer(.fac(sum(x))/(prod(c(.fac(x),.fac(table(x))))))
    out <- .C("c_wrap",
              as.integer(x),
              as.integer(length(x)),
              ans = integer(sum(x)*num.of.parts),
              PACKAGE="partitions"
              )$ans
    dim(out) <- c(sum(x),num.of.parts)
  }
  as.partition(out)
}

"listParts" <- function(x,do.set=FALSE) {
    jj <- setparts(x)
    if(do.set){
        out <- do.call(sets::set,apply(jj,2,vec_to_set))
    } else {
        out <- apply(jj, 2, vec_to_eq)
    }
    return(out)
}

"print.equivalence" <- function(x,sep=getOption("separator"), ...){
  if(is.null(sep)){sep <- ","}
  f <- function(x){paste(c("(",paste(x,collapse=sep),")"),collapse="")}
  out <- paste(unlist(lapply(x,f)),collapse="")
  x <-  unclass(x)
  return(invisible(print(noquote(out))))
}

"print.partition" <- function(x, mat=getOption("matrixlike"), h=getOption("horiz"), ...){
  x <- as.matrix(unclass(x))
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
      out <- .C("c_allparts",
                as.integer(n),
                as.integer(n*nn),
                ans = integer(n*nn),
                PACKAGE="partitions"
                )$ans
      dim(out) <- c(n,nn)
    } else {
      out <- matrix(0,0,0)
    }
    return(as.partition(out))
  }

"nextpart" <- function(part,check=TRUE){
  if(check){
    stopifnot(part==round(part))
    stopifnot(all(part >= 0))
    stopifnot(sum(part)==length(part))
    jj <- diff(rev(part))
    stopifnot(all(jj >= 0))
    stopifnot(any(jj >  0))
  }
  .C("c_nextpart",
     ans=as.integer(part),
     PACKAGE="partitions")$ans
}

"islastpart" <- function(part){all(part==1)}

"firstpart" <- function(n){
  as.integer(c(n,integer(n-1)))
}

".tri" <- function(n){floor( (sqrt(1+8*n)-1)/2)}

"diffparts" <-
function(n){
  if(length(n)>1){
    stop("argument must be of length 1")
  }
  nn <- Q(n)
  n.tri <- .tri(n)
  out <- .C("c_alldiffparts",
            as.integer(n),
            as.integer(n.tri*nn),
            as.integer(n.tri),
            ans = integer(n.tri*nn+1),
            PACKAGE="partitions"
            )$ans[-nn*n.tri-1]
  dim(out) <- c(n.tri,nn)
  return(as.partition(out))
}

"nextdiffpart" <- function(part,check=TRUE){
  n.tri <- .tri(sum(part))
  if(check){
    stopifnot(all(part==round(part)))
    stopifnot(all(part >= 0))
    stopifnot(length(part) == n.tri)
    jj <- diff(rev(part))
    stopifnot(all(jj >= 0))
    stopifnot( !(all( (jj==1) | (jj==2)) & sum(jj==2)==1))
    jj2 <- diff(rev(part[part>0]))
    stopifnot(all(jj2>0))
  }

  .C("c_nextdiffpart",
     ans   = as.integer(part),
     n.tri = as.integer(n.tri),
     PACKAGE="partitions")$ans
}

"islastdiffpart" <- function(part){
  jj <- diff(rev(part))
  all( (jj==1) | (jj==2)) & (sum(jj==2)==1)
}

"firstdiffpart" <- function(n){
  as.integer(c(n,integer(.tri(n)-1)))
}

"restrictedparts" <- function(n, m, include.zero=TRUE, decreasing=TRUE){
  if(m>n){  #NB: strict
    if(!include.zero){
      stop("m>n and include.zero=FALSE are incompatible")
    }
    jj <- Recall(n,n,include.zero=include.zero,decreasing=decreasing)
    if(decreasing){
      jj <- rbind(jj, matrix(0L,m-n,ncol(jj)))
    } else {
      jj <- rbind(matrix(0L,m-n,ncol(jj)),jj)
    }
    return(as.partition(jj))
  }

  jj.n <- R(m,n,include.zero=include.zero)
  len <- m*jj.n

  out <- .C("c_allrestrictedparts",
           as.integer(m),
           as.integer(n),
           as.integer(len),
           as.integer(include.zero),
           ans=integer(len),
           PACKAGE="partitions"
           )$ans
  dim(out) <- c(m,jj.n)
  if(decreasing){
    out <- out[rev(seq_len(m)),,drop=FALSE]
  }
  return(as.partition(out))
}

"nextrestrictedpart" <- function(part,check=TRUE){
  if(check){
    stopifnot(all(part==round(part)))
    stopifnot(all(part>=0))
    stopifnot(max(part)-min(part)>=2)
    stopifnot(all(diff(part)<=0))
  }
    rev(.C("c_nextrestrictedpart",
       ans=as.integer(rev(part)),
       ignore=as.integer(length(part)),
       PACKAGE = "partitions")$ans)
}

"islastrestrictedpart" <- function(part){
  max(part)-min(part) <= 1
}

"firstrestrictedpart" <- function(n, m, include.zero=TRUE){
  if(include.zero){
    return(as.integer(c(n,integer(m-1))))
  } else {
    return(as.integer(c(n-m+1,rep(1,m-1))))
  }
}

"blockparts" <- function(f,n=NULL,include.fewer=FALSE){
  s <- sum(f)
  if(is.null(n)){
    return(as.partition(Recall(c(s,f),s)[-1,]))
  }
  if(s<n){stop("too many blocks: n<sum(f)")}
  if(include.fewer){
    return(as.partition(Recall(c(s,f),n)[-1,]))
  }
  nf <- names(f)
  f <- as.vector(f)
  fnz <- f[f>0]
  nb <- S(fnz,n)
  lfnz <- length(fnz)
  lf <- length(f)
  out <- .C("c_allblockparts",
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
  return(as.partition(out))
}

"islastblockpart" <- function(part, f, n=NULL, include.fewer=FALSE){
  if(all(part==0)){
    if(is.null(n)){
      return(FALSE)
    } else {
      if(n==0){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }

  if(is.null(n)){
    return(all(part==f))
  }

  jj <- sum(cumprod(part==0)) + sum(cumprod(rev(part==f)))
  ans <- (jj==length(part)) | (jj==length(part)-1)

  if(include.fewer){
    return(ans & sum(part)==n)
  } else {
    return(ans)
  }

}

"nextblockpart" <- function(part, f, n=sum(part), include.fewer=FALSE, check=TRUE){
  if(check){
    stopifnot(all(part==round(part)))
    stopifnot(all(f==round(f)))
    stopifnot(all(part>=0))
    stopifnot(all(f>0))
    stopifnot(length(part) == length(f))
    stopifnot(all(part <= f))
  }
  s <- sum(f)

  if(is.null(n)){
    return(Recall(
                  part = c(s-sum(part),part),
                  f    = c(s,f),
                  n    = -1
                  )[-1]
           )
  }

  if(include.fewer){
    return(Recall(
                  part = c(n-sum(part),part),
                  f    = c(s,f),
                  n    = n
                  )[-1]
           )
  }

  .C("c_nextblockpart",
     ans=as.integer(part),
     as.integer(f),
     as.integer(length(part)),
     PACKAGE = "partitions")$ans
}

"firstblockpart" <- function(f, n=NULL, include.fewer=FALSE){
  if(is.null(n)){
    return(as.integer(f*0))
  }
  if(include.fewer){
    return(as.integer(f*0))
  }
  stopifnot(sum(f) >= n)
  out <- f*0
  jj <- cumsum(f)<n
  if(any(jj)){
    out[jj] <- f[jj]
    out[max(which(jj))+1] <- n-sum(out)
    return(as.integer(out))
  } else {
    return(as.integer(c(n,rep(0,length(f)-1))))
  }
}

"compositions" <-
function(n, m=NULL, include.zero=TRUE){
    if(!is.null(m)){
      if(include.zero){
        return(blockparts(rep(n,m),n))
      } else {
        return(1L + blockparts(rep(n-1,m),n-m))
      }
    }
    pad <- function(i,n){c(i,integer(n-length(i)))}
    f <- function(x){diff(c(as.integer(0),which(c(x,TRUE))))}
    jj <- apply(as.matrix(expand.grid(rep(list(c(FALSE,TRUE)),n-1))),1, f)
    out <-  do.call(cbind,lapply(jj, pad, n))
    rownames(out) <- NULL
    return(as.partition(out))
}

"islastcomposition" <- function(comp, restricted, include.zero=TRUE){
  if(restricted){
    m <- length(comp)
    if(include.zero){
      return(all(comp == c(integer(m-1),sum(comp))))
    } else {
      return(all(comp == c(rep(1,m-1),sum(comp)-m+1)))
    }
  } else {
    return(all(comp==1))
  }
}

"nextcomposition" <- function(comp, restricted, include.zero=TRUE, check=TRUE){
  if(check){
    stopifnot(!islastcomposition(comp, restricted=restricted, include.zero=include.zero))
    stopifnot(all(comp == round(comp)))
    if(restricted){
      if(include.zero){
        stopifnot(all(comp >= 0))
      } else {
        stopifnot(all(comp > 0))
      }
    } else {
      stopifnot(all(comp >= 0))
    }
  }
  n <- sum(comp)
  if(restricted){
    m <- length(comp)
    if(include.zero){
      return(nextblockpart(comp,rep(n,m),n,check=check))
    } else {
      return(1L+nextblockpart(comp-1L,rep(n-1,m),n-m,check=check))
    }
  } else {
    comp <- c(comp , integer(sum(comp)-length(comp)))
    return(rev(bintocomp(tobin(todec(rev(comptobin(comp)))+1,sum(comp)-1))))
  }
}

"firstcomposition" <-
function(n, m=NULL, include.zero=TRUE){
    if(!is.null(m)){
      if(include.zero){
        return(as.integer(c(n,rep(0,m-1))))
      } else {
        return(as.integer(c(n-m+1,rep(1,m-1))))
      }
    }
    return(as.integer(c(n,rep(0,n-1))))
}

"P" <-
function(n, give=FALSE){
  stopifnot(length(n)==1)
  n <- n+1
  jj <- .C("c_numbparts",
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
  stopifnot(length(m)==1)
  stopifnot(length(n)==1)
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
  stopifnot(length(n)==1)
  n <- n+1
  jj <- .C("c_numbdiffparts",
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

"conjugate" <- function(x, sorted = TRUE){
  if (!length(x))
    return(integer(0))
  x <- as.matrix(x)
  # if (!sorted)
  #   x <- apply(x, 2, sort, decreasing = TRUE)
  mx <- max(x)
  nc <- ncol(x)
  out <- .C("c_conjugate",
           as.integer(x),
           as.integer(nrow(x)),
           as.integer(nc),
           as.integer(mx),
           as.integer(sorted),
           ans=integer(mx*nc),
           PACKAGE = "partitions"
           )$ans
  dim(out) <- c(mx,nc)
  if(nc>1){
    out <- as.partition(out)
  }
  return(drop(out))
}

"durfee_sorted" <- function(x){
  x <- as.matrix(x)
    .C("c_durfee",
       as.integer(x),
       as.integer(nrow(x)),
       as.integer(ncol(x)),
       ans=integer(ncol(x)),
       PACKAGE="partitions")$ans
}

"durfee" <- function(x, sorted = TRUE){
  x <- as.matrix(x)
  if (sorted){
      return(durfee_sorted(x))
  } else {
      return(durfee_sorted(apply(x,2,sort,decreasing=TRUE)))
  }
}

"perms" <- function(n){
  stopifnot(length(n) ==1)
  stopifnot(n == round(n))
  nc <- factorial(n)  # nc = number of columns
  out <- .C("c_allperms",
            as.integer(seq_len(n)),
            as.integer(n),
            as.integer(nc),
            ans=integer(n*nc),
            PACKAGE="partitions"
            )$ans
  dim(out) <- c(n,nc)
  return(as.partition(out))
}

"tobin" <- function(n,len,check=TRUE){
  if(check){
    stopifnot(n == round(n))
    stopifnot(n >= 0)
    stopifnot(len == round(len))
    stopifnot(len >= 0)
  }

  .C("c_tobin",
     as.integer(n),
     ans=integer(len),
     as.integer(len),
     PACKAGE="partitions")$ans
}

"todec" <- function(bin){
  sum(bin * 2^(rev(seq_along(bin)-1)))
}

"comptobin" <- function(comp, check=TRUE){
  if(check){
    stopifnot(all(comp==round(comp)))
  }
  comp <- comp[comp>0]
  s <- sum(comp)
  .C("c_comptobin",
     as.integer(comp),
     as.integer(length(comp)),
     ans=integer(s),
     PACKAGE = "partitions")$ans[-s]
}

"bintocomp" <- function(bin, use.C=TRUE, check=TRUE){
  if(use.C){
    .C("c_bintocomp",
       as.integer(bin),
       as.integer(length(bin)),
       ans=as.integer(rep(1,1+sum(bin != 0))),
       PACKAGE="partitions")$ans
  } else {
    diff(c(0L,which(c(bin,1L)>0)))
  }
}

"plainperms" <- function(n){
  fn <- factorial(n)
  kk <- integer(n*fn)

  out <- .C("c_plainperms",
            ans = kk,
            as.integer(n),
            as.integer(fn),
            PACKAGE="partitions"
            )$ans

  dim(out) <- c(n,fn)
  return(as.partition(out))
}

`mset` <-  function(v){
  v <- sort(v)
  stopifnot(all(v==round(v)))
  n <- length(v)
  nn <- round(exp(lfactorial(n)-sum(lfactorial(table(v)))))

  out <- .C("c_allperms",
            as.integer(v),
            as.integer(n),
            as.integer(nn),
            ans = integer(n*nn),
            PACKAGE="partitions"
            )$ans
  dim(out) <- c(n,nn)
  return(as.partition(out))
}

`multiset` <- function(v,n=length(v)){
  v <- sort(v)
  if(n==length(v)){return(mset(v))} # unnecessary, function works if this line is commented out
  if(n==1){return(as.partition(rbind(sort(unique(v)))))}
  m <- blockparts(table(v),n)
  m <- lapply(split(m,col(m)),function(u){rep(unique(v),u)})
  m <- lapply(m,mset)
  as.partition(do.call("cbind",m))
}

`vec_to_set` <- function(vec){
    jj <- sort(unique(vec))
    M <- outer(jj,vec,`==`)
    out <- lapply(split(M, seq_len(nrow(M))),which)
    names(out) <- NULL
    return(do.call(sets::set,lapply(out,sets::as.set)))
}

`vec_to_eq` <- function(vec){
    out <- split(seq_along(vec),vec)
    class(out) <- c(class(out),"equivalence")
    return(out)
  }

`multinomial` <- function(v){
    jj <- rep(seq_along(v),v)
    out <- as.partition(apply(multiset(jj),2,order))
    rownames(out) <- rep(names(v),v)
    return(out)
}

`allbinom` <- function(n,k){as.partition(multinomial(c(k,n-k))[seq_len(k),,drop=FALSE])}

`genrif` <- function(v){

  f <- function(x){
    out <- x
    n <- 0
    for(i in seq_along(v)){
      out[x==i] <- n + seq_len(v[i])
      n <- n + v[i]
    }
    return(out)
  }
  as.partition(apply(multiset(rep(seq_along(v),times=v)),2,f))
}

`riffle` <- function(p,q=p){genrif(c(p,q))}

