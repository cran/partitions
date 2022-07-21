### R code from vignette source 'partitionspaper.Rnw'

###################################################
### code chunk number 1: partitionspaper.Rnw:123-123
###################################################



###################################################
### code chunk number 2: partitionspaper.Rnw:124-125
###################################################
require(partitions)


###################################################
### code chunk number 3: partitionspaper.Rnw:133-134
###################################################
P(100)


###################################################
### code chunk number 4: partitionspaper.Rnw:138-139
###################################################
diffparts(10)


###################################################
### code chunk number 5: partitionspaper.Rnw:171-177
###################################################
f <- function(x){prod(factorial(x),factorial(tabulate(x)))}
prob <- function(a,n){
  jj <- restrictedparts(n,a,include.zero=FALSE)
  N <- factorial(a)*factorial(n)*sum(1/apply(jj,2,f))
  return(N/a^n)
}


