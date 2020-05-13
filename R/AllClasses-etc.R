### Really want  S4 methods for  all the binary operations.
### Otherwise we "never" make
###       <bigz> o  <bigq>
### or    <bigq> o  <bigz>    working --
##
## But unfortunately the above seems "impossible", see
##  see also  setMethod() in ./matrix-prods.R

## OTOH: This *still* helps to define single-dispatch methods for  asNumeric() :
##       {why does it work there ??}

setOldClass("bigz")#, prototype=as.bigz(integer()))
setOldClass("bigq")#, prototype=as.bigq(integer()))
##                cannot use as.bigz() yet which is only defined in ./bigz.R

## diff() method for these: this is just base::diff.default()
## ---- with 2 lines commented out: '##>>'
.diff.big <- function(x, lag = 1L, differences = 1L, ...)
{
    ismat <- is.matrix(x)
    xlen <- if(ismat) dim(x)[1L] else length(x)
    if (length(lag) != 1L || length(differences) > 1L ||
        lag < 1L || differences < 1L)
	stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen)
	return(x[0L]) # empty, but of proper mode
    r <- x ##>>  r <- unclass(x)  # don't want class-specific subset methods
    i1 <- -seq_len(lag)
    if (ismat)
	for (i in seq_len(differences))
	    r <- r[i1, , drop = FALSE] -
                r[-nrow(r):-(nrow(r)-lag+1L), , drop = FALSE]
    else
        for (i in seq_len(differences))
            r <- r[i1] - r[-length(r):-(length(r)-lag+1L)]
 ##>>  class(r) <- oldClass(x)
    r
}
##--> and entries in ../NAMESPACE
