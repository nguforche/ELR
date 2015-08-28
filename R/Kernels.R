#' Functions to compute randomized feature space kernels.
#' 
#' Though these functions are exported, the user will not necessary use them directly. 
#' They are used internally by other functions and packages 
#'
#' @name kernels 
#' @param dat  A matrix/data frame 
#' @param d,p Dimension of data/random feature space 
#' @param ken kernel function 
#' @param seed The value of .Random.seed for which ELM was trained. Set to NULL during training    
#' @param  W Random Weight matrix as computed by the function \code{\link{Weights.ELR}}
#' @param nCols,nRows number of columns and rows of random weight matrix 
#' @param minValue,maxValue minimum and maximum value for uniform distribution
#' @param ratio mixture proportion for hybrid kernel  
#' @return random feature space kernel matrix or random weight matrix  
NULL 
#' @rdname kernels
#' @export 
sigmoid.kernel.ELR <- function(dat, p, W){
dat <- cbind(dat, b = rep(1, nrow(dat)))
Q <- dat%*%t(W)
Q[] <- sigmoid(Q)
return(Q)
}
#' @rdname kernels 
#' @export 
wavelet.hyperbolic.sine.kernel <- function(dat, p, W, ratio=0.3){
a <- W$w
b = W$b 
b <- matrix(rep(b, nrow(dat)), nrow = p, ncol = nrow(dat), byrow = F)
Q <- dat%*%t(a) + t(b)
Q[] <- ratio*invHypSine(Q) + (1-ratio)*wavelet(Q) 
return(Q)
}
#' @rdname kernels 
#' @export 
rbf <- function(dat, p, W){
w <- W$w
b = W$b 
Q =  sapply(1:p, function(yy) apply(dat, 1, function(x) (b[yy]*sum((x-w[yy, ])^2))))
Q[] <- exp(-Q)
return(Q)
}
#' @rdname kernels
#' @export 
HypTan <- function(dat, p, W){
dat <- cbind(dat, b = rep(1, nrow(dat)))
Q <- dat%*%t(W)
Q[] <- tanh(Q)
return(Q)
}
#' @rdname kernels 
#' @export 
Fourier <- function(dat, p, W){
dat <- cbind(dat, b = rep(1, nrow(dat)))
Q <- dat%*%t(W)
Q[] <- cos(Q)
return(Q)
}
#' @rdname kernels 
#' @export 
HardLimit <- function(dat, p, W){
dat <- cbind(dat, b = rep(1, nrow(dat)))
Q <- dat%*%t(W)
Q <- (Q <= 0) + 0 
return(Q)
}
#' @rdname kernels  
#' @export 
MultiQuad <- function(dat, p, W){
w <- W$w
b = W$b 
Q =  sapply(1:p, function(yy) apply(dat, 1, function(x) 
( sqrt(sqrt(sum((x-w[yy, ])^2)) + b[yy]^2) )))
Q[] <- exp(-Q)
return(Q)
}
#' @rdname kernels  
#' @export 
Kernel.ELR <- function(dat, p, W, ken="sigmoid"){
switch(ken, 
sigmoid = {Q <- sigmoid.kernel.ELR(dat, p, W)},
HypTan = {Q <- HypTan(dat, p, W)},
Fourier = {Q <- Fourier(dat, p, W)},
HardLimit = {Q <- HardLimit(dat, p, W)},
rbf = {Q <- rbf(dat, p, W)},
MultiQuad = {Q <- MultiQuad(dat, p, W)},
hybrid = {Q <- wavelet.hyperbolic.sine.kernel(dat, p, W)},
stop("Wrong Kernel Choice"))
return(Q)
}
#' @rdname kernels 
#' @export 
Weights.ELR <- function(d, p, ken, seed = NULL){
if (exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE))
     temp <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)  ### get existing random seed 
else temp<- NULL

if(is.null(seed)){    
    runif(1)  ## generates new seed 
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) ## get it 
} else assign(".Random.seed",  value  = seed, envir=.GlobalEnv)  ## assign seed from a previous session 
   
switch(ken,
sigmoid = {     
	W.s <- replicate(d, runif(p, -1/sqrt(d), 1/sqrt(d)))
	b.s <- apply(W.s, 1, function(x) runif(1, min(x), abs(max(x))))
	W = cbind(W.s, b.s)},
HypTan = {     
	W.s <- replicate(d, runif(p, -1/sqrt(d), 1/sqrt(d)))
	b.s <- apply(W.s, 1, function(x) runif(1, min(x), abs(max(x))))
	W = cbind(W.s, b.s)},
Fourier = {   
	W.s <- replicate(d, runif(p, -1/sqrt(d), 1/sqrt(d)))
	b.s <- apply(W.s, 1, function(x) runif(1, min(x), abs(max(x))))
	W = cbind(W.s, b.s)},
	
HardLimit = {      
	W.s <- replicate(d, runif(p, -1/sqrt(d), 1/sqrt(d)))
	b.s <- apply(W.s, 1, function(x) runif(1, min(x), abs(max(x))))
	W = cbind(W.s, b.s)},
	
rbf = {
        w <- randomMatrix.ELR(d, p, -1, 1)
	b <- runif(p, min = 0, max = 1)
	W = list(w=w, b=b)},
MultiQuad = {
        w <- randomMatrix.ELR(d, p, -1, 1)
	b <- runif(p, min = 0, max = 1)
	W = list(w=w, b=b)},
	
hybrid = {
        w <- randomMatrix.ELR(d, p, -1, 1)
	b <- runif(p, min = -1, max = 1)
	W = list(w=w, b=b)},
	
	
stop("Wrong Kernel Choice"))

if (!is.null(temp)) assign(".Random.seed", temp, envir=.GlobalEnv)  ## restore the old random seed 
else if(!is.null(seed)) rm(.Random.seed, pos=1)

res <- list(W = W, seed = seed) 
return(res)
	}
#' @rdname kernels 
#' @export 
randomMatrix.ELR <- function (nCols, nRows, minValue, maxValue) 
{
    myMat <- matrix(runif(nCols * nRows, min = minValue, max = maxValue), 
        ncol = nCols)
    myMat
}



























