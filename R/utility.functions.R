#' Utility functions.
#' 
#' These functions are not meant to be directly called by the user. They 
#' are used internally by other functions and packages that may depend on 
#' ELR: sigmoid is the sigmoid function, 
#' tr computes the trace of a matrix, hprod computes the hadamard product 
#'
#' @name utility.functions 
#' @param x,y A numeric or character vector depending on the function 
#' @return  A numeric or character vector depending on the function 
NULL 
#' @rdname  utility.functions 
#' @export 
sigmoid <- function(x) 1.0 / (1 + exp(-x))
#' @rdname  utility.functions 
sigminus <- function(x) 1.0 + exp(-x) 
#' @rdname  utility.functions 
sigplus <- function(x) 1.0 + exp(x)
#' @rdname  utility.functions 
#' @export 
tr <- function(x) matrix.trace(x)
#' @rdname  utility.functions 
#' @export 
hprod <- function(x,y){hadamard.prod(x, y)}
#' @rdname  utility.functions 
#' @export 
wavelet <- function(x) {
hprod(cos(1.5*x), exp(-0.05*x^2))
}
#' @rdname  utility.functions 
#' @export 
invHypSine <- function(x) log(abs(x + sqrt(x^2 +1)))





