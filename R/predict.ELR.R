#' predict.ELR 
#' 
#' Make predictions using a fitted ELR model
#' 
#'  @name predict.ELR
#'  
#' @param object Fitted model from \code{\link{ELR}}.
#' @param newdata A new input data frame.
#' @param type of prediction: "prop" for probabilities and "class" for 
#' class labels.
#' @param ... Further arguments passed to or from other methods.
#' @return A list with items 
#' \item{prob}{predicted class probabilities}
#' \item{class}{predicted class memberships obtained by thresholding 
#' class probabilities at the prevalence rate of the positive class} 
#'
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#'
NULL
#' 
#'  @name predict.ELR
#' @export 
predict.ELR <- function(object, newdata, type = c("prob", "class")[1], ...){
if (!inherits(object, "ELR")) stop("Object must be a \"ELR \"'")
beta = object$beta 
b = object$b 
p = object$p  
if(is.null(object$form))
newdata <- data.matrix(newdata)    
else 
newdata <- data.matrix(newdata[, rhs.form(object$form), drop = FALSE])
W <- Weights.ELR(dim(newdata)[2], p, object$ken, seed = object$seed)
Q = Kernel.ELR(newdata, p, W$W, object$ken)
p = Q%*%beta + b 
prob = 1/(1 + exp(-p)) 
if(type == "class") 
pred = ifelse(prob >= object$prior, 1, 0)
else if(type == "prob") {
pred =  cbind(1-prob, prob)
colnames(pred) = c("0", "1")
} else stop("type unknown") 
return(pred)
}
#'  @name predict.ELR
#' @export 
predict.ELM <- function(object, newdata, ...){
  if (!inherits(object, "ELM")) stop("Object must be a \"ELM \"'")
  beta = object$beta 
  b = object$b 
  p = object$p 
  
  if(is.null(object$form))
  newdata <- data.matrix(newdata)    
  else 
  newdata <- data.matrix(newdata[, rhs.form(object$form), drop = FALSE])
  
  W <- Weights.ELR(dim(newdata)[2], p, object$ken, seed = object$seed)
  Q = Kernel.ELR(newdata, p, W$W, object$ken)
  pred = Q%*%beta + b 
  return(pred)
}

         

