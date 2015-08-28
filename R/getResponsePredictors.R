#' These functions extract the response or predictors of a formular 
#' 
#' @name getLHS.RHS
#' @param formula R formula 
#' @return response or predictor variables 
NULL 
#' @rdname  getLHS.RHS
#' @examples 
#' form <- as.formula("y ~ V1 + V2")
#' resp <- lhs.form(form)
#'  
#' @export
lhs.form <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
    response <- attr(tt, "response") # index of response var
    vars[response] 
}
#' @rdname  getLHS.RHS
#' @examples 
#' rhs.vars <- rhs.form(form)
#' @export
rhs.form <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
    response <- attr(tt, "response") # index of response var
    vars[-response] 
}

