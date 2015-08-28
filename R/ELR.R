#' ELR   
#' 
#' Train and fit LS-ELR and IRLS-ELR  on the training set. Also performs grid search for parameter 
#' selection.
#' @name ELR
#' @param Y True class; a factor or binary 0,1 numeric vector
#' @param X A matrix of predictors   
#' @param form A formula
#' @param dat Matrix/data frame 
#' @param para A named parameter list with required  entries:
#' @param fit (Logical) make predictions on training set? 
#' \itemize{
#' \item{p: }{Dimension of random feature space }
#' \item{gamma: }{Numeric regularization parameter }
#' \item{ken: }{Type of kernel function}
#' \item{max.iter: }{Maximum number of iterations for IRLS-ELR}
#' \item{tol: }{Convergence tolerance for IRLS-ELR}
#' \item{grid: }{Optional dataframe with possible combinations of p and gamma 
#'  to perform grid search}
#' }
#' @param model Type of model to train: "LS-ELR" or "IRLS-ELR". Default 
#' to "LS-ELR".
#' @param \dots Further arguments passed to or from other methods.
#' @return An object of class \code{ELR}; a list with items 
#' \item{beta}{estimated parameters}
#' \item{p}{Dimension of random feature space}
#' \item{W}{Random weight matrix} 
#' \item{ken}{kernel function}
#' \item{pred}{A list with fitted class memberships and class probabilities 
#' on the training set}
#' \item{para}{A named parameter list} 
#' \item{prior}{prevalence rate of the positive class in the training set}
#'
NULL 
#' @rdname ELR 
#' @export
ELR <- function(X, ...) UseMethod("ELR")
#' @rdname ELR 
#' @export
ELR.default <- function(X, Y, para, model=c("LS-ELR","IRLS-ELR")[1], fit = FALSE, ...){   
Y = matrix(ifelse((as.numeric(factor(Y))-1) == 1, 1, 0))
X = data.matrix(X) 
if(sum(!c("ken","p","gamma")%in%names(para)) > 0) stop("Missing parameter items") 
ken <- para$ken
p <- para$p
gamma <- para$gamma
if(model == "LS-ELR"){
mod <- TrainELR(X, Y, p, gamma, ken)
} else if(model == "IRLS-ELR") {
if(sum(!c("tol","max.iter")%in%names(para)) > 0) stop("Missing parameter items") 
tol = para$tol
max.iter = para$max.iter
tol <- para$tol
mod <- TrainIRLSELR(X, Y, p, gamma, ken, max.iter, tol)
} else stop("wrong model type") 
mod$prior <- sum(Y[,1]==1)/nrow(X)  
if(fit) mod$pred <- predict(mod, X)       
mod$para = para  
class(mod) <- "ELR"
return(mod)
}
#' @rdname ELR 
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' dat <- SynData()
#' ix = sample(nrow(dat), floor(nrow(dat)*0.75))
#' dat.trn = dat[ix, ]
#' dat.tst = dat[-ix, ]
#' form <- as.formula(paste("key ~ ", paste(names(dat)[!names(dat)%in%"key"], collapse = "+")))
#' para <- list( ken = "sigmoid", p = 100, gamma = 10.01)
#' mod <- ELR(form, dat.trn, para, model="LS-ELR")
#' prob <- predict(mod, dat.tst)
#' perf <- Performance(prob[,2], dat.tst$key)
#' 
#' para <- list( ken = "sigmoid", p = 100, gamma = 10.01, tol = 1e-6, max.iter = 100)
#' mod2 <- ELR(form, dat.trn, para, model="IRLS-ELR")
#' prob2 <- predict(mod2, dat.tst)
#' perf2 <- Performance(prob2[,2], dat.tst$key)
#'}
ELR.formula <- function(form, dat, para, model=c("LS-ELR","IRLS-ELR")[1], fit = FALSE, ...){  
resp =  lhs.form(form) 
Y = matrix(ifelse(dat[, resp]==1, 1, 0))
rhs.vars = rhs.form(form)
X = dat[, rhs.vars, drop = FALSE]
X = data.matrix(X) 
ken <- para$ken
p <- para$p
gamma <- para$gamma
if(model == "LS-ELR"){
mod <- TrainELR(X, Y, p, gamma, ken)
} else if(model == "IRLS-ELR") {
tol = para$tol
max.iter = para$max.iter
tol <- para$tol
mod <- TrainIRLSELR(X, Y, p, gamma, ken, max.iter, tol)
} else stop("wrong model type") 
mod$form <- form
mod$prior <- sum(Y[,1]==1)/nrow(dat) 
if(fit) mod$pred <- predict(mod, X)       
mod$para = para  
class(mod) <- "ELR"
return(mod)
}
#' @rdname ELR 
#' @export
ELRgridsearch <- function(X, ...) UseMethod("ELRgridsearch")
#' @rdname ELR 
#' @export
ELRgridsearch.default <- function(X, Y, para, model = c("LS-ELR","IRLS-ELR")[1],...){
ix = sample(nrow(X), floor(nrow(X)*0.75))
X.trn = X[ix, ]
X.val = X[-ix, ]
Y.trn = Y[ix, ]
Y.val =   Y[-ix, ]
grid <- para$grid 
n = nrow(grid) 
g.mean = c()
for(ii in 1:n){
para$p <- grid[ii, "p"]
para$gamma <- grid[ii, "gamma"] 
elr <-  ELR(X.trn, Y.trn, para, model)
pred <- predict(elr, X.val)
threshold <- opt.thresh(pred[, 2], Y.val)
cls <-ifelse(pred[, 2] >= threshold, 1, 0)
g.mean <-  c(g.mean, as.numeric(Eval(cls, Y.val, lambda = threshold)["G.mean"]))
}
ix = which.max(g.mean)
para$p <- grid[ix, "p"]
para$gamma <- grid[ix, "gamma"] 
para$grid <- NULL 
res <- ELR(X,Y, para, model)
return(res) 
}
#' @rdname ELR 
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' dat <- SynData()
#' ix = sample(nrow(dat), floor(nrow(dat)*0.75))
#' dat.trn = dat[ix, ]
#' dat.tst = dat[-ix, ]
#' form <- as.formula(paste("key ~ ", paste(names(dat)[!names(dat)%in%"key"], collapse = "+")))
#' para <- list( ken = "sigmoid", grid = expand.grid(p = c(50, 100), gamma = c(0.05,10.01)))
#' mod <- ELRgridsearch(form, dat.trn, para, model="LS-ELR")
#' pred <- predict(mod, dat.tst)
#' perf <- Performance(pred$prob[,2], dat.tst$key)
#' 
#' para <- list( ken = "sigmoid", p = 100, gamma = 10.01, tol = 1e-6, max.iter = 100)
#' mod2 <- ELR(form, dat.trn, para, model="IRLS-ELR")
#' pred2 <- predict(mod2, dat.tst)
#' perf2 <- Performance(pred2$prob[,2], dat.tst$key)
#'}
ELRgridsearch.formula <- function(form, dat, para, model = c("LS-ELR","IRLS-ELR")[1],...){
resp <- lhs.form(form)
ix = sample(nrow(dat), floor(nrow(dat)*0.75))
dat.trn = dat[ix, ]
dat.val = dat[-ix, ]
y.val =   dat[-ix, resp]
grid <- para$grid 
n = nrow(grid) 
g.mean = c()
for(ii in 1:n){
para$p <- grid[ii, "p"]
para$gamma <- grid[ii, "gamma"] 
elr <-  ELR(form, dat.trn, para, model)
pred <- predict(elr, dat.val[, rhs.form(form)])
threshold <- opt.thresh(pred[, 2], y.val)
cls <-ifelse(pred[, 2] >= threshold, 1, 0)
g.mean <-  c(g.mean, as.numeric(Eval(cls, y.val, lambda = threshold)["G.mean"]))
}
ix = which.max(g.mean)
para$p <- grid[ix, "p"]
para$gamma <- grid[ix, "gamma"] 
para$grid <- NULL 
res <- ELR(form, dat, para, model)
return(res) 
}

#' @rdname ELR 
#' @export
ELM <- function(X, ...) UseMethod("ELM")

#' @rdname ELR 
#' @export
ELM.default <- function(X, Y, para, fit = FALSE, ...){  
  Y = matrix(Y)
  X = data.matrix(X) 
  if(sum(!c("ken","p","gamma")%in%names(para)) > 0) stop("Missing parameter items") 
  ken <- para$ken
  p <- para$p
  gamma <- para$gamma
  mod <- TrainELM(X, Y, p, gamma, ken)  
  if(fit) mod$pred <- predict(mod, X) 
  mod$para = para  
  class(mod) <- "ELM"
  return(mod)
}
#' @rdname ELR 
#' @export
ELM.formula <- function(form, dat, para, fit = FALSE, ...){  
  resp =  lhs.form(form) 
  Y = data.matrix(dat[, resp, drop = FALSE])
  rhs.vars = rhs.form(form)
  X = dat[, rhs.vars, drop = FALSE]
  X = data.matrix(X) 
  ken <- para$ken
  p <- para$p
  gamma <- para$gamma

  mod <- TrainELM(X, Y, p, gamma, ken) 
  mod$form <- form  
  if(fit) mod$pred <- predict(mod, X)
  mod$para = para  
  class(mod) <- "ELM"
  return(mod)
}

#' @rdname ELR 
#' @export
ELMgridsearch <- function(X, ...) UseMethod("ELMgridsearch")
#' @rdname ELR 
#' @export
ELMgridsearch.default <- function(X, Y, para,...){
  ix = sample(nrow(X), floor(nrow(X)*0.75))
  X.trn = data.matrix(X[ix, ,drop=FALSE])
  X.val = data.matrix(X[-ix,,drop = FALSE])
  Y.trn = data.matrix(Y[ix, , drop = FALSE])
  Y.val =   data.matrix(Y[-ix, ,drop = FALSE])
  grid <- para$grid 
  n = nrow(grid) 
  mae = c()
  for(ii in 1:n){
    para$p <- grid[ii, "p"]
    para$gamma <- grid[ii, "gamma"] 
    elr <-  ELM(X.trn, Y.trn, para)
    pred <- predict(elr, X.val)    
    mae = c(mae, as.numeric(regr.eval(trues=Y.val[,1],preds= as.numeric(pred),stats= "mae", 
                                      train.y = Y.trn[,1])))                                                  
    }
  rm(elr)
  ix = which.min(mae)
  para$p <- grid[ix, "p"]
  para$gamma <- grid[ix, "gamma"] 
  para$grid <- NULL 
  res <- ELM(X,Y, para)
  return(res) 
}
#' @rdname ELR 
#' @export
ELMgridsearch.formula <- function(form, dat, para,...){
  resp <- lhs.form(form)
  rhs.vars <- rhs.form(form)
  ix = sample(nrow(dat), floor(nrow(dat)*0.75))
  Y.trn <- data.matrix(dat[ix, resp, drop = FALSE])
  X.trn = data.matrix(dat[ix, rhs.vars, drop=FALSE])
  X.val = data.matrix(dat[-ix, rhs.vars, drop = FALSE])
  Y.val =   data.matrix(dat[-ix, resp, drop = FALSE])
  grid <- para$grid 
  n = nrow(grid) 
  mae = c()
  for(ii in 1:n){
    para$p <- grid[ii, "p"]
    para$gamma <- grid[ii, "gamma"] 
    elr <-  ELM(X.trn, Y.trn, para)
    pred <- predict(elr, X.val)    
    mae = c(mae, as.numeric(regr.eval(trues= Y.val[, 1],preds= as.numeric(pred),stats= "mae", 
                                      train.y =  Y.trn[, 1])))   
                                                     
  }
  rm(elr)
  ix = which.min(mae)  
  para$p <- grid[ix, "p"]
  para$gamma <- grid[ix, "gamma"] 
  para$grid <- NULL 
  res <- ELM(form, dat, para)
  return(res) 
}













