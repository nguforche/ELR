#' ELR   
#' 
#'  Train least squares ELR (LS-ELR) and iterative re-weighted least squares ELR 
#'  (IRLS-ELR). 
#' @name Train.ELR
#' @param X Matrix of predictors 
#' @param Y 0,1 binary response matrix 
#' @param p Dimension of random feature space 
#' @param gamma Numeric regularization parameter 
#' @param ken Type of kenel function 
#' @param max.iter Maximum number of iteration 
#' @param tol Convergence tolerance  
#' @return A list with items 
#' \item{beta}{estimated parameters}
#' \item{p}{Dimension of random feature space}
#' \item{W}{Random weight matrix} 
#' \item{ken}{kernel function}
#' @param ... Further arguments passed to or from other methods.
#' @import matrixcalc MASS
NULL 
#' @rdname Train.ELR 
#' @export
TrainELR <- function(X, Y, p, gamma, ken,...){ 
W <- Weights.ELR(ncol(X), p, ken)
Q <- Kernel.ELR(X, p, W$W, ken)
s <- 4.0/gamma
f <- 4.0*t(Q)%*%(Y-0.5) 
B <- apply(Q, 2, sum)
V <- t(Q)%*%Q 
C <- Q%*%ginv(V)
C <- matrix(apply(C, 2, sum), nrow = 1) 
diag(V) <- diag(V) + s  
A.inv <- ginv(V) 
Z = A.inv%*%f
S.i = 1/as.numeric(C%*%A.inv%*%B )
b <- as.numeric(S.i*C%*%Z)
beta = Z - A.inv%*%B*b
res =  list(beta = beta, b = b, p = p, seed = W$seed, ken = ken)
class(res) = "ELR"
return(res)
}
#' @rdname Train.ELR 
#' @export
TrainIRLSELR <- function(X, Y, p, gamma, ken, max.iter, tol){
W <- Weights.ELR(ncol(X), p, ken)
Q <- Kernel.ELR(X, p, W$W, ken)
s <- 1.0/gamma
beta =  matrix(rep(0, ncol(Q)))
b = 0
u.old <- rbind(beta, b)
r.old = 0 
B <- apply(Q, 2, sum)
V <- t(Q)%*%Q + 0.000001
C <- Q%*%ginv(V)
C.1 <- matrix(apply(C, 2, sum), nrow = 1) 

for(ii in 1:max.iter){
w = Q%*%beta + b 
prob = 1/(1 + exp(-w))
sig = as.numeric(prob*(1-prob))
Wj <- s*t(Q)%*%matrix(t(sapply(1:nrow(Q), function(x) (1/sig[x])*C[x, ])), 
         nrow=nrow(Q), ncol= ncol(C))
#diag(V) <- diag(V) + diag(Wj) 
A = V + Wj + 1e-8 
#V = V + 0.0000001 
A.inv <- ginv(A) 
z <- t(Q)%*%(w  + (1/sig)*(Y - prob)) 
Z = A.inv%*%z
S.i = 1/as.numeric(C.1%*%A.inv%*%B )
b <- as.numeric(S.i*C.1%*%Z)
beta = Z - A.inv%*%B*b
u = rbind(beta, b)
r <- as.numeric(sqrt(t(u - u.old)%*%(u-u.old)))
r.err = abs(r - r.old) 
if(!is.nan(r.err) & (r.err <=  1e-30)){
print("restart with elr")
p = 50
gamma = 10.05 
s =  1000
mod <- TrainELR(X, Y, p, gamma, ken)
beta <- mod$beta 
b = mod$b 
W <- Weights.ELR(ncol(X), p, ken, seed = mod$seed)
Q <- Kernel.ELR(X, p, W$W, ken)
V <- t(Q)%*%Q + 0.000001
C <- Q%*%ginv(V)
C.1 <- matrix(apply(C, 2, sum), nrow = 1) 
u.old <- rbind(beta, b)
}
if(is.nan(r)) {
print("restart at random")
p = 50
s <- 1000
W <- Weights.ELR(ncol(X), p, ken)
Q <- Kernel.ELR(X, p, W$W, ken)
B <- apply(Q, 2, sum)
V <- t(Q)%*%Q + 0.000001
C <- Q%*%ginv(V)
C.1 <- matrix(apply(C, 2, sum), nrow = 1) 
beta =  matrix(rep(0, ncol(Q)))
b = 0
u.old <- rbind(beta, b)
} else if(r < tol) break
u.old <- rbind(beta, b)
r.old <-  r 
}
if(r > tol) print("IRWLS did not converge")
res =  list(beta = beta, b = b, p = p, seed = W$seed, ken = ken)
class(res) = "ELR"
return(res)
}

#' @rdname Train.ELR 
#' @export
### ELM regression 
TrainELM <- function(X, Y, p, gamma, ken){ 
  W <- Weights.ELR(ncol(X), p, ken)
  Q <- Kernel.ELR(X, p, W$W, ken)
  s <- 1/gamma  
  f <- t(Q)%*%Y 
  B <- apply(Q, 2, sum)
  V <- t(Q)%*%Q 
  C <- Q%*%ginv(V)
  C <- matrix(apply(C, 2, sum), nrow = 1) 
  diag(V) <- diag(V) + s  
  A.inv <- ginv(V) 
  Z = A.inv%*%f
  S.i = 1/as.numeric(C%*%A.inv%*%B )
  b <- as.numeric(S.i*C%*%Z)
  beta = Z - A.inv%*%B*b
  res =  list(beta = beta, b = b, p = p, seed = W$seed, ken = ken)
  class(res) = "ELM"
  return(res)
}




         
