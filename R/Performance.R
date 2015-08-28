#' Performance.
#' 
#' Compute performance measures for classification/regression models: opt.thresh computes 
#' optimal classification threshold; Eval computes several performance measures 
#' from confusion matrix; Performance similarly computes performance measures 
#'  using optimal or a specified classification threshold.  
#'
#' @name Performance 
#' @param pred,prob A numeric predicted response vector or class decisions 
#' scores or class probabilities
#' @param obs True class; a factor or binary 0,1 numeric vector  
#' @param opt.methods what methods should be used to optimize thresholds.
#' commonly used methods include:
#' \itemize{
#' \item{opt.methods = 9: }{minimum distance between the point (0,1) 
#' and the ROC curve. Default}
#' \item{ opt.methods = 3 : }{maximum sensitivity + specificity} 
#' \item{opt.methods = 2 : }{sensitivity = specificity}
#' }
#' @param prevalence  prevalence rate of the positive class 
#' @param threshold classification threshold 
#' See \code{\link[PresenceAbsence]{optimal.thresholds}}
#' @param lambda numeric weight for weighted accuracy 
#' @param ... Further arguments passed to or from other methods.
#' @return  A numeric or character vector depending on the function 
NULL 
#' @importFrom PresenceAbsence optimal.thresholds
#' @importFrom PresenceAbsence presence.absence.accuracy
#' @importFrom PresenceAbsence predicted.prevalence
#' @importFrom caret confusionMatrix
#' @importFrom DMwR regr.eval
#'
#' @rdname  Performance
#' @export 
opt.thresh <- function(prob, obs, opt.methods = 9){
thresh = 0.5 
if(length(unique(obs)) > 1){
obs <- as.numeric(as.factor(obs))-1 
SIMDATA = cbind.data.frame(plotID = 1:length(obs), Observed = obs, Predicted = prob)
thresh <- optimal.thresholds(SIMDATA, threshold = 101, which.model = 1, opt.methods = opt.methods)
thresh <- ifelse(length(thresh["Predicted"]) >= 1,as.numeric(thresh["Predicted"]), 0.5)
}
return(thresh)
}
#' @rdname  Performance 
#' @export
Eval <- function(prob, obs, lambda = 0.5){
obs = factor(obs, levels = c("1", "0"))
pred = factor(prob,  levels = c("1", "0"))
tab = ftable(obs = obs, pred = pred)
a = tab[1,1];  b = tab[1,2]
c = tab[2,1];  d = tab[2,2]
TNr = ifelse(d+c > 0, d/(d+c), NA)
TPr = ifelse(a+b > 0, a/(a+b), NA)
G.mean = sqrt(TPr*TNr)
W.acc = lambda*TPr + (1-lambda)*TNr 
acc = (a+d)/sum(tab)
P = ifelse(a+c > 0, a/(a+c), NA)
R = ifelse(a+b > 0, a/(a+b), NA)
F = 2*P*R/(P+R)
res = cbind.data.frame(F1= F, Recall = R, Precision = P, Spec = TNr, acc = acc, W.acc = W.acc, G.mean = G.mean)
res 
}
#' @rdname  Performance 
#' @export
Performance <- function(pred, obs, prevalence = NULL, threshold=NULL, ...){
if(length(unique(obs)) == 2) {  ## classification task 
obs <- as.numeric(factor(obs))-1 
## get best cut-off 
if(is.null(threshold))
threshold <- opt.thresh(pred, obs, ...)
### get these performance measures
nme = c("PCC", "PCC.sd", "AUC", "AUC.sd", "sensitivity", "sensitivity.sd", 
"specificity", "specificity.sd")
xx = cbind.data.frame(plotID = 1:length(pred), Observed = obs, Predicted = pred)
accuracy <- presence.absence.accuracy(xx, threshold = threshold, st.dev = TRUE)[, nme]
pred.prev <- predicted.prevalence(DATA=xx, threshold = threshold)[, 
c("Obs.Prevalence", "Predicted")]
nme <- c("Pos Pred Value", "Neg Pred Value", "Balanced Accuracy")
if(is.null(prevalence)) prevalence = as.numeric(pred.prev$Obs.Prevalence)
obs <- factor(ifelse(obs == 1, "Yes", "No"), levels = c("Yes", "No"))
pred <- factor(ifelse(pred >= threshold, "Yes", "No"), levels = c("Yes", "No"))
cmx <- confusionMatrix(data=pred, reference=obs,  prevalence = prevalence)$byClass[nme]
res <- cbind.data.frame(accuracy, t(cmx), pred.prev, threshold = threshold)
} else 
res = data.frame(t(regr.eval(obs, pred, stats= c("mae","mse","rmse"))))
return(res)
}

 














