% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Performance.R
\name{Performance}
\alias{Eval}
\alias{Performance}
\alias{opt.thresh}
\title{Performance.}
\usage{
opt.thresh(prob, obs, opt.methods = 9)

Eval(prob, obs, lambda = 0.5)

Performance(pred, obs, prevalence = NULL, threshold = NULL, ...)
}
\arguments{
\item{obs}{True class; a factor or binary 0,1 numeric vector}

\item{opt.methods}{what methods should be used to optimize thresholds.
commonly used methods include:
\itemize{
\item{opt.methods = 9: }{minimum distance between the point (0,1)
and the ROC curve. Default}
\item{ opt.methods = 3 : }{maximum sensitivity + specificity}
\item{opt.methods = 2 : }{sensitivity = specificity}
}}

\item{lambda}{numeric weight for weighted accuracy}

\item{pred,prob}{A numeric predicted response vector or class decisions
scores or class probabilities}

\item{prevalence}{prevalence rate of the positive class}

\item{threshold}{classification threshold
See \code{\link[PresenceAbsence]{optimal.thresholds}}}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
A numeric or character vector depending on the function
}
\description{
Compute performance measures for classification/regression models: opt.thresh computes
optimal classification threshold; Eval computes several performance measures
from confusion matrix; Performance similarly computes performance measures
 using optimal or a specified classification threshold.
}

