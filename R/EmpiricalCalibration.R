# @file EmpiricalCalibration.R
#
# Copyright 2014 Observational Health Data Sciences and Informatics
#
# This file is part of EmpiricalCalibration.R
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Fit the null distribution
#'
#' @description
#' \code{fitNull} fits the null distribution to a set of negative controls
#'
#' @details
#' This function fits a Gaussian function to the negative control estimates as described in Schuemie et al (2014).
#' 
#' 
#' @param logRr     A numeric vector of effect estimates on the log scale
#' @param seLogRr    The standard error of the log of the effect estimates. Hint: often the standard 
#' error = (log(<lower bound 95 percent confidence interval>) - log(<effect estimate>))/qnorm(0.025) 
#'  
#' @return An object of type \code{null} containing the mean and standard deviation (both on the log scale) of the
#' null distribution.
#' 
#' @examples 
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0,]
#' null <- fitNull(negatives$logRr,negatives$seLogRr)
#' 
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why empirical
#' calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
#' 
#' @export
fitNull <- function(logRr,seLogRr){
  
  gaussianProduct<-function(mu1,mu2,sd1,sd2) {
    (2*pi)^(-1/2)*(sd1^2+sd2^2)^(-1/2)*exp(-(mu1-mu2)^2/(2*(sd1^2+sd2^2)))
  }
  
  LL <- function(theta,estimate,se){
    result <- 0
    for (i in 1:length(estimate))
      result = result - log(gaussianProduct(estimate[i], theta[1],se[i],exp(theta[2])))
    if (is.infinite(result))
      result = 99999
    result
  }
  theta <- c(0,0)
  fit <- optim(theta,LL, estimate = logRr, se = seLogRr,hessian=TRUE)
  fisher_info <- solve(fit$hessian)
  prop_sigma <- sqrt(diag(fisher_info))
  null <- fit$par
  null[2] <- exp(null[2])
  names(null) = c("mean","sd")
  attr(null,"LB95CI") <- c(fit$par[1] + qnorm(0.025) * prop_sigma[1],exp(fit$par[2] + qnorm(0.025) * prop_sigma[2]))
  attr(null,"UB95CI") <- c(fit$par[1] + qnorm(0.975) * prop_sigma[1],exp(fit$par[2] + qnorm(0.975) * prop_sigma[2]))
  attr(null,"CovarianceMatrix") <- fisher_info
  class(null) <- "null"
  null
}

#' @export
print.null <- function(x, ...) {
  writeLines("Estimated null distribution\n")
  output <- data.frame(Estimate = c(x[1], x[2]), 
                       lb95 = attr(x,"LB95CI"),
                       ub95 = attr(x,"UB95CI"))
  colnames(output) <- c("Estimate", "lower .95", "upper .95")
  rownames(output) <- c("Mean", "SD")
  printCoefmat(output)
}


#' Calibrate the p-value
#'
#' @description
#' \code{calibrateP} computes calibrated p-values using the fitted null distribution
#'
#' @details
#' This function computes a calibrated two-sided p-value as described in Schuemie et al (2014).
#' 
#' @param logRr     A numeric vector of one or more effect estimates on the log scale
#' @param seLogRr    The standard error of the log of the effect estimates. Hint: often the standard 
#' error = (log(<lower bound 95 percent confidence interval>) - log(<effect estimate>))/qnorm(0.025) 
#' @param null      An object of class \code{null} created using the \code{fitNull} function
#' @param pValueConfidenceInterval    If true, computes the 95 percent confidence interval of the calibrated p-value
#'  
#' @return A two-sided calibrated p-value.
#' 
#' @examples 
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0,]
#' null <- fitNull(negatives$logRr,negatives$seLogRr)
#' positive <- sccs[sccs$groundTruth == 1,]
#' calibrateP(positive$logRr,positive$seLogRr, null)
#' 
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why empirical
#' calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
#' 
#' @export
calibrateP <- function(logRr, seLogRr, null, pValueConfidenceInterval = FALSE){
  
  oneAdjustedP <- function(logRR, se, null){
    P_upper_bound = pnorm((null[1]-logRR)/sqrt(null[2]^2+se^2)) 
    P_lower_bound = pnorm((logRR-null[1])/sqrt(null[2]^2+se^2)) 
    2*min(P_upper_bound,P_lower_bound)
  }
  
  adjustedP <- vector(length=length(logRr))
  for (i in 1:length(logRr))
    adjustedP[i] = oneAdjustedP(logRr[i],seLogRr[i],null)
  
  if (pValueConfidenceInterval){
    adjustedP <- data.frame(p = adjustedP, lb95ci = 0,ub95ci = 0)
    rand <- MASS::mvrnorm(10000,c(null[1], log(null[2])),attr(null,"CovarianceMatrix"))
    for (i in 1:length(logRr)){
      P_upper_bound = pnorm((rand[,1]-logRr[i])/sqrt(exp(rand[,2])^2+seLogRr[i]^2)) 
      P_lower_bound = pnorm((logRr[i]-rand[,1])/sqrt(exp(rand[,2])^2+seLogRr[i]^2)) 
      #Take min:
      p <- P_upper_bound
      p[P_lower_bound < p] <- P_lower_bound[P_lower_bound < p]
      p <- p * 2
      
      adjustedP$lb95ci[i] <- quantile(p,0.025)
      adjustedP$ub95ci[i] <- quantile(p,0.975)
    }
  }
  adjustedP
}

#' Compute the (traditional) p-value
#'
#' @description
#' \code{computeTraditionalP} computes the traditional two-sided p-value based on the log of the 
#' relative risk and the standerd error of the log of the relative risk.
#'
#' @param logRr     A numeric vector of one or more effect estimates on the log scale
#' @param seLogRr    The standard error of the log of the effect estimates. Hint: often the standard 
#' error = (log(<lower bound 95 percent confidence interval>) - log(<effect estimate>))/qnorm(0.025) 
#'
#' @return A two-sided (traditional) p-value.
#' 
#' @examples 
#' data(sccs)
#' positive <- sccs[sccs$groundTruth == 1,]
#' computeTraditionalP(positive$logRr, positive$seLogRr)
#' 
#' @export
computeTraditionalP <- function(logRr, seLogRr) {
  z <- logRr / seLogRr
  return(2*pmin(pnorm(z),1-pnorm(z))) # 2-sided p-value  
}