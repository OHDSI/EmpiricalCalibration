# @file EmpiricalCalibrationUsingMcmc.R
#
# Copyright 2015 Observational Health Data Sciences and Informatics
#
# This file is part of EmpiricalCalibration
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

proposalFunction <- function(param, scale = c(0.1, 0.1), fixed = NULL) {
  dim <- length(param)
  draw <- rnorm(dim, mean = param, sd = scale)
  if (!is.null(fixed)) {
    draw[fixed] <- param[fixed]
  }
  draw[2] <- abs(draw[2])
  return(draw)
}

runMetropolisMcmc <- function(startValue, ll, iterations, scale, logRr, seLogRr, fixed = NULL) {
  dim <- length(startValue)
  chain <- array(dim = c(iterations + 1, dim))
  logLik <- array(dim = c(iterations + 1, 1))
  acc <- array(dim = c(iterations + 1, 1))
  
  logLik[1] <- -ll(startValue, logRr, seLogRr)
  chain[1, ] <- c(startValue)
  acc[1] <- 1
  
  for (i in 1:iterations) {
    # print(paste('itr =', i))
    proposal <- proposalFunction(chain[i, ], scale = scale, fixed = fixed)
    
    newLogLik <- tryCatch(-ll(proposal, logRr, seLogRr), error = function(e) {
      -1e+10
    })
    
    prob <- exp(newLogLik - logLik[i])
    if (runif(1) < prob) {
      chain[i + 1, ] <- proposal
      logLik[i + 1] <- newLogLik
      acc[i + 1] <- 1
    } else {
      chain[i + 1, ] <- chain[i, ]
      logLik[i + 1] <- logLik[i]
      acc[i + 1] <- 0
    }
  }
  result <- list(logLik = logLik, chain = chain, acc = acc)
  return(result)
}

gaussianProduct <- function(mu1, mu2, sd1, sd2) {
  (2 * pi)^(-1/2) * (sd1^2 + sd2^2)^(-1/2) * exp(-(mu1 - mu2)^2/(2 * (sd1^2 + sd2^2)))
}

# Note: This function uses precision instead of standard deviation:
logLikelihood <- function(theta, estimate, se) {
  result <- 0
  for (i in 1:length(estimate)) {
    result <- result + log(gaussianProduct(estimate[i], theta[1], se[i], 1/sqrt(theta[2])))
  }
  if (is.infinite(result)) {
    result <- -99999
  }
  # Add weak prior for when precision becomes very large:
  result <- result + dgamma(theta[2], shape = 0.001, rate = 0.001, log = TRUE)
  return(-result)
}

binarySearchMu <- function(modeMu, modeSigma, alpha = 0.1, logRrNegatives = logRrNegatives, seLogRrNegatives = seLogRrNegatives, precision = 0.0000001){
  q <- qchisq(1-alpha, 1)/2
  L <- modeMu
  H <- 10
  llMode <- -logLikelihood(c(modeMu, modeSigma), estimate = logRrNegatives, se = seLogRrNegatives)
  while (H >= L) {
    M <- L + (H - L)/2
    llM <- -logLikelihood(c(M, modeSigma), estimate = logRrNegatives, se = seLogRrNegatives)
    metric <- llMode - llM - q
    #writeLines(paste("M =", M, "Metric = ",metric))
    if (metric > precision){
      H <- M 
    } else if (-metric > precision) {
      L <- M 
    } else {
      return(abs(M - modeMu))
    }
    if (M == modeMu || M == 10)
      return(0)
  }
}

binarySearchSigma <- function(modeMu, modeSigma, alpha = 0.1, logRrNegatives = logRrNegatives, seLogRrNegatives = seLogRrNegatives, precision = 0.0000001){
  q <- qchisq(1-alpha, 1)/2
  L <- modeSigma
  H <- 1500
  llMode <- -logLikelihood(c(modeMu, modeSigma), estimate = logRrNegatives, se = seLogRrNegatives)
  while (H >= L) {
    M <- L + (H - L)/2
    llM <- -logLikelihood(c(modeMu, M), estimate = logRrNegatives, se = seLogRrNegatives)
    metric <- llMode - llM - q
    #writeLines(paste("M =", M, "Metric = ",metric))
    if (metric > precision){
      H <- M 
    } else if (-metric > precision) {
      L <- M 
    } else {
      return(abs(M - modeSigma))
    }
    if (M == modeSigma || M == 1500)
      return(0)
  }
}

#' Fit the null distribution using MCMC
#'
#' @description
#' \code{fitNull} fits the null distribution to a set of negative controls using Markov Chain Monte Carlo (MCMC).
#'
#' @details
#' This is an experimental function for computing the 95 percent confidence interval of a calibrated
#' p-value using Markov-Chain Monte Carlo (MCMC). This should give better estimates than the default
#' function when the standard deviation of the error distribution is close to zero.
#'
#' @param logRr     A numeric vector of effect estimates on the log scale
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#' @param iter               Number of iterations of the MCMC.
#'
#' @return
#' An object of type \code{mcmcNull} containing the mean and standard deviation (both on the log scale) of
#' the null distribution, as well as the MCMC trace.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitMcmcNull(negatives$logRr, negatives$seLogRr)
#' null
#' plotMcmcTrace(null)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' calibrateP(null, positive$logRr, positive$seLogRr, pValueConfidenceInterval = TRUE)
#'
#' @export
fitMcmcNull <- function(logRrNegatives,
                        seLogRrNegatives,
                        logRrPositives,
                        seLogRrPositives,
                        iter = 10000) {
  fit <- optim(c(0, 0.1), logLikelihood, estimate = logRrNegatives, se = seLogRrNegatives)
  
  # Profile likelihood for roughly correct scale:
  scale <- binarySearchMu(fit$par[1], fit$par[2], logRrNegatives = logRrNegatives, seLogRrNegatives = seLogRrNegatives)
  scale <- c(scale,binarySearchSigma(fit$par[1], fit$par[2], logRrNegatives = logRrNegatives, seLogRrNegatives = seLogRrNegatives))
 
  #writeLines(paste("Scale:", paste(scale,collapse=",")))
  mcmc <- runMetropolisMcmc(fit$par, logLikelihood, iterations = iter, scale, logRrNegatives, seLogRrNegatives)
  result <- c(mean(mcmc$chain[, 1]), mean(mcmc$chain[, 2]))
  attr(result,"mcmc") <- mcmc
  class(result) <- "mcmcNull"
  return(result)
}

#' @export
print.mcmcNull <- function(x, ...) {
  writeLines("Estimated null distribution (using MCMC)\n")
  mcmc <-attr(x,"mcmc")  
  lb95Mean <- quantile(mcmc$chain[, 1], 0.025)
  ub95Mean <- quantile(mcmc$chain[, 1], 0.975)
  lb95Precision <- quantile(mcmc$chain[, 2], 0.025)
  ub95Precision <- quantile(mcmc$chain[, 2], 0.975)
  output <- data.frame(Estimate = c(x[1], x[2]), lb95 = c(lb95Mean, lb95Precision) , ub95 = c(ub95Mean, ub95Precision))
  colnames(output) <- c("Estimate", "lower .95", "upper .95")
  rownames(output) <- c("Mean", "Precision")
  printCoefmat(output)
  writeLines(paste("\nAcceptance rate:",mean(mcmc$acc)))
}

#' @export
calibrateP.mcmcNull <- function(null, logRr, seLogRr, pValueConfidenceInterval = FALSE) {
  mcmc <-attr(null,"mcmc")  
  adjustedP <- data.frame(p = rep(1, length(logRr)), lb95ci = 0, ub95ci = 0)
  for (i in 1:length(logRr)) {
    P_upper_bound <- pnorm((mcmc$chain[,
                                       1] - logRr[i])/sqrt((1/sqrt(mcmc$chain[,
                                                                                       2]))^2 + seLogRr[i]^2))
    P_lower_bound <- pnorm((logRr[i] - mcmc$chain[,
                                                           1])/sqrt((1/sqrt(mcmc$chain[, 2]))^2 + seLogRr[i]^2))
    p <- P_upper_bound
    p[P_lower_bound < p] <- P_lower_bound[P_lower_bound < p]
    p <- p * 2
    adjustedP$p[i] <- quantile(p, 0.5)
    adjustedP$lb95ci[i] <- quantile(p, 0.025)
    adjustedP$ub95ci[i] <- quantile(p, 0.975)
  }
  attr(adjustedP, "mcmc") <- mcmc
  if (pValueConfidenceInterval) {
    return(adjustedP)
  } else {
    return(adjustedP$p)
  }
}
