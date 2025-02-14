# @file EmpiricalCalibrationUsingMcmc.R
#
# Copyright 2025 Observational Health Data Sciences and Informatics
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

proposalFunction <- function(param, scale) {
  dim <- length(param)
  draw <- rnorm(dim, mean = param, sd = scale)
  
  # Precision cannot be negative:
  draw[2] <- abs(draw[2])
  # draw[2] <- max(0, draw[2])
  return(draw)
}

runMetropolisMcmc <- function(startValue, iterations, scale, logRr, seLogRr) {
  dim <- length(startValue)
  chain <- array(dim = c(iterations + 1, dim))
  logLik <- array(dim = c(iterations + 1, 1))
  acc <- array(dim = c(iterations + 1, 1))
  
  logLik[1] <- -logLikelihoodNullMcmc(startValue, logRr, seLogRr)
  chain[1, ] <- c(startValue)
  acc[1] <- 1
  
  for (i in 1:iterations) {
    # print(paste('itr =', i))
    proposal <- proposalFunction(chain[i, ], scale = scale)
    newLogLik <- tryCatch(-logLikelihoodNullMcmc(proposal, logRr, seLogRr), error = function(e) {
      -1e+10
    })
    
    # print(paste(paste(proposal, collapse = ","), newLogLik))
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

binarySearchMu <- function(modeMu,
                           modeSigma,
                           alpha = 0.1,
                           logRrNegatives = logRrNegatives,
                           seLogRrNegatives = seLogRrNegatives,
                           precision = 1e-07) {
  q <- qchisq(1 - alpha, 1) / 2
  L <- modeMu
  H <- 10
  llMode <- -logLikelihoodNullMcmc(c(modeMu, modeSigma), logRr = logRrNegatives, seLogRr = seLogRrNegatives)
  while (H >= L) {
    M <- L + (H - L) / 2
    llM <- -logLikelihoodNullMcmc(c(M, modeSigma), logRr = logRrNegatives, seLogRr = seLogRrNegatives)
    metric <- llMode - llM - q
    # writeLines(paste('M =', M, 'Metric = ',metric))
    if (metric > precision) {
      H <- M
    } else if (-metric > precision) {
      L <- M
    } else {
      return(abs(M - modeMu))
    }
    if (M == modeMu || M == 10) {
      return(0)
    }
  }
}

binarySearchSigma <- function(modeMu,
                              modeSigma,
                              alpha = 0.1,
                              logRrNegatives = logRrNegatives,
                              seLogRrNegatives = seLogRrNegatives,
                              precision = 1e-07) {
  q <- qchisq(1 - alpha, 1) / 2
  llMode <- -logLikelihoodNullMcmc(c(modeMu, modeSigma), logRr = logRrNegatives, seLogRr = seLogRrNegatives)
  L <- modeSigma
  for (i in 1:10) {
    H <- modeSigma + exp(i)
    llM <- -logLikelihoodNullMcmc(c(modeMu, H), logRr = logRrNegatives, seLogRr = seLogRrNegatives)
    metric <- llMode - llM - q
    if (metric > 0) {
      break
    }
  }
  if (i == 10) {
    return(0)
  }
  while (H >= L) {
    M <- L + (H - L) / 2
    llM <- -logLikelihoodNullMcmc(c(modeMu, M), logRr = logRrNegatives, seLogRr = seLogRrNegatives)
    metric <- llMode - llM - q
    # writeLines(paste('M =', M, 'Metric = ',metric))
    if (metric > precision) {
      H <- M
    } else if (-metric > precision) {
      L <- M
    } else {
      return(abs(M - modeSigma))
    }
    if (M == modeSigma) {
      return(0)
    }
  }
}

#' Fit the null distribution using MCMC
#'
#' @description
#' \code{fitNull} fits the null distribution to a set of negative controls using Markov Chain Monte
#' Carlo (MCMC).
#'
#' @details
#' This is an experimental function for computing the 95 percent credible interval of a calibrated
#' p-value using Markov-Chain Monte Carlo (MCMC).
#'
#' @param logRr     A numeric vector of effect estimates on the log scale
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#' @param iter      Number of iterations of the MCMC.
#'
#' @return
#' An object of type \code{mcmcNull} containing the mean and standard deviation (both on the log
#' scale) of the null distribution, as well as the MCMC trace.
#'
#' @examples
#' \dontrun{
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitMcmcNull(negatives$logRr, negatives$seLogRr)
#' null
#' plotMcmcTrace(null)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' calibrateP(null, positive$logRr, positive$seLogRr)
#' }
#' @export
fitMcmcNull <- function(logRr, seLogRr, iter = 100000) {
  if (any(is.infinite(seLogRr))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting null distribution")
    logRr <- logRr[!is.infinite(seLogRr)]
    seLogRr <- seLogRr[!is.infinite(seLogRr)]
  }
  if (any(is.infinite(logRr))) {
    warning("Estimate(s) with infinite logRr detected. Removing before fitting null distribution")
    seLogRr <- seLogRr[!is.infinite(logRr)]
    logRr <- logRr[!is.infinite(logRr)]
  }
  if (any(is.na(seLogRr))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting null distribution")
    logRr <- logRr[!is.na(seLogRr)]
    seLogRr <- seLogRr[!is.na(seLogRr)]
  }
  if (any(is.na(logRr))) {
    warning("Estimate(s) with NA logRr detected. Removing before fitting null distribution")
    seLogRr <- seLogRr[!is.na(logRr)]
    logRr <- logRr[!is.na(logRr)]
  }
  if (any(abs(logRr) > log(100))) {
    warning("Estimate(s) with extreme logRr detected: abs(logRr) > log(100). Removing before fitting null distribution")
    seLogRr <- seLogRr[abs(logRr) <= log(100)]
    logRr <- logRr[abs(logRr) <= log(100)]
  }
  if (length(logRr) == 0) {
    warning("No valid estimates left. Returning undefined null distribution")
    mcmc <- list(chain = matrix(c(NA,NA), nrow = 1, ncol = 2))
  } else {
    fit <- optim(c(0, 1), logLikelihoodNullMcmc, logRr = logRr, seLogRr = seLogRr)
    
    # Profile likelihood for roughly correct scale:
    scale <- binarySearchMu(fit$par[1],
                            fit$par[2],
                            logRrNegatives = logRr,
                            seLogRrNegatives = seLogRr
    )
    scale <- c(scale, binarySearchSigma(fit$par[1],
                                        fit$par[2],
                                        logRrNegatives = logRr,
                                        seLogRrNegatives = seLogRr
    ))
    
    # writeLines(paste('Scale:', paste(scale,collapse=',')))
    mcmc <- runMetropolisMcmc(fit$par, iterations = iter, scale, logRr, seLogRr)
  }
  result <- c(median(mcmc$chain[, 1]), median(mcmc$chain[, 2]))
  attr(result, "mcmc") <- mcmc
  class(result) <- "mcmcNull"
  return(result)
}

#' @export
print.mcmcNull <- function(x, ...) {
  writeLines("Estimated null distribution (using MCMC)\n")
  if (is.na(x[1])) {
    writeLines("Undefined")
  } else {
    
    mcmc <- attr(x, "mcmc")
    lb95Mean <- quantile(mcmc$chain[, 1], 0.025)
    ub95Mean <- quantile(mcmc$chain[, 1], 0.975)
    lb95Precision <- quantile(mcmc$chain[, 2], 0.025)
    ub95Precision <- quantile(mcmc$chain[, 2], 0.975)
    output <- data.frame(
      Estimate = c(x[1], x[2]),
      lb95 = c(lb95Mean, lb95Precision),
      ub95 = c(ub95Mean, ub95Precision)
    )
    colnames(output) <- c("Estimate", "lower .95", "upper .95")
    rownames(output) <- c("Mean", "Precision")
    printCoefmat(output)
    writeLines(paste("\nAcceptance rate:", mean(mcmc$acc)))
  }
}

#' @describeIn
#' calibrateP Computes the calibrated P-value and 95 percent credible interval using Markov Chain
#' Monte Carlo (MCMC).
#'
#' @param pValueOnly   If true, will return only the calibrated P-value itself, not the credible
#'                     interval.
#'
#' @export
calibrateP.mcmcNull <- function(null, logRr, seLogRr, twoSided = TRUE, upper = TRUE, pValueOnly, ...) {
  mcmc <- attr(null, "mcmc")
  adjustedP <- data.frame(p = rep(as.numeric(NA), length(logRr)), lb95ci = as.numeric(NA), ub95ci = as.numeric(NA))
  if (!is.na(null[1])) {
    for (i in 1:length(logRr)) {
      if (is.na(logRr[i]) || is.infinite(logRr[i]) || is.na(seLogRr[i]) || is.infinite(seLogRr[i])) {
        adjustedP$p[i] <- NA
        adjustedP$lb95ci[i] <- NA
        adjustedP$ub95ci[i] <- NA
      } else {
        pUpperBound <- pnorm((mcmc$chain[, 1] - logRr[i]) / sqrt((1 / sqrt(mcmc$chain[, 2]))^2 + seLogRr[i]^2))
        pLowerBound <- pnorm((logRr[i] - mcmc$chain[, 1]) / sqrt((1 / sqrt(mcmc$chain[, 2]))^2 + seLogRr[i]^2))
        if (twoSided) {
          p <- pUpperBound
          p[pLowerBound < p] <- pLowerBound[pLowerBound < p]
          p <- p * 2
        } else if (upper) {
          p <- pUpperBound
        } else {
          p <- pLowerBound
        }
        adjustedP$p[i] <- quantile(p, 0.5)
        adjustedP$lb95ci[i] <- quantile(p, 0.025)
        adjustedP$ub95ci[i] <- quantile(p, 0.975)
      }
    }
  }
  if (missing(pValueOnly) || pValueOnly == FALSE) {
    attr(adjustedP, "mcmc") <- mcmc
    return(adjustedP)
  } else {
    return(adjustedP$p)
  }
}
