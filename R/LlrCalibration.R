# Copyright 2021 Observational Health Data Sciences and Informatics
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

#' Calibrate the log likelihood ratio
#'
#' @description
#' \code{calibrateP} computes calibrated p-values using the fitted null distribution
#'
#' @details
#' This function computes a calibrated two-sided p-value as described in Schuemie et al (2014).
#'
#' @param likelihoodApproximation    A data frame containing either normal, skew-normal, custom parametric, or grid
#'                  likelihood data. 
#' @param null      An object of class \code{null} created using the \code{fitNull} function or an
#'                  object of class \code{mcmcNull} created using the \code{fitMcmcNull} function.
#' @param twoSided  Compute two-sided (TRUE) or one-sided (FALSE) p-value?
#' @param upper     If one-sided: compute p-value for upper (TRUE) or lower (FALSE) bound?
#'
#' @return
#' The calibrated og likelihood ratio.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' calibrateLlr(null, positive)
#'

#' @export
calibrateLlr <- function(null, likelihoodApproximation, twoSided = FALSE, upper = TRUE) {
  UseMethod("calibrateLlr")
}


#' @describeIn
#' calibrateLlr Computes the calibrated log likelihood ratio using asymptotic assumptions.
#' 
#' @export
calibrateLlr.null <- function(null, likelihoodApproximation, twoSided = FALSE, upper = TRUE) {
  if (twoSided || !upper)
    stop("Currently only one-sided upper LLRs are supported")
  
  if ("logRr" %in% colnames(likelihoodApproximation)) {
    message("Detected data following normal distribution")
    logLikelihood <- function(x, row) {
      dnorm(x, mean = row$logRr, sd = row$seLogRr, log = TRUE)
    }
  } else if ("gamma" %in% colnames(likelihoodApproximation)) {
    message("Detected data following custom parameric distribution")
    logLikelihood <- function(x, row) {
      ((exp(row$gamma * (x - row$mu)))) * ((-(x - row$mu)^2)/(2 * row$sigma^2))
    }
  } else if ("alpha" %in% colnames(likelihoodApproximation)) {
    message("Detected data following skew normal distribution")
    logLikelihood <- function(x, row) {
      if (is.infinite(row$sigma)) {
        return(rep(0, length(x)))
      }
      return(log(2) + dnorm(x, row$mu, row$sigma, log = TRUE) + pnorm(row$alpha * (x - row$mu), 0, row$sigma, log.p = TRUE))
    }
  } else {
    message("Detected data following grid distribution")
    # TODO: implement efficient grid likelihood function
  }
  calibrateOneLlr <- function(i) {
    optimum <-  suppressWarnings(optim(0, function(x) -logLikelihood(x, row = likelihoodApproximation[i, ])))
    mle <- optimum$par
    ml <- -optimum$value
    if (null[2] < 0.001) {
      if (mle < null[1]) {
        calibratedLlr <- 0
      } else {
        calibratedLlr <- ml - logLikelihood(null[1], likelihoodApproximation[i, ])
      }
    } else {
      calibratedP <- integrate(f = computePAtNull, 
                               logLikelihood = logLikelihood, 
                               row = likelihoodApproximation[i, ],
                               mle = mle,
                               ml = ml,
                               null = null,
                               lower = null[1] - 10 * null[2], 
                               upper = null[1] + 10 * null[2])$value
      calibratedLlr <- computeLlrFromP(calibratedP)
    }
    
    return(calibratedLlr)
  }
  
  calibratedLlr <- sapply(1:nrow(likelihoodApproximation), calibrateOneLlr)
  names(calibratedLlr) <- NULL
  return(calibratedLlr)
}

# gridLikelihood <- function(x, row) {
#   gridX <- as.numeric(colnames(row))
#   if (any(is.na(gridX))) {
#     stop("Expecting grid data, but not all column names are numeric")
#   }
#   lowIdx <- 
#   
# }

computePFromLlr <- function(llr, mle, nullLogRr = 0) {
  ifelse(mle < nullLogRr, 1 - (1 - pchisq(2 * llr, df = 1))/2, (1 - pchisq(2 * llr, df = 1))/2)
}

computeLlrFromP <- function(p) {
  if (p > 0.5) {
    return(0)
  } else {
    return(qchisq(1 - 2 * p, df = 1) / 2)
  }
}

computePAtNull <- function(x, logLikelihood, row, mle, ml, null) {
  llr <- ml - logLikelihood(x, row)
  p <- computePFromLlr(llr, mle = mle, nullLogRr = x)
  p * dnorm(x, null[1], null[2])
}
