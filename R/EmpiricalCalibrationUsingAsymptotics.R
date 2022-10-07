# @file EmpiricalCalibrationUsingAsymptotics.R
#
# Copyright 2022 Observational Health Data Sciences and Informatics
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

#' Fit the null distribution
#'
#' @description
#' \code{fitNull} fits the null distribution to a set of negative controls
#'
#' @details
#' This function fits a Gaussian function to the negative control estimates as described in Schuemie
#' et al (2014).
#'
#'
#' @param logRr     A numeric vector of effect estimates on the log scale
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#'
#' @return
#' An object containing the parameters of the null distribution.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' null
#'
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why
#' empirical calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
#'
#' @export
fitNull <- function(logRr, seLogRr) {
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
  if (length(logRr) == 0) {
    warning("No estimates remaining")
    null <- c(NA, NA)
  } else {
    theta <- c(0, 100)
    fit <- optim(theta, logLikelihoodNull, logRr = logRr, seLogRr = seLogRr)
    null <- fit$par
    null[2] <- 1 / sqrt(null[2])
  }
  names(null) <- c("mean", "sd")
  class(null) <- "null"
  return(null)
}

#' @export
print.null <- function(x, ...) {
  writeLines("Estimated null distribution\n")
  output <- data.frame(Estimate = c(x[1], x[2]))
  colnames(output) <- c("Estimate")
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
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#' @param null      An object of class \code{null} created using the \code{fitNull} function or an
#'                  object of class \code{mcmcNull} created using the \code{fitMcmcNull} function.
#' @param twoSided  Compute two-sided (TRUE) or one-sided (FALSE) p-value?
#' @param upper     If one-sided: compute p-value for upper (TRUE) or lower (FALSE) bound?
#' @param ...       Any additional parameters (currently none).
#'
#' @return
#' The calibrated p-value.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' calibrateP(null, positive$logRr, positive$seLogRr)
#'
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why
#' empirical calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
#'
#' @export
calibrateP <- function(null, logRr, seLogRr, twoSided = TRUE, upper = TRUE, ...) {
  UseMethod("calibrateP")
}


#' @describeIn
#' calibrateP Computes the calibrated P-value using asymptotic assumptions.
#'
#' @export
calibrateP.null <- function(null, logRr, seLogRr, twoSided = TRUE, upper = TRUE, ...) {
  if (length(logRr) != length(seLogRr)) {
    stop("The logRr and seLogRr arguments must be of equal length")
  }
  
  calibrateOneP <- function(i) {
    if (is.na(logRr[i]) || is.infinite(logRr[i]) || is.na(seLogRr[i]) || is.infinite(seLogRr[i])) {
      return(NA)
    } else {
      pUpperBound <- pnorm((null[1] - logRr[i]) / sqrt(null[2]^2 + seLogRr[i]^2))
    }
    pLowerBound <- pnorm((logRr[i] - null[1]) / sqrt(null[2]^2 + seLogRr[i]^2))
    if (twoSided) {
      return(2 * min(pUpperBound, pLowerBound))
    } else if (upper) {
      return(pUpperBound)
    } else {
      return(pLowerBound)
    }
  }
  
  calibratedP <- sapply(1:length(logRr), calibrateOneP)
  names(calibratedP) <- NULL
  return(calibratedP)
}

#' Compute the (traditional) p-value
#'
#' @description
#' \code{computeTraditionalP} computes the traditional two-sided p-value based on the log of the
#' relative risk and the standard error of the log of the relative risk.
#'
#' @param logRr     A numeric vector of one or more effect estimates on the log scale
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#' @param twoSided  Compute two-sided (TRUE) or one-sided (FALSE) p-value?
#' @param upper     If one-sided: compute p-value for upper (TRUE) or lower (FALSE) bound?
#'
#' @return
#' The (traditional) p-value.
#'
#' @examples
#' data(sccs)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' computeTraditionalP(positive$logRr, positive$seLogRr)
#'
#' @export
computeTraditionalP <- function(logRr, seLogRr, twoSided = TRUE, upper = TRUE) {
  z <- logRr / seLogRr
  pUpperBound <- 1 - pnorm(z)
  pLowerBound <- pnorm(z)
  if (twoSided) {
    return(2 * pmin(pUpperBound, pLowerBound))
  } else if (upper) {
    return(pUpperBound)
  } else {
    return(pLowerBound)
  }
}


#' Fit the null distribution using non-normal log-likelihood approximations
#'
#' @description
#' \code{fitNullNonNormalLl} fits the null distribution to a set of negative controls
#'
#' @details
#' This function fits a Gaussian function to the negative control estimates, using non-normal
#' approximations of the per-negative control log likelihood.
#'
#'
#' @param likelihoodApproximations Either a data frame containing normal, skew-normal, or custom parametric likelihood
#'                                 approximations, or a list of (adaptive) grid likelihood profiles.
#'
#' @return
#' An object containing the parameters of the null distribution.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNullNonNormalLl(negatives)
#' null
#'
#' @export
fitNullNonNormalLl <- function(likelihoodApproximations) {
  if (is.data.frame(likelihoodApproximations)) {
    if ("logRr" %in% colnames(likelihoodApproximations)) {
      message("Detected data following normal distribution")
      return(fitNull(logRr = likelihoodApproximations$logRr, seLogRr = likelihoodApproximations$seLogRr))
    } else if ("gamma" %in% colnames(likelihoodApproximations)) {
      message("Detected data following custom parameric distribution")
      type <- "custom"
      llApproximationFunction <- customLlApproximation
      idx <- is.na(likelihoodApproximations$mu) | is.na(likelihoodApproximations$sigma) | is.na(likelihoodApproximations$gamma)
      if (any(idx)) {
        warning("Approximations with NA parameters detected. Removing before fitting null distribution")
        likelihoodApproximations <- likelihoodApproximations[!idx, ]
      }
      likelihoodApproximations <- split(likelihoodApproximations, 1:nrow(likelihoodApproximations))
    } else if ("alpha" %in% colnames(likelihoodApproximations)) {
      message("Detected data following skew normal distribution")
      type <- "skewNormal"
      llApproximationFunction <- skewNormalLlApproximation
      idx <- is.na(likelihoodApproximations$mu) | is.na(likelihoodApproximations$sigma) | is.na(likelihoodApproximations$alpha)
      if (any(idx)) {
        warning("Approximations with NA parameters detected. Removing before fitting null distribution")
        likelihoodApproximations <- likelihoodApproximations[!idx, ]
      }
      likelihoodApproximations <- split(likelihoodApproximations, 1:nrow(likelihoodApproximations))
    } else {
      message("Detected data following grid distribution")
      type <- "grid"
      llApproximationFunction <- gridLlApproximation
      suppressWarnings(point <- as.numeric(colnames(likelihoodApproximations)))
      if (any(is.na(point))) {
        stop("Expecting grid data, but not all column names are numeric")
      }
      convertToDataFrame <- function(i) {
        data.frame(
          value = as.numeric(likelihoodApproximations[i, ]),
          point = point
        )
      }
      likelihoodApproximations <- lapply(1:nrow(likelihoodApproximations), convertToDataFrame)
    }
  } else {
    message("Detected data following grid distribution")
    type <- "grid"
    llApproximationFunction <- gridLlApproximation
    for (i in 1:length(likelihoodApproximations)) {
      likelihoodApproximation <- likelihoodApproximations[[i]]
      maxValue <- max(likelihoodApproximation$value)
      likelihoodApproximation$value <- likelihoodApproximation$value - maxValue
      likelihoodApproximations[[i]] <- as.data.frame(likelihoodApproximation)
    }
  }
  
  theta <- c(0, 100)
  fit <- optim(
    par = theta,
    fn = logLikelihoodNullNonNormalLl,
    likelihoodApproximations = likelihoodApproximations,
    llApproximationFunction = llApproximationFunction
  )
  null <- fit$par
  null[2] <- 1 / sqrt(null[2])
  names(null) <- c("mean", "sd")
  class(null) <- "null"
  return(null)
}
