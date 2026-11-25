# @file EmpiricalCalibrationUsingAsymptotics.R
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

#' Fit the null distribution
#'
#' @description
#' `fitNull` fits the null distribution to a set of negative controls
#'
#' @details
#' This function fits a Gaussian function to the negative control estimates as described in Schuemie
#' et al (2014).
#'
#'
#' @param estimate  A numeric vector of effect estimates. On the log relative risk scale these are
#'                  log relative risks.
#' @param se        The standard error of the effect estimates. Hint: on the log relative risk scale
#'                  often the standard error = (log(<lower bound 95 percent confidence interval>) -
#'                  log(<effect estimate>))/qnorm(0.025)
#' @param logRr     DEPRECATED. Use `estimate` instead.
#' @param seLogRr   DEPRECATED. Use `se` instead.
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
fitNull <- function(estimate, se, logRr = NULL, seLogRr = NULL) {
  if (!is.null(logRr)) {
    warnDeprecatedArguments(c("logRr", "seLogRr"), c("estimate", "se"))
    estimate <- logRr
    se <- seLogRr
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimate, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (any(is.infinite(se))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting null distribution")
    estimate <- estimate[!is.infinite(se)]
    se <- se[!is.infinite(se)]
  }
  if (any(is.infinite(estimate))) {
    warning("Estimate(s) with infinite estimate detected. Removing before fitting null distribution")
    se <- se[!is.infinite(estimate)]
    estimate <- estimate[!is.infinite(estimate)]
  }
  if (any(is.na(se))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting null distribution")
    estimate <- estimate[!is.na(se)]
    se <- se[!is.na(se)]
  }
  if (any(is.na(estimate))) {
    warning("Estimate(s) with NA estimate detected. Removing before fitting null distribution")
    se <- se[!is.na(estimate)]
    estimate <- estimate[!is.na(estimate)]
  }
  if (any(abs(estimate) > log(100))) {
    warning("Extreme estimate(s) detected: abs(estimate) > log(100). Removing before fitting null distribution")
    se <- se[abs(estimate) <= log(100)]
    estimate <- estimate[abs(estimate) <= log(100)]
  }
  if (length(estimate) == 0) {
    warning("No estimates remaining")
    null <- c(NA, NA)
  } else {
    theta <- c(0, 1)
    fit <- optim(theta, logLikelihoodNull, logRr = estimate, seLogRr = se)
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
#' `calibrateP` computes calibrated p-values using the fitted null distribution
#'
#' @details
#' This function computes a calibrated two-sided p-value as described in Schuemie et al (2014).
#'
#' @param estimate  A numeric vector of one or more effect estimates, on the same scale as used when
#'                  fitting the null (for example log relative risks).
#' @param se        The standard error of the effect estimates. Hint: on the log relative risk scale
#'                  often the standard error = (log(<lower bound 95 percent confidence interval>) -
#'                  log(<effect estimate>))/qnorm(0.025)
#' @param null      An object of class `null` created using the `fitNull` function or an
#'                  object of class `mcmcNull` created using the `fitMcmcNull` function.
#' @param twoSided  Compute two-sided (TRUE) or one-sided (FALSE) p-value?
#' @param upper     If one-sided: compute p-value for upper (TRUE) or lower (FALSE) bound?
#' @param ...       Any additional parameters (currently none).
#' @param logRr     DEPRECATED. Use `estimate` instead.
#' @param seLogRr   DEPRECATED. Use `se` instead.
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
calibrateP <- function(null, estimate, se, twoSided = TRUE, upper = TRUE, ..., logRr = NULL, seLogRr = NULL) {
  if (!is.null(logRr)) {
    warnDeprecatedArguments(c("logRr", "seLogRr"), c("estimate", "se"))
    estimate <- logRr
    se <- seLogRr
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(null, classes = "null", add = errorMessages)
  checkmate::assertNumeric(estimate, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se, min.len = 1, add = errorMessages)
  checkmate::assertLogical(twoSided, len = 1, add = errorMessages)
  checkmate::assertLogical(upper, len = 1, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  UseMethod("calibrateP")
}


#' @describeIn
#' calibrateP Computes the calibrated P-value using asymptotic assumptions.
#'
#' @export
calibrateP.null <- function(null, estimate, se, twoSided = TRUE, upper = TRUE, ..., logRr = NULL, seLogRr = NULL) {
  if (!is.null(logRr)) {
    warnDeprecatedArguments(c("logRr", "seLogRr"), c("estimate", "se"))
    estimate <- logRr
    se <- seLogRr
  }
  if (length(estimate) != length(se)) {
    stop("The estimate and se arguments must be of equal length")
  }
  
  calibrateOneP <- function(i) {
    if (is.na(estimate[i]) || is.infinite(estimate[i]) || is.na(se[i]) || is.infinite(se[i])) {
      return(NA)
    } else {
      pUpperBound <- pnorm((null[1] - estimate[i]) / sqrt(null[2]^2 + se[i]^2))
    }
    pLowerBound <- pnorm((estimate[i] - null[1]) / sqrt(null[2]^2 + se[i]^2))
    if (twoSided) {
      return(2 * min(pUpperBound, pLowerBound))
    } else if (upper) {
      return(pUpperBound)
    } else {
      return(pLowerBound)
    }
  }
  
  calibratedP <- sapply(1:length(estimate), calibrateOneP)
  names(calibratedP) <- NULL
  return(calibratedP)
}

#' Compute the (traditional) p-value
#'
#' @description
#' `computeTraditionalP` computes the traditional two-sided p-value based on the effect estimate
#' (for example the log of the relative risk) and its standard error.
#'
#' @param estimate  A numeric vector of one or more effect estimates (for example log relative risks).
#' @param se        The standard error of the effect estimates. Hint: on the log relative risk scale
#'                  often the standard error = (log(<lower bound 95 percent confidence interval>) -
#'                  log(<effect estimate>))/qnorm(0.025)
#' @param twoSided  Compute two-sided (TRUE) or one-sided (FALSE) p-value?
#' @param upper     If one-sided: compute p-value for upper (TRUE) or lower (FALSE) bound?
#' @param logRr     DEPRECATED. Use `estimate` instead.
#' @param seLogRr   DEPRECATED. Use `se` instead.
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
computeTraditionalP <- function(estimate, se, twoSided = TRUE, upper = TRUE, logRr = NULL, seLogRr = NULL) {
  if (!is.null(logRr)) {
    warnDeprecatedArguments(c("logRr", "seLogRr"), c("estimate", "se"))
    estimate <- logRr
    se <- seLogRr
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimate, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se, min.len = 1, add = errorMessages)
  checkmate::assertLogical(twoSided, len = 1, add = errorMessages)
  checkmate::assertLogical(upper, len = 1, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  z <- estimate / se
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
#' `fitNullNonNormalLl` fits the null distribution to a set of negative controls
#'
#' @details
#' This function fits a Gaussian function to the negative control estimates, using non-normal
#' approximations of the per-negative control log likelihood.
#'
#' Note: for non-normal approximations the optimization starts at a null with mean 0 and standard deviation 1, and
#' is not rescaled, so these may not be reliable for estimates on a scale with very small standard errors, such as
#' risk differences.
#'
#' @param likelihoodApproximations Either a data frame containing normal, skew-normal, or custom parametric likelihood
#'                                 approximations, or a list of (adaptive) grid likelihood profiles. Normal
#'                                 approximations are identified by the `estimate` and `se` columns (or the
#'                                 deprecated `logRr` and `seLogRr` columns).
#' @param scale                    The scale of the effect estimates. Can be `"log relative risk"` or
#'                                 `"risk difference"`. Currently only used when the approximations are
#'                                 normal (see details).
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
fitNullNonNormalLl <- function(likelihoodApproximations, scale = "log relative risk") {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assert(is.data.frame(likelihoodApproximations) ||
                     is.matrix(likelihoodApproximations) ||
                     is.list(likelihoodApproximations), add = errorMessages)
  checkmate::assertString(scale, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.data.frame(likelihoodApproximations)) {
    likelihoodApproximations <- renameLegacyNormalColumns(likelihoodApproximations)
    if ("estimate" %in% colnames(likelihoodApproximations)) {
      message("Detected data following normal distribution")
      return(fitNull(
        estimate = likelihoodApproximations$estimate,
        se = likelihoodApproximations$se
      ))
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
  
  theta <- c(0, 1)
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
