# @file ConfidenceIntervalCalibration.R
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

#' Fit a systematic error model
#'
#' @details
#' Fit a model of the systematic error as a function of true effect size. This model is an extension
#' of the method for fitting the null distribution. The mean and log(standard deviations) of the error
#' distributions are assumed to be linear with respect to the true effect size, and each component is
#' therefore represented by an intercept and a slope.
#'
#' @param logRr                      A numeric vector of effect estimates on the log scale.
#' @param seLogRr                    The standard error of the log of the effect estimates. Hint: often
#'                                   the standard error = (log(<lower bound 95 percent confidence
#'                                   interval>) - log(<effect estimate>))/qnorm(0.025).
#' @param trueLogRr                  A vector of the true effect sizes.
#' @param estimateCovarianceMatrix   Should a covariance matrix be computed? If so, confidence
#'                                   intervals for the model parameters will be available.
#' @param legacy                     If true, a legacy error model will be fitted, meaning standard
#'                                   deviation is linear on the log scale. If false, standard deviation
#'                                   is assumed to be simply linear.
#'
#' @return
#' An object of type \code{systematicErrorModel}.
#'
#' @examples
#' controls <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' model <- fitSystematicErrorModel(controls$logRr, controls$seLogRr, controls$trueLogRr)
#' model
#'
#' @export
fitSystematicErrorModel <- function(logRr,
                                    seLogRr,
                                    trueLogRr,
                                    estimateCovarianceMatrix = FALSE,
                                    legacy = FALSE) {
  if (any(is.infinite(seLogRr))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.infinite(seLogRr)]
    logRr <- logRr[!is.infinite(seLogRr)]
    seLogRr <- seLogRr[!is.infinite(seLogRr)]
  }
  if (any(is.infinite(logRr))) {
    warning("Estimate(s) with infinite logRr detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.infinite(logRr)]
    seLogRr <- seLogRr[!is.infinite(logRr)]
    logRr <- logRr[!is.infinite(logRr)]
  }
  if (any(is.na(seLogRr))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.na(seLogRr)]
    logRr <- logRr[!is.na(seLogRr)]
    seLogRr <- seLogRr[!is.na(seLogRr)]
  }
  if (any(is.na(logRr))) {
    warning("Estimate(s) with NA logRr detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.na(logRr)]
    seLogRr <- seLogRr[!is.na(logRr)]
    logRr <- logRr[!is.na(logRr)]
  }
  
  if (legacy) {
    theta <- c(0, 1, -2, 0)
    logLikelihood <- minLogLikelihoodErrorModelLegacy
    parscale <- c(1, 1, 10, 10)
  } else {
    theta <- c(0, 1, 0.1, 0)
    logLikelihood <- minLogLikelihoodErrorModel
    parscale <- c(1, 1, 1, 1)
  }
  
  fit <- optim(theta,
               logLikelihood,
               logRr = logRr,
               seLogRr = seLogRr,
               trueLogRr = trueLogRr,
               method = "BFGS",
               hessian = TRUE,
               control = list(parscale = parscale)
  )
  model <- fit$par
  if (legacy) {
    names(model) <- c("meanIntercept", "meanSlope", "logSdIntercept", "logSdSlope")
  } else {
    names(model) <- c("meanIntercept", "meanSlope", "sdIntercept", "sdSlope")
  }
  if (estimateCovarianceMatrix) {
    fisher_info <- solve(fit$hessian)
    prop_sigma <- sqrt(diag(fisher_info))
    attr(model, "CovarianceMatrix") <- fisher_info
    attr(model, "LB95CI") <- fit$par + qnorm(0.025) * prop_sigma
    attr(model, "UB95CI") <- fit$par + qnorm(0.975) * prop_sigma
  }
  class(model) <- "systematicErrorModel"
  model
}

#' Calibrate confidence intervals
#'
#' @details
#' Compute calibrated confidence intervals based on a model of the systematic error.
#'
#' @param logRr     A numeric vector of effect estimates on the log scale.
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025).
#' @param ciWidth   The width of the confidence interval. Typically this would be .95, for the 95
#'                  percent confidence interval.
#' @param model     An object of type \code{systematicErrorModel} as created by the
#'                  \code{\link{fitSystematicErrorModel}} function.
#'
#' @return
#' A data frame with calibrated confidence intervals and point estimates.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' model <- fitSystematicErrorModel(data$logRr, data$seLogRr, data$trueLogRr)
#' newData <- simulateControls(n = 15, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' result <- calibrateConfidenceInterval(newData$logRr, newData$seLogRr, model)
#' result
#'
#' @export
calibrateConfidenceInterval <- function(logRr, seLogRr, model, ciWidth = 0.95) {
  opt <- function(x,
                  z,
                  logRr,
                  se,
                  interceptMean,
                  slopeMean,
                  interceptSd,
                  slopeSd,
                  legacy) {
    mean <- interceptMean + slopeMean * x
    if (legacy) {
      sd <- exp(interceptSd + slopeSd * x)
    } else {
      sd <- interceptSd + slopeSd * abs(x)
    }
    numerator <- mean - logRr
    denominator <- sqrt((sd)^2 + (se)^2)
    if (is.infinite(denominator)) {
      if (numerator > 0) {
        return(z + 1e-10)
      } else {
        return(z - 1e-10)
      }
    }
    return(z + numerator / denominator)
  }
  
  logBound <- function(ciWidth,
                       lb = TRUE,
                       logRr,
                       se,
                       interceptMean,
                       slopeMean,
                       interceptSd,
                       slopeSd,
                       legacy) {
    z <- qnorm((1 - ciWidth) / 2)
    if (lb) {
      z <- -z
    }
    if (legacy) {
      # Simple grid search for upper bound where opt is still positive:
      if (slopeSd > 0) {
        lower <- -100
        upper <- lower
        while (opt(
          x = upper,
          z = z,
          logRr = logRr,
          se = se,
          interceptMean = interceptMean,
          slopeMean = slopeMean,
          interceptSd = interceptSd,
          slopeSd = slopeSd,
          legacy = legacy
        ) < 0 && upper < 100) {
          upper <- upper + 1
        }
        if (upper == lower | upper == 100) {
          return(NA)
        }
      } else {
        upper <- 100
        lower <- upper
        while (opt(
          x = lower,
          z = z,
          logRr = logRr,
          se = se,
          interceptMean = interceptMean,
          slopeMean = slopeMean,
          interceptSd = interceptSd,
          slopeSd = slopeSd,
          legacy = legacy
        ) > 0 && lower > -100) {
          lower <- lower - 1
        }
        if (upper == lower | lower == -100) {
          return(NA)
        }
      }
    } else {
      lower <- -100
      upper <- 100
      lowerValue <- opt(
        x = lower,
        z = z,
        logRr = logRr,
        se = se,
        interceptMean = interceptMean,
        slopeMean = slopeMean,
        interceptSd = interceptSd,
        slopeSd = slopeSd,
        legacy = legacy
      )
      upperValue <- opt(
        x = upper,
        z = z,
        logRr = logRr,
        se = se,
        interceptMean = interceptMean,
        slopeMean = slopeMean,
        interceptSd = interceptSd,
        slopeSd = slopeSd,
        legacy = legacy
      )
      if ((lowerValue < 0 && upperValue < 0) || (lowerValue > 0 && upperValue > 0)) {
        return(NA)
      }
    }
    return(uniroot(
      f = opt,
      interval = c(lower, upper),
      z = z,
      logRr = logRr,
      se = se,
      interceptMean = interceptMean,
      slopeMean = slopeMean,
      interceptSd = interceptSd,
      slopeSd = slopeSd,
      legacy = legacy
    )$root)
  }
  
  legacy <- (names(model)[3] == "logSdIntercept")
  result <- data.frame(logRr = rep(as.numeric(NA), length(logRr)), logLb95Rr = as.numeric(NA), logUb95Rr = as.numeric(NA))
  if (!is.na(model[1])) {
    for (i in 1:nrow(result)) {
      if (is.infinite(logRr[i]) || is.na(logRr[i]) || is.infinite(seLogRr[i]) || is.na(seLogRr[i])) {
        result$logRr[i] <- NA
        result$logLb95Rr[i] <- NA
        result$logUb95Rr[i] <- NA
      } else {
        result$logRr[i] <- logBound(
          0,
          TRUE,
          logRr[i],
          seLogRr[i],
          model[1],
          model[2],
          model[3],
          model[4],
          legacy
        )
        result$logLb95Rr[i] <- logBound(
          ciWidth,
          TRUE,
          logRr[i],
          seLogRr[i],
          model[1],
          model[2],
          model[3],
          model[4],
          legacy
        )
        result$logUb95Rr[i] <- logBound(
          ciWidth,
          FALSE,
          logRr[i],
          seLogRr[i],
          model[1],
          model[2],
          model[3],
          model[4],
          legacy
        )
      }
    }
  }
  result$seLogRr <- (result$logLb95Rr - result$logUb95Rr) / (2 * qnorm((1 - ciWidth) / 2))
  return(result)
}

#' Compute the (traditional) confidence interval
#'
#' @description
#' \code{computeTraditionalCi} computes the traditional confidence interval based on the log of the
#' relative risk and the standard error of the log of the relative risk.
#'
#' @param logRr     A numeric vector of one or more effect estimates on the log scale
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#' @param ciWidth   The width of the confidence interval. Typically this would be .95, for the 95
#'                  percent confidence interval.
#'
#' @return
#' The point estimate and confidence interval
#'
#' @examples
#' data(sccs)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' computeTraditionalCi(positive$logRr, positive$seLogRr)
#'
#' @export
computeTraditionalCi <- function(logRr, seLogRr, ciWidth = .95) {
  return(data.frame(
    rr = exp(logRr),
    lb = exp(logRr + qnorm((1 - ciWidth) / 2) * seLogRr),
    ub = exp(logRr - qnorm((1 - ciWidth) / 2) * seLogRr)
  ))
}

#' Convert empirical null distribution to systematic error model
#'
#' @description
#' This function converts an empirical null distribution, fitted using estimates only for negative controls,
#' into a systematic error distribution that can be used to calibrate confidence intervals in addition to
#' p-values.
#'
#' Whereas the \code{\link{fitSystematicErrorModel}} uses positive controls to determine how the error
#' distribution changes with true effect size, this function requires the user to make an assumption. The
#' default assumption, \code{meanSlope = 1} and \code{sdSlope = 0}, specify a belief that the error
#' distribution is the same for all true effect sizes. In many cases this assumption is likely to be correct,
#' however, if an estimation method is biased towards the null this assumption will be violated, causing the
#' calibrated confidence intervals to have lower than nominal coverage.
#'
#' @param null         The empirical null distribution fitted using either the \code{\link{fitNull}}
#'                     or the \code{\link{fitMcmcNull}} function.
#' @param meanSlope    The slope for the mean of the error distribution. A slope of 1 means the error
#'                     is the same for different values of the true relative risk.
#' @param sdSlope   The slope for the log of the standard deviation of the error distribution. A slope
#'                     of 0 means the standard deviation is the same for different values of the true
#'                     relative risk.
#'
#' @return An object of type \code{systematicErrorModel}.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' model <- convertNullToErrorModel(null)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' calibrateConfidenceInterval(positive$logRr, positive$seLogRr, model)
#'
#' @export
convertNullToErrorModel <- function(null, meanSlope = 1, sdSlope = 0) {
  if (is(null, "null")) {
    model <- c(null[1], meanSlope, null[2], sdSlope)
  } else if (is(null, "mcmcNull")) {
    model <- c(null[1], meanSlope, 1 / sqrt(null[2]), sdSlope)
  } else {
    stop("Null argument should be of type 'null' or 'mcmcNull'")
  }
  names(model) <- c("meanIntercept", "meanSlope", "sdIntercept", "sdSlope")
  class(model) <- "systematicErrorModel"
  model
}
