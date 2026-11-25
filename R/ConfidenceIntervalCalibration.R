# @file ConfidenceIntervalCalibration.R
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

#' Fit a systematic error model
#'
#' @details
#' Fit a model of the systematic error as a function of true effect size. This model is an extension
#' of the method for fitting the null distribution. The mean and log(standard deviations) of the error
#' distributions are assumed to be linear with respect to the true effect size, and each component is
#' therefore represented by an intercept and a slope.
#'
#' @param estimate                   A numeric vector of effect estimates. On the log relative risk scale
#'                                   these are log relative risks.
#' @param se                         The standard error of the effect estimates. Hint: on the log relative
#'                                   risk scale often the standard error = (log(<lower bound 95 percent
#'                                   confidence interval>) - log(<effect estimate>))/qnorm(0.025).
#' @param trueEffect                 A vector of the true effect sizes (on the same scale as `estimate`).
#' @param estimateCovarianceMatrix   Should a covariance matrix be computed? If so, confidence
#'                                   intervals for the model parameters will be available.
#' @param legacy                     If true, a legacy error model will be fitted, meaning standard
#'                                   deviation is linear on the log scale. If false, standard deviation
#'                                   is assumed to be simply linear. The legacy model is not supported for the
#'                                   risk difference scale.
#' @param logRr                      DEPRECATED. Use `estimate` instead.
#' @param seLogRr                    DEPRECATED. Use `se` instead.
#' @param trueLogRr                  DEPRECATED. Use `trueEffect` instead.
#'
#' @return
#' An object of type `systematicErrorModel`.
#'
#' @examples
#' controls <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' model <- fitSystematicErrorModel(controls$estimate, controls$se, controls$trueEffect)
#' model
#'
#' @export
fitSystematicErrorModel <- function(estimate,
                                    se,
                                    trueEffect,
                                    estimateCovarianceMatrix = FALSE,
                                    legacy = FALSE,
                                    logRr = NULL,
                                    seLogRr = NULL,
                                    trueLogRr = NULL) {
  if (!is.null(logRr)) {
    warnDeprecatedArguments(c("logRr", "seLogRr", "trueLogRr"), c("estimate", "se", "trueEffect"))
    estimate <- logRr
    se <- seLogRr
    trueEffect <- trueLogRr
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimate, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(trueEffect, min.len = 1, add = errorMessages)
  checkmate::assertLogical(estimateCovarianceMatrix, len = 1, add = errorMessages)
  checkmate::assertLogical(legacy, len = 1, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(trueLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (any(is.infinite(se))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting error model")
    trueEffect <- trueEffect[!is.infinite(se)]
    estimate <- estimate[!is.infinite(se)]
    se <- se[!is.infinite(se)]
  }
  if (any(is.infinite(estimate))) {
    warning("Estimate(s) with infinite estimate detected. Removing before fitting error model")
    trueEffect <- trueEffect[!is.infinite(estimate)]
    se <- se[!is.infinite(estimate)]
    estimate <- estimate[!is.infinite(estimate)]
  }
  if (any(is.na(se))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting error model")
    trueEffect <- trueEffect[!is.na(se)]
    estimate <- estimate[!is.na(se)]
    se <- se[!is.na(se)]
  }
  if (any(is.na(estimate))) {
    warning("Estimate(s) with NA estimate detected. Removing before fitting error model")
    trueEffect <- trueEffect[!is.na(estimate)]
    se <- se[!is.na(estimate)]
    estimate <- estimate[!is.na(estimate)]
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
  
  # Mean intercept and sd intercept are in units of the estimates:
  fit <- optim(theta,
               logLikelihood,
               logRr = estimate,
               seLogRr = se,
               trueLogRr = trueEffect,
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
    fisher_info <- fisher_info 
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
#' @param estimate  A numeric vector of effect estimates. On the log relative risk scale these are
#'                  log relative risks.
#' @param se        The standard error of the effect estimates. Hint: on the log relative risk scale
#'                  often the standard error = (log(<lower bound 95 percent confidence interval>) -
#'                  log(<effect estimate>))/qnorm(0.025).
#' @param ciWidth   The width of the confidence interval. Typically this would be .95, for the 95
#'                  percent confidence interval.
#' @param model     An object of type `systematicErrorModel` as created by the
#'                  [fitSystematicErrorModel()] function.
#' @template Scale
#' @param logRr     DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr   DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#'
#' @return
#' A data frame with calibrated confidence intervals and point estimates, on the same scale as the
#' provided estimates (so for example on the log scale for log relative risks). The data frame has columns
#' `estimate`, `lb`, `ub`, and `se`. For backwards compatibility it also contains the deprecated columns
#' `logRr`, `logLb95Rr`, `logUb95Rr`, and `seLogRr` with identical values.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' model <- fitSystematicErrorModel(data$estimate, data$se, data$trueEffect)
#' newData <- simulateControls(n = 15, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' result <- calibrateConfidenceInterval(newData$estimate, newData$se, model)
#' result
#'
#' @export
calibrateConfidenceInterval <- function(estimate,
                                        se,
                                        model,
                                        ciWidth = 0.95,
                                        scale = "log relative risk",
                                        logRr = NULL,
                                        seLogRr = NULL) {
  if (!is.null(logRr)) {
    warnDeprecatedArguments(c("logRr", "seLogRr"), c("estimate", "se"))
    estimate <- logRr
    se <- seLogRr
    scale <- LOG_RELATIVE_RISK
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimate, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se, min.len = 1, add = errorMessages)
  checkmate::assertClass(model, classes = "systematicErrorModel", add = errorMessages)
  checkmate::assertNumber(ciWidth, add = errorMessages)
  checkmate::assertChoice(scale, choices = c(LOG_RELATIVE_RISK, RISK_DIFFERENCE), add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

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
    # Root finding tolerance is in units of the estimates, so needs to follow the magnitude of the
    # standard error for absolute scales:
    tol <- ifelse(scale == RISK_DIFFERENCE, 1e-4 * max(se, 1e-10), .Machine$double.eps^0.25)
    return(uniroot(
      f = opt,
      interval = c(lower, upper),
      tol = tol,
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
  
  if (names(model)[3] == "logSdIntercept" && scale == RISK_DIFFERENCE) {
    stop("The legacy error model is only supported for the log relative risk scale")
  }
  legacy <- (names(model)[3] == "logSdIntercept")
  result <- data.frame(estimate = rep(as.numeric(NA), length(estimate)), lb = as.numeric(NA), ub = as.numeric(NA))
  if (!is.na(model[1])) {
    for (i in 1:nrow(result)) {
      if (is.infinite(estimate[i]) || is.na(estimate[i]) || is.infinite(se[i]) || is.na(se[i])) {
        result$estimate[i] <- NA
        result$lb[i] <- NA
        result$ub[i] <- NA
      } else {
        result$estimate[i] <- logBound(
          0,
          TRUE,
          estimate[i],
          se[i],
          model[1],
          model[2],
          model[3],
          model[4],
          legacy
        )
        result$lb[i] <- logBound(
          ciWidth,
          TRUE,
          estimate[i],
          se[i],
          model[1],
          model[2],
          model[3],
          model[4],
          legacy
        )
        result$ub[i] <- logBound(
          ciWidth,
          FALSE,
          estimate[i],
          se[i],
          model[1],
          model[2],
          model[3],
          model[4],
          legacy
        )
      }
    }
  }
  result$se <- (result$lb - result$ub) / (2 * qnorm((1 - ciWidth) / 2))
  result$logRr <- result$estimate
  result$logLb95Rr <- result$lb
  result$logUb95Rr <- result$ub
  result$seLogRr <- result$se
  return(result)
}

#' Compute the (traditional) confidence interval
#'
#' @description
#' `computeTraditionalCi` computes the traditional confidence interval based on the effect estimate
#' (for example the log of the relative risk) and its standard error.
#'
#' @param estimate  A numeric vector of one or more effect estimates. On the log relative risk scale these
#'                  are log relative risks.
#' @param se        The standard error of the effect estimates. Hint: on the log relative risk scale
#'                  often the standard error = (log(<lower bound 95 percent confidence interval>) -
#'                  log(<effect estimate>))/qnorm(0.025)
#' @param ciWidth   The width of the confidence interval. Typically this would be .95, for the 95
#'                  percent confidence interval.
#' @template Scale
#' @param logRr     DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr   DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#'
#' @return
#' The point estimate and confidence interval, in a data frame with columns `rr`, `lb`, and `ub` (all
#' exponentiated, so on the relative risk scale) when `scale = "log relative risk"`, or columns `estimate`, `lb`,
#' and `ub` when `scale = "risk difference"`.
#'
#' @examples
#' data(sccs)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' computeTraditionalCi(positive$logRr, positive$seLogRr)
#'
#' @export
computeTraditionalCi <- function(estimate,
                                 se,
                                 ciWidth = .95,
                                 scale = "log relative risk",
                                 logRr = NULL,
                                 seLogRr = NULL) {
  if (!is.null(logRr)) {
    warnDeprecatedArguments(c("logRr", "seLogRr"), c("estimate", "se"))
    estimate <- logRr
    se <- seLogRr
    scale <- LOG_RELATIVE_RISK
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimate, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se, min.len = 1, add = errorMessages)
  checkmate::assertNumber(ciWidth, add = errorMessages)
  checkmate::assertChoice(scale, choices = c(LOG_RELATIVE_RISK, RISK_DIFFERENCE), add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  lb <- estimate + qnorm((1 - ciWidth) / 2) * se
  ub <- estimate - qnorm((1 - ciWidth) / 2) * se
  if (scale == LOG_RELATIVE_RISK) {
    return(data.frame(rr = exp(estimate), lb = exp(lb), ub = exp(ub)))
  } else {
    return(data.frame(estimate = estimate, lb = lb, ub = ub))
  }
}

#' Convert empirical null distribution to systematic error model
#'
#' @description
#' This function converts an empirical null distribution, fitted using estimates only for negative controls,
#' into a systematic error distribution that can be used to calibrate confidence intervals in addition to
#' p-values.
#'
#' Whereas the [fitSystematicErrorModel()] uses positive controls to determine how the error
#' distribution changes with true effect size, this function requires the user to make an assumption. The
#' default assumption, `meanSlope = 1` and `sdSlope = 0`, specify a belief that the error
#' distribution is the same for all true effect sizes. In many cases this assumption is likely to be correct,
#' however, if an estimation method is biased towards the null this assumption will be violated, causing the
#' calibrated confidence intervals to have lower than nominal coverage.
#'
#' @param null         The empirical null distribution fitted using either the [fitNull()]
#'                     or the [fitMcmcNull()] function.
#' @param meanSlope    The slope for the mean of the error distribution. A slope of 1 means the error
#'                     is the same for different values of the true relative risk.
#' @param sdSlope   The slope for the log of the standard deviation of the error distribution. A slope
#'                     of 0 means the standard deviation is the same for different values of the true
#'                     relative risk.
#'
#' @return An object of type `systematicErrorModel`.
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
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumber(meanSlope, add = errorMessages)
  checkmate::assertNumber(sdSlope, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is(null, "mcmcNull")) {
    model <- c(null[1], meanSlope, 1 / sqrt(null[2]), sdSlope)
  } else if (is(null, "null")) {
    model <- c(null[1], meanSlope, null[2], sdSlope)
  } else {
    stop("Null argument should be of type 'null' or 'mcmcNull'")
  }
  names(model) <- c("meanIntercept", "meanSlope", "sdIntercept", "sdSlope")
  class(model) <- "systematicErrorModel"
  model
}
