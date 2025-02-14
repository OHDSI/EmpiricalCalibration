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

#' Calibrate the log likelihood ratio
#'
#' @description
#' \code{calibrateLlr} computes calibrated log likelihood ratio using the fitted null distribution
#'
#' @param likelihoodApproximation Either a data frame containing normal, skew-normal, or custom parametric likelihood
#'                                approximations, or a list of (adaptive) grid likelihood profiles.
#' @param null      An object of class \code{null} created using the \code{fitNull} function or an
#'                  object of class \code{mcmcNull} created using the \code{fitMcmcNull} function.
#' @param twoSided  Compute two-sided (TRUE) or one-sided (FALSE) p-value?
#' @param upper     If one-sided: compute p-value for upper (TRUE) or lower (FALSE) bound?
#'
#' @return
#' The calibrated log likelihood ratio.
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
  if (twoSided || !upper) {
    stop("Currently only one-sided upper LLRs are supported")
  }

  if (is.data.frame(likelihoodApproximation) || is.numeric(likelihoodApproximation)) {
    if ("logRr" %in% colnames(likelihoodApproximation)) {
      message("Detected data following normal distribution")
      type <- "normal"
      logLikelihood <- normalLlApproximaton
      likelihoodApproximation <- split(likelihoodApproximation, 1:nrow(likelihoodApproximation))
    } else if ("gamma" %in% colnames(likelihoodApproximation)) {
      message("Detected data following custom parameric distribution")
      type <- "custom"
      logLikelihood <- customLlApproximation
      likelihoodApproximation <- split(likelihoodApproximation, 1:nrow(likelihoodApproximation))
    } else if ("alpha" %in% colnames(likelihoodApproximation)) {
      message("Detected data following skew normal distribution")
      type <- "skewNormal"
      logLikelihood <- skewNormalLlApproximation
      likelihoodApproximation <- split(likelihoodApproximation, 1:nrow(likelihoodApproximation))
    } else if ("point" %in% colnames(likelihoodApproximation)) {
      message("Detected data following grid distribution")
      type <- "grid"
      logLikelihood <- gridLlApproximation
      likelihoodApproximation <- list(likelihoodApproximation)
    } else {
      message("Detected data following grid distribution")
      type <- "grid"
      logLikelihood <- gridLlApproximation
      if (is.numeric(likelihoodApproximation)) {
        likelihoodApproximation <- as.data.frame(t(likelihoodApproximation))
      }
      point <- as.numeric(colnames(likelihoodApproximation))
      if (any(is.na(point))) {
        stop("Expecting grid data, but not all column names are numeric")
      }
      convertToDataFrame <- function(i) {
        value <- as.numeric(likelihoodApproximation[i, ])
        maxValue <- max(value)
        data.frame(
          value = value - maxValue,
          point = point
        )
      }
      likelihoodApproximation <- lapply(1:nrow(likelihoodApproximation), convertToDataFrame)
    }
  } else {
    message("Detected data following grid distribution")
    type <- "grid"
    logLikelihood <- gridLlApproximation
  }
  useMcmc <- is(null, "mcmcNull")

  calibrateOneLlr <- function(i) {
    parameters <- likelihoodApproximation[[i]]
    if (type == "grid") {
      idx <- which(parameters$value == max(parameters$value))[1]
      mle <- parameters$point[idx]
      ml <- parameters$value[idx]
    } else if (type == "normal") {
      if (is.na(parameters$logRr) ||
        is.na(parameters$seLogRr) ||
        is.infinite(parameters$logRr) ||
        is.infinite(parameters$seLogRr)) {
        return(NA)
      }
      mle <- parameters$logRr
      ml <- logLikelihood(mle, parameters = parameters)
    } else {
      optimum <- suppressWarnings(optim(0, function(x) -logLikelihood(x = x, parameters = parameters)))
      mle <- optimum$par
      ml <- -optimum$value
    }
    if (useMcmc) {
      chain <- attr(null, "mcmc")$chain
      calibratedLlr <- mapply(
        FUN = calibrateOneLlrOneNull,
        nullMu = chain[, 1],
        nullSigma = 1 / sqrt(chain[, 2]),
        MoreArgs = list(
          logLikelihood = logLikelihood,
          parameters = parameters,
          mle = mle,
          ml = ml
        )
      )
      result <- quantile(calibratedLlr, c(0.5, 0.025, 0.975))
      return(result)
    } else {
      result <- calibrateOneLlrOneNull(
        nullMu = null[1],
        nullSigma = null[2],
        logLikelihood = logLikelihood,
        parameters = parameters,
        mle = mle,
        ml = ml
      )
      return(result)
    }
  }

  calibratedLlr <- sapply(1:length(likelihoodApproximation), calibrateOneLlr)
  if (useMcmc) {
    calibratedLlr <- as.data.frame(t(calibratedLlr))
    colnames(calibratedLlr) <- c("llr", "ci95Lb", "ci95Ub")
  } else {
    names(calibratedLlr) <- NULL
  }
  return(calibratedLlr)
}

computePFromLlr <- function(llr, mle, nullLogRr = 0) {
  # For very large llr values pchisq returns 0.
  # Behaves roughly log linear in that region, so using linear extrapolation.
  # llr <- 20:33
  # p <- (1 - pchisq(2 * llr, df = 1))/2
  # plot(llr, log(p))
  # fit <- lm(log(p) ~ llr)
  # # Recompute intercept at cut point for smooth integration:
  # formatC(log(p[length(p)]) - llr[length(llr)] * coef(fit)[2], digits = 20)
  # formatC(coef(fit)[2], digits = 20)
  p <- ifelse(llr > 33,
    exp(-2.4039960592000753081 - 1.0193835554520327413 * llr),
    (1 - pchisq(2 * llr, df = 1)) / 2
  )
  ifelse(mle < nullLogRr, 1 - p, p)
}

computeLlrFromP <- function(p) {
  if (p >= 0.5) {
    return(0)
  } else {
    if (p < 1e-16) {
      # For very small p values qchisq returns Inf.
      # Behaves roughly log linear in that region, so using linear extrapolation.
      # p <- exp(-10:-25)
      # llr <- qchisq(1 - 2 * p, df = 1) / 2
      # plot(log(p), llr)
      # fit <- lm(llr ~ log(p))
      # coef(fit)
      return(-2.0604312 - 0.9677458 * log(p))
    } else {
      return(qchisq(1 - 2 * p, df = 1) / 2)
    }
  }
}

computePAtNull <- function(x, logLikelihood, parameters, mle, ml, null) {
  llr <- ml - logLikelihood(x = x, parameters = parameters)
  p <- computePFromLlr(llr, mle = mle, nullLogRr = x)
  p * dnorm(x, null[1], null[2])
}

calibrateOneLlrOneNull <- function(nullMu, nullSigma, logLikelihood, parameters, mle, ml) {
  if (nullSigma < 0.001) {
    if (mle < nullMu) {
      calibratedLlr <- 0
    } else {
      calibratedLlr <- ml - logLikelihood(x = nullMu, parameters = parameters)
    }
  } else {
    calibratedP <- integrate(
      f = computePAtNull,
      logLikelihood = logLikelihood,
      parameters = parameters,
      mle = mle,
      ml = ml,
      null = c(nullMu, nullSigma),
      lower = nullMu - 10 * nullSigma,
      upper = nullMu + 10 * nullSigma,
      stop.on.error = FALSE
    )$value
    calibratedLlr <- computeLlrFromP(calibratedP)
  }
  return(calibratedLlr)
}
