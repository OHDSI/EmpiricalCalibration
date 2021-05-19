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
#' \code{calibrateLlr} computes calibrated log likelihood ratio using the fitted null distribution
#'
#' @param likelihoodApproximation    A data frame containing either normal, skew-normal, custom parametric, or grid
#'                  likelihood data. 
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
  if (twoSided || !upper)
    stop("Currently only one-sided upper LLRs are supported")
  
  if ("logRr" %in% colnames(likelihoodApproximation)) {
    message("Detected data following normal distribution")
    type <- "normal"
    logLikelihood <- function(x, row, ...) {
      dnorm(x, mean = row$logRr, sd = row$seLogRr, log = TRUE)
    }
    gridX <- NULL
  } else if ("gamma" %in% colnames(likelihoodApproximation)) {
    message("Detected data following custom parameric distribution")
    type <- "custom"
    logLikelihood <- function(x, row, ...) {
      ((exp(row$gamma * (x - row$mu)))) * ((-(x - row$mu)^2)/(2 * row$sigma^2))
    }
    gridX <- NULL
  } else if ("alpha" %in% colnames(likelihoodApproximation)) {
    message("Detected data following skew normal distribution")
    type <- "skewNormal"
    logLikelihood <- function(x, row, ...) {
      if (is.infinite(row$sigma)) {
        return(rep(0, length(x)))
      }
      return(log(2) + dnorm(x, row$mu, row$sigma, log = TRUE) + pnorm(row$alpha * (x - row$mu), 0, row$sigma, log.p = TRUE))
    }
    gridX <- NULL
  } else {
    message("Detected data following grid distribution")
    type <- "grid"
    logLikelihood <- gridLikelihood
    if (!is.data.frame(likelihoodApproximation)) {
      likelihoodApproximation <- as.data.frame(t(likelihoodApproximation))
    }
    gridX <- as.numeric(colnames(likelihoodApproximation))
    if (any(is.na(gridX))) {
      stop("Expecting grid data, but not all column names are numeric")
    }
  }
  useMcmc <- is(null, "mcmcNull")
  
  calibrateOneLlr <- function(i) {
    row <- likelihoodApproximation[i, ]
    if (type == "grid") {
      row <- as.numeric(row)
      idx <- which(row == max(row))[1]
      mle <- gridX[idx]
      ml <- row[idx]
    } else if (type == "normal") {
      mle <- row$logRr
      ml <- logLikelihood(row$logRr, row = row, gridX = gridX)
    } else {
      optimum <-  suppressWarnings(optim(0, function(x) -logLikelihood(x = x, row = row, gridX = gridX)))
      mle <- optimum$par
      ml <- -optimum$value
    }
    if (useMcmc) {
      chain <- attr(null, "mcmc")$chain
      calibratedLlr <- mapply(FUN = calibrateOneLlrOneNull, 
                              nullMu = chain[, 1], 
                              nullSigma = 1/sqrt(chain[, 2]),
                              MoreArgs = list(
                                logLikelihood = logLikelihood,
                                row = row,
                                gridX = gridX,
                                mle = mle,
                                ml = ml))
      result <- quantile(calibratedLlr, c(0.5, 0.025, 0.975))
      return(result)
    } else {
      result <- calibrateOneLlrOneNull(nullMu = null[1],
                                       nullSigma = null[2],
                                       logLikelihood = logLikelihood,
                                       row = row,
                                       gridX = gridX,
                                       mle = mle,
                                       ml = ml)
      return(result    )
    }
  }
  
  calibratedLlr <- sapply(1:nrow(likelihoodApproximation), calibrateOneLlr)
  if (useMcmc) {
    calibratedLlr <- as.data.frame(t(calibratedLlr))
    colnames(calibratedLlr) <- c("llr", "ci95Lb", "ci95Ub")
  } else {
    names(calibratedLlr) <- NULL
  }
  return(calibratedLlr)
}

computePFromLlr <- function(llr, mle, nullLogRr = 0) {
  ifelse(mle < nullLogRr, 1 - (1 - pchisq(2 * llr, df = 1))/2, (1 - pchisq(2 * llr, df = 1))/2)
}

computeLlrFromP <- function(p) {
  if (p >= 0.5) {
    return(0)
  } else {
    return(qchisq(1 - 2 * p, df = 1) / 2)
  }
}

computePAtNull <- function(x, logLikelihood, row, mle, ml, null, gridX) {
  llr <- ml - logLikelihood(x = x, row = row, gridX = gridX)
  p <- computePFromLlr(llr, mle = mle, nullLogRr = x)
  p * dnorm(x, null[1], null[2])
}

calibrateOneLlrOneNull <- function(nullMu, nullSigma, logLikelihood, row, gridX, mle, ml) {
  if (nullSigma < 0.001) {
    if (mle < nullMu) {
      calibratedLlr <- 0
    } else {
      calibratedLlr <- ml - logLikelihood(x = nullMu, row = row, gridX = gridX)
    }
  } else {
    calibratedP <- integrate(f = computePAtNull, 
                             logLikelihood = logLikelihood, 
                             row = row,
                             mle = mle,
                             ml = ml,
                             gridX = gridX,
                             null = c(nullMu, nullSigma),
                             lower = nullMu - 10 * nullSigma, 
                             upper = nullMu + 10 * nullSigma)$value
    calibratedLlr <- computeLlrFromP(calibratedP)
  }
  return(calibratedLlr)
}
