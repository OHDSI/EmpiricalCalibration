# @file ExpectedSystematicError.R
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

closedFormIntegral <- function(x, mu, sigma) {
  mu * pnorm(x, mu, sigma) - 1 - sigma^2 * dnorm(x, mu, sigma)
}

closedFormIntegeralAbsolute <- function(mu, sigma) {
  closedFormIntegral(Inf, mu = mu, sigma = sigma) - 2 * closedFormIntegral(0, mu = mu, sigma = sigma) + closedFormIntegral(-Inf, mu = mu, sigma = sigma)
}

#' Compute the expected absolute systematic error
#'
#' @description
#' For a random study estimate, what is the expected value of the absolute systematic error?
#' Provides a single summary value for a null distribution. The expected systematic error of a null
#' distribution is equal to its mean (mu), and is insensitive to the spread of the null distribution (sigma).
#'
#' Taking the absolute value of the expected systematic error we can express both mean and spread of the
#' estimated null distribution.
#'
#' Note that the expected absolute systematic error is expressed on the scale of the estimates used to fit the
#' null (for example log relative risk or risk difference), so values are not comparable across scales.
#'
#' @param null      An object of class `null` created using the `fitNull` function or an
#'                  object of class `mcmcNull` created using the `fitMcmcNull` function.
#' @param alpha     The expected type I error for computing the credible interval.
#'
#' @return
#' The expected absolute systematic error. If the provided `null` argument is of type `mcmcNull`,
#' the credible interval (defined by `alpha`) is also returned.
#'
#' @seealso [compareEase()] for comparing the expected absolute systematic error of two sets of estimates for the same negative controls.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' computeExpectedAbsoluteSystematicError(null)
#' @export
computeExpectedAbsoluteSystematicError <- function(null, alpha = 0.05) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(null, classes = "null", add = errorMessages)
  checkmate::assertNumber(alpha, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  UseMethod("computeExpectedAbsoluteSystematicError")
}

#' @export
computeExpectedAbsoluteSystematicError.null <- function(null, alpha = 0.05) {
  if (null[1] == 0 && null[2] == 0) {
    return(0)
  }
  result <- closedFormIntegeralAbsolute(null[1], null[2])
  names(result) <- NULL
  return(result)
}

#' @export
computeExpectedAbsoluteSystematicError.mcmcNull <- function(null, alpha = 0.05) {
  if (is.na(null[1])) {
    result <- data.frame(
      ease = as.numeric(NA),
      ciLb = as.numeric(NA),
      ciUb = as.numeric(NA)
    )
  } else {
    chain <- attr(null, "mcmc")$chain
    dist <- apply(chain, 1, function(x) closedFormIntegeralAbsolute(x[1], 1 / sqrt(x[2])))
    result <- quantile(dist, c(0.5, alpha / 2, 1 - (alpha / 2)))
    result <- data.frame(
      ease = result[1],
      ciLb = result[2],
      ciUb = result[3]
    )
    row.names(result) <- NULL
  }
  return(result)
}

computeEaseBoostrap <- function(estimate, se, alpha, sampleSize) {
  boostrap <- function(i) {
    idx <- sample.int(length(estimate), length(estimate), replace = TRUE)
    null <- fitNull(estimate[idx], se[idx])
    return(computeExpectedAbsoluteSystematicError(null))
  }
  sample <- sapply(1:sampleSize, boostrap)
  result <- quantile(sample, c(0.5, 0.025, 0.975))
  result <- data.frame(
    ease = result[1],
    ciLb = result[2],
    ciUb = result[3]
  )
  row.names(result) <- NULL
  return(result)
}

#' Compare EASE of correlated sets of estimates
#'
#' @details
#' Compare the expected absolute systematic error (EASE) of two sets of estimates for the same set of negative controls.
#'
#' Important: the two sets of estimates (estimate1 + se1 and estimate2 + se2) should be in identical order, so that for
#' example the first item in each vector corresponds to the same negative control.
#'
#' @param estimate1  A numeric vector of effect estimates generated using the first method (for example log relative risks).
#' @param se1        The standard error of the effect estimates generated using the first method.
#' @param estimate2  A numeric vector of effect estimates generated using the second method (on the same scale as `estimate1`).
#' @param se2        The standard error of the effect estimates generated using the second method.
#' @param alpha      The expected type I error for computing confidence intervals and p-values.
#' @param sampleSize The number of samples in the bootstraps.
#' @param logRr1     DEPRECATED. Use `estimate1` instead.
#' @param seLogRr1   DEPRECATED. Use `se1` instead.
#' @param logRr2     DEPRECATED. Use `estimate2` instead.
#' @param seLogRr2   DEPRECATED. Use `se2` instead.
#'
#' @examples
#' # Simulate results of first method:
#' ncs1 <- simulateControls(n = 50)
#'
#' # Simulate second method to be more biased:
#' ncs2 <- ncs1
#' ncs2$estimate <- ncs2$estimate + rnorm(nrow(ncs2), mean = 0.1, sd = 0.1)
#'
#' delta <- compareEase(
#'   estimate1 = ncs1$estimate,
#'   se1 = ncs1$se,
#'   estimate2 = ncs2$estimate,
#'   se2 = ncs2$se
#' )
#' delta
#' attr(delta, "ease1")
#' attr(delta, "ease2")
#' @return
#' A data frame with 4 columns: the point estimate, confidence interval lower bound, and upper bound for the difference
#' between EASE in the two sets of negative controls, and a p value against the null hypothesis that the EASE is the
#' same for the two sets.
#'
#' The data frame has two attributes: ease1 and ease2, providing the EASE estimates (and confidence intervals) for the
#' two sets, computed using bootstrapping. Note that these estimates may somewhat different from those generated using
#' [computeExpectedAbsoluteSystematicError()], because a different approach is used to compute the confidence
#' interval. The approach used here will more closely align with the computation of the difference in EASE.
#'
#' @export
compareEase <- function(estimate1,
                        se1,
                        estimate2,
                        se2,
                        alpha = 0.05,
                        sampleSize = 1000,
                        logRr1 = NULL,
                        seLogRr1 = NULL,
                        logRr2 = NULL,
                        seLogRr2 = NULL) {
  if (!is.null(logRr1)) {
    warnDeprecatedArguments(c("logRr1", "seLogRr1", "logRr2", "seLogRr2"), c("estimate1", "se1", "estimate2", "se2"))
    estimate1 <- logRr1
    se1 <- seLogRr1
    estimate2 <- logRr2
    se2 <- seLogRr2
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimate1, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se1, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(estimate2, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se2, min.len = 1, add = errorMessages)
  checkmate::assertNumber(alpha, add = errorMessages)
  checkmate::assertCount(sampleSize, add = errorMessages)
  checkmate::assertNumeric(logRr1, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr1, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRr2, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr2, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (length(unique(c(length(estimate1), length(se1), length(estimate2), length(se2)))) != 1) {
    stop("Arguments estimate1, se1, estimate2, and se2 should be of equal length.")
  }
  
  ease1 <- computeEaseBoostrap(estimate = estimate1, se = se1, alpha = alpha, sampleSize = sampleSize)
  ease2 <- computeEaseBoostrap(estimate = estimate2, se = se2, alpha = alpha, sampleSize = sampleSize)
  delta <- ease1$ease - ease2$ease
  sampleDeltaEase <- function(i) {
    idx <- sample.int(length(estimate1), length(estimate1), replace = TRUE)
    null1 <- fitNull(estimate1[idx], se1[idx],)
    ease1 <- computeExpectedAbsoluteSystematicError(null1)
    null2 <- fitNull(estimate2[idx], se2[idx])
    ease2 <- computeExpectedAbsoluteSystematicError(null2)
    return(ease1 - ease2)
  }
  sample <- sapply(1:sampleSize, sampleDeltaEase)
  deltaCi <- quantile(sample, c(0.025, 0.975))
  # Tolerance is in units of the estimates, so needs to follow the magnitude of the standard errors:
  epsilon <- 0.0001 
  if (delta > 0) {
    deltaP <- mean(sample < epsilon)
  } else {
    deltaP <- mean(sample > -epsilon)
  }
  delta <- data.frame(
    delta = delta,
    ciLb = deltaCi[1],
    ciUb = deltaCi[2],
    p = deltaP
  )
  rownames(delta) <- NULL
  attr(delta, "ease1") <- ease1
  attr(delta, "ease2") <- ease2
  return(delta)
}
