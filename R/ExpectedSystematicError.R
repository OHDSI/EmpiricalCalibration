# @file ExpectedSystematicError.R
#
# Copyright 2020 Observational Health Data Sciences and Informatics
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

closedFormIntegral <- function(x, mu ,sigma) {
  mu * pnorm(x, mu, sigma) - 1 - sigma^2 * dnorm(x, mu, sigma)
}

closedFormIntegeralAbsolute <- function(mu, sigma) {
  closedFormIntegral(Inf, mu = mu, sigma = sigma) - 2*closedFormIntegral(0, mu = mu, sigma = sigma) + closedFormIntegral(-Inf, mu = mu, sigma = sigma)
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
#' @param null      An object of class \code{null} created using the \code{fitNull} function or an
#'                  object of class \code{mcmcNull} created using the \code{fitMcmcNull} function.
#' @param alpha     The expected type I error for computing the credible interval.
#'
#' @return
#' The expected absolute systematic error. If the provided \code{null} argument is of type \code{mcmcNull},
#' the credible interval (defined by \code{alpha}) is also returned.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' computeExpectedSystematicError(null)
#'
#' @export
computeExpectedAbsoluteSystematicError <- function(null, alpha = 0.05) {
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
  chain <- attr(null, "mcmc")$chain
  dist <- apply(chain, 1, function(x) closedFormIntegeralAbsolute(x[1], 1 / sqrt(x[2])))
  result <- quantile(dist, c(0.5, alpha / 2, 1 - (alpha / 2)))
  result <- data.frame(ease = result[1],
                       ciLb = result[2],
                       ciUb = result[3])
  row.names(result) <- NULL
  return(result)
}
