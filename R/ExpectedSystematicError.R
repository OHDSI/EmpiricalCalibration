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

#' Compute the expected systematic error
#'
#' @description
#' \code{calibrateP} computes calibrated p-values using the fitted null distribution
#'
#' @details
#' This function computes a calibrated two-sided p-value as described in Schuemie et al (2014).
#'
#' @param null      An object of class \code{null} created using the \code{fitNull} function or an
#'                  object of class \code{mcmcNull} created using the \code{fitMcmcNull} function.
#' @param ...       Any additional parameters (currently none).
#'
#' @return
#' The two-sided calibrated p-value.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' computeExpectedSystematicError(null)
#'
#' @export
computeExpectedSystematicError <- function(null) {
  UseMethod("computeExpectedSystematicError") 
}

#' @export
computeExpectedSystematicError.null <- function(null) {
  if (null[1] == 0 && null[2] == 0) {
    return(0)
  }
  result <- closedFormIntegeralAbsolute(null[1], null[2])
  names(result) <- NULL
  return(result)
}

#' @export
computeExpectedSystematicError.mcmcNull <- function(null) {
  chain <- attr(null, "mcmc")$chain
  dist <- apply(chain, 1, function(x) closedFormIntegeralAbsolute(x[1], 1 / sqrt(x[2])))
  # dist <- apply(chain, 1, function(x) closedFormIntegeralAbsolute(x[1], x[2]))
  result <- quantile(dist, c(0.5, 0.025, 0.975))
  result <- data.frame(expectedSystematicError = result[1],
                       lb95ci = result[2],
                       ub95ci = result[3])
  row.names(result) <- NULL
  return(result)
}
