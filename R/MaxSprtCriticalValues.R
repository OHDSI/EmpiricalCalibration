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


#' Compute critical values for Poisson data
#'
#' @description
#' Obtains critical values for the group and continuous sequential MaxSPRT test with Poisson data,
#' using a Wald type upper boundary, which is flat with respect to the likelihood ratio function,
#' and with a pre-specified upper limit on the sample size.
#'
#' It is often not possible to select a critical value that corresponds to the exact alpha specified.
#' Instead, this function will select the least conservative critical value having an alpha below
#' the one specified, so the sequential analysis is conservative.
#'
#' This function is a re-implementation of the \code{CV.Poisson} function in the \code{Sequential}
#' package, using Monte-Carlo.
#'
#' @param groupSizes    Vector containing the expected number of events under H0 for each test.
#' @param minimumEvents The minimum number of events needed before the null hypothesis can be rejected.
#' @param alpha         The significance level, or the type 1 error probability, which is the probability
#'                      of rejecting the null hypothesis when it is true.
#' @param sampleSize    Sample size for the Monte-Carlo simulations.
#' @param nullMean      The mean of the empirical null distribution.
#' @param nullSd        The standard deviation of the empirical null distribution.
#'
#' @return
#' The computed critical value. The 'alpha' attribute of the result indicates the selected alpha.
#'
#' @examples
#' groupSizes <- rep(1, 10)
#' computeCvPoisson(groupSizes)
#'
#' @export
computeCvPoisson <- function(groupSizes, minimumEvents = 1, alpha = 0.05, sampleSize = 1e6, nullMean = 0, nullSd = 0) {
  if (any(groupSizes < 0)) {
    stop("Group sizes should be positive")
  }
  if (minimumEvents < 0 || minimumEvents != round(minimumEvents)) {
    stop("Minimum events should be a positive integer")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("Alpha should be a number between 0 and 1")
  }
  if (sampleSize < 0 || sampleSize != round(sampleSize)) {
    stop("sampleSize should be a positive integer")
  }

  values <- samplePoissonMaxLrr(
    groupSizes = groupSizes,
    minimumEvents = minimumEvents,
    sampleSize = sampleSize,
    nullMean = nullMean,
    nullSd = nullSd
  )
  if (max(values) == 0) {
    warning("Unable to find a critical value that could lead to rejection of the null")
    result <- Inf
    attr(result, "alpha") <- 0
    return(result)
  }
  values <- values[order(-values)]
  alphas <- 1:sampleSize / sampleSize
  idx <- !duplicated(values)
  alphas <- alphas[idx]
  values <- values[idx]
  pick <- max(which(alphas < alpha))
  message(sprintf("Selected alpha: %0.3f (least conservative value below %s)", alphas[pick], alpha))
  result <- values[pick]
  attr(result, "alpha") <- alphas[pick]
  return(result)
}

#' Compute critical values for Binomial data
#'
#' @description
#' Obtains critical values for the group and continuous sequential MaxSPRT test with Binomial data,
#' using a Wald type upper boundary, which is flat with respect to the likelihood ratio function,
#' and with a pre-specified upper limit on the sample size.
#'
#' It is often not possible to select a critical value that corresponds to the exact alpha specified.
#' Instead, this function will select the least conservative critical value having an alpha below
#' the one specified, so the sequential analysis is conservative.
#'
#' This function is a re-implementation of the \code{CV.Binomial} function in the \code{Sequential}
#' package, using Monte-Carlo.
#'
#' @param groupSizes    Vector containing the expected number of events under H0 for each test.
#' @param z             For a matched case-control analysis, z is the number of controls matched to
#'                      each case under the null hypothesis. For a self-controlled analysis, z is the
#'                      control time divided by the time at risk.
#' @param minimumEvents The minimum number of events needed before the null hypothesis can be rejected.
#' @param alpha         The significance level, or the type 1 error probability, which is the probability
#'                      of rejecting the null hypothesis when it is true.
#' @param sampleSize    Sample size for the Monte-Carlo simulations.
#' @param nullMean      The mean of the empirical null distribution.
#' @param nullSd        The standard deviation of the empirical null distribution.
#'
#' @return
#' The computed critical value. The 'alpha' attribute of the result indicates the selected alpha.
#'
#' @examples
#' groupSizes <- rep(1, 10)
#' computeCvBinomial(groupSizes, z = 4)
#'
#' @export
computeCvBinomial <- function(groupSizes, z, minimumEvents = 1, alpha = 0.05, sampleSize = 1e6, nullMean = 0, nullSd = 0) {
  if (any(groupSizes < 0)) {
    stop("Group sizes should be positive")
  }
  if (z < 0) {
    stop("Z should be positive")
  }
  if (minimumEvents < 0 || minimumEvents != round(minimumEvents)) {
    stop("Minimum events should be a positive integer")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("Alpha should be a number between 0 and 1")
  }
  if (sampleSize < 0 || sampleSize != round(sampleSize)) {
    stop("sampleSize should be a positive integer")
  }

  p <- 1 / (1 + z)
  values <- sampleBinomialMaxLrr(
    groupSizes = groupSizes,
    p = p,
    minimumEvents = minimumEvents,
    sampleSize = sampleSize,
    nullMean = nullMean,
    nullSd = nullSd
  )
  if (max(values) == 0) {
    warning("Unable to find a critical value that could lead to rejection of the null")
    result <- Inf
    attr(result, "alpha") <- 0
    return(result)
  }
  values <- values[order(-values)]
  alphas <- 1:sampleSize / sampleSize
  idx <- !duplicated(values)
  alphas <- alphas[idx]
  values <- values[idx]
  pick <- max(which(alphas < alpha))
  message(sprintf("Selected alpha: %0.3f (least conservative value below %s)", alphas[pick], alpha))
  result <- values[pick]
  attr(result, "alpha") <- alphas[pick]
  return(result)
}

#' Compute critical values for Poisson regression data
#'
#' @description
#' Obtains critical values for the group and continuous sequential MaxSPRT test with Poisson regression data,
#' using a Wald type upper boundary, which is flat with respect to the likelihood ratio function,
#' and with a pre-specified upper limit on the sample size.
#'
#' It is often not possible to select a critical value that corresponds to the exact alpha specified.
#' Instead, this function will select the least conservative critical value having an alpha below
#' the one specified, so the sequential analysis is conservative.
#'
#' @param groupSizes    Vector containing the expected number of events under H0 for each test.
#' @param z             The control time divided by the time at risk.
#' @param minimumEvents The minimum number of events needed before the null hypothesis can be rejected.
#' @param alpha         The significance level, or the type 1 error probability, which is the probability
#'                      of rejecting the null hypothesis when it is true.
#' @param sampleSize    Sample size for the Monte-Carlo simulations.
#'
#' @return
#' The computed critical value. The 'alpha' attribute of the result indicates the selected alpha.
#'
#' @examples
#' groupSizes <- rep(1, 10)
#' computeCvPoissonRegression(groupSizes, z = 4)
#'
#' @export
computeCvPoissonRegression <- function(groupSizes, z, minimumEvents = 1, alpha = 0.05, sampleSize = 1e6) {
  if (any(groupSizes < 0)) {
    stop("Group sizes should be positive")
  }
  if (minimumEvents < 0 || minimumEvents != round(minimumEvents)) {
    stop("Minimum events should be a positive integer")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("Alpha should be a number between 0 and 1")
  }
  if (sampleSize < 0 || sampleSize != round(sampleSize)) {
    stop("sampleSize should be a positive integer")
  }
  if (z < 0) {
    stop("Z should be positive")
  }

  values <- samplePoissonRegressionMaxLrr(
    groupSizes = groupSizes,
    z = z,
    minimumEvents = minimumEvents,
    sampleSize = sampleSize
  )
  values <- values[order(-values)]
  alphas <- 1:sampleSize / sampleSize
  idx <- !duplicated(values)
  alphas <- alphas[idx]
  values <- values[idx]
  pick <- max(which(alphas < alpha))
  message(sprintf("Selected alpha: %0.3f (least conservative value below %s)", alphas[pick], alpha))
  result <- values[pick]
  attr(result, "alpha") <- alphas[pick]
  return(result)
}
