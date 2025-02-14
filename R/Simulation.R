# @file Simulation.R
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

#' Simulate (negative) controls
#'
#' @details
#' Generate point estimates given known true effect sizes and standard errors
#'
#' @param n           Number of controls to simulate.
#' @param mean        The mean of the error distribution (on the log RR scale).
#' @param sd          The standard deviation of the error distribution (on the log RR scale).
#' @param seLogRr     The standard error of the log of the relative risk. This is recycled for the
#'                    controls. The default is to sample these from a uniform distribution.
#' @param trueLogRr   The true relative risk (on the log scale) used to generate these controls.  This
#'                    is recycled for the controls.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' plotTrueAndObserved(data$logRr, data$seLogRr, data$trueLogRr)
#'
#' @export
simulateControls <- function(n = 50,
                             mean = 0,
                             sd = 0.1,
                             seLogRr = runif(n, min = 0.01, max = 0.2),
                             trueLogRr = 0) {
  theta <- rnorm(n, mean = mean, sd = sd)
  logRr <- rnorm(n, mean = trueLogRr + theta, sd = seLogRr)
  data <- data.frame(logRr = logRr, seLogRr = seLogRr, trueLogRr = trueLogRr)
  return(data)
}


#' Simulate survival data for MaxSPRT computation
#'
#' @param n                 Number of subjects.
#' @param pExposure         Probability of being in target cohort.
#' @param backgroundHazard  Background hazard (risk of the outcome per day).
#' @param tar               Time at risk for each exposure
#' @param nullMu            Null distribution mean (at log HR scale)
#' @param nullSigma         Null distribution SD (at log HR scale)
#' @param maxT              Maximum time to simulate.
#' @param looks             Number of (evenly spaced) looks at the data.
#' @param numberOfNegativeControls  Number of negative controls to simulate.
#' @param numberOfPositiveControls  Number of positive controls to simulate.
#' @param positiveControlEffectSize The true effect size of the positive controls.
#'
#' @details
#' Simulate survival data for negative and positive controls. The data provides multiple looks at data accruing over time, with
#' each look having more data than the one before. Systematic error for each outcome is drawn from the prespecified null distribution.
#'
#' The outcome IDs are assigned sequentially starting at 1, with the first IDs used for the negative controls, and the latter IDs used
#' for the positive controls.
#'
#' @examples
#' data <- simulateMaxSprtData(n = 1000)
#' head(data)
#'
#' @return
#' A data frame with 5 variables: \describe{ \item{time}{Time from index date to either the event or
#' end of observation, whichever came first} \item{outcome}{Whether the outcome occurred (1) or not (0)} \item{exposure}{Whether
#' the subject was exposed (TRUE) or not (FALSE)} \item{lookTime}{The time point when the look occurred. } \item{outcomeId}{A unique
#' identifier for data corresponding to a single outcome. Lower IDs indicate negative controls, higher IDs indicate the positive control} }
#'
#' @export
simulateMaxSprtData <- function(n = 10000,
                                pExposure = 0.5,
                                backgroundHazard = 0.001,
                                tar = 10,
                                nullMu = 0.2,
                                nullSigma = 0.2,
                                maxT = 100,
                                looks = 10,
                                numberOfNegativeControls = 50,
                                numberOfPositiveControls = 1,
                                positiveControlEffectSize = 4) {
  simulateOutcome <- function(trueEffectSize) {
    computeAtT <- function(t) {
      truncatedTime <- time
      idxTruncated <- tIndex + time > t
      truncatedTime[idxTruncated] <- t - tIndex[idxTruncated]
      truncatedOutcome <- outcome
      truncatedOutcome[idxTruncated] <- 0
      data <- data.frame(
        time = truncatedTime,
        outcome = truncatedOutcome,
        exposure = exposure
      )
      data <- data[data$time > 0, ]
      data$lookTime <- t
      return(data)
    }

    tIndex <- runif(n, 0, maxT)
    exposure <- runif(n) < pExposure
    systematicError <- rnorm(n = 1, mean = nullMu, sd = nullSigma)
    hazard <- ifelse(exposure, backgroundHazard * trueEffectSize * exp(systematicError), backgroundHazard)
    tOutcome <- rexp(n, hazard)
    outcome <- tOutcome < tar
    time <- rep(tar, n)
    time[outcome] <- tOutcome[outcome]
    t <- seq(0, maxT, length.out = looks + 1)[-1]
    results <- lapply(t, computeAtT)
    results <- do.call(rbind, results)
    return(results)
  }

  dataSets <- list()

  # Negative controls
  for (i in 1:numberOfNegativeControls) {
    dataSet <- simulateOutcome(1)
    dataSet$outcomeId <- i
    dataSets[[i]] <- dataSet
  }

  # Positive controls
  for (i in numberOfNegativeControls + (1:numberOfPositiveControls)) {
    dataSet <- simulateOutcome(positiveControlEffectSize)
    dataSet$outcomeId <- i
    dataSets[[i]] <- dataSet
  }

  maxSprtSimulationData <- do.call(rbind, dataSets)
  return(maxSprtSimulationData)
}
