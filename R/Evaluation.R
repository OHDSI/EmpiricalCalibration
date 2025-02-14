# @file Evaluation.R
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

#' Evaluate confidence interval calibration
#'
#' @description
#' \code{evaluateCiCalibration} performs a leave-one-out cross-validation to evaluate the calibration
#' confidence intervals.
#'
#' @details
#' The empirical calibration is performed using a leave-one-out design: The confidence interval of an
#' effect is computed by fitting a null using all other controls.
#'
#' @param logRr                  A numeric vector of effect estimates on the log scale.
#' @param seLogRr                The standard error of the log of the effect estimates. Hint: often the
#'                               standard error = (log(<lower bound 95 percent confidence interval>) -
#'                               log(<effect estimate>))/qnorm(0.025).
#' @param trueLogRr              The true log relative risk.
#' @param strata                 Variable used to stratify the plot. Set \code{strata = NULL} for no
#'                               stratification.
#' @param crossValidationGroup   What should be the unit for the cross-validation? By default the unit
#'                               is a single control, but a different grouping can be provided, for
#'                               example linking a negative control to synthetic positive controls
#'                               derived from that negative control.
#' @param legacy                 If true, a legacy error model will be fitted, meaning standard
#'                               deviation is linear on the log scale. If false, standard deviation
#'                               is assumed to be simply linear.
#'
#' @return
#' A data frame specifying the coverage per strata (usually true effect size) for a wide range of widths
#' of the confidence interval. The result also includes the fraction of estimates that was below and above
#' the confidence interval.
#'
#' @examples
#' \dontrun{
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' eval <- evaluateCiCalibration(data$logRr, data$seLogRr, data$trueLogRr)
#' }
#' @export
evaluateCiCalibration <- function(logRr,
                                  seLogRr,
                                  trueLogRr,
                                  strata = as.factor(trueLogRr),
                                  crossValidationGroup = 1:length(logRr),
                                  legacy = FALSE) {
  if (!is.null(strata) && !is.factor(strata)) {
    stop("Strata argument should be a factor (or null)")
  }
  if (is.null(strata)) {
    strata <- as.factor(-1)
  }
  data <- data.frame(
    logRr = logRr,
    seLogRr = seLogRr,
    trueLogRr = trueLogRr,
    strata = strata,
    crossValidationGroup = crossValidationGroup
  )
  if (any(is.infinite(data$seLogRr))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting error model")
    data <- data[!is.infinite(seLogRr), ]
  }
  if (any(is.infinite(data$logRr))) {
    warning("Estimate(s) with infinite logRr detected. Removing before fitting error model")
    data <- data[!is.infinite(logRr), ]
  }
  if (any(is.na(data$seLogRr))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting error model")
    data <- data[!is.na(seLogRr), ]
  }
  if (any(is.na(data$logRr))) {
    warning("Estimate(s) with NA logRr detected. Removing before fitting error model")
    data <- data[!is.na(logRr), ]
  }
  computeCoverage <- function(j, subResult, dataLeftOut, model) {
    subset <- dataLeftOut[dataLeftOut$strata == subResult$strata[j], ]
    if (nrow(subset) == 0) {
      return(0)
    }
    # writeLines(paste0("ciWidth: ", subResult$ciWidth[j], ", strata: ", subResult$strata[j], ", model: ", paste(model, collapse = ",")))
    # writeLines(paste0("ciWidth: ", subResult$ciWidth[j], ", logRr: ", paste(subset$logRr, collapse = ","), ", seLogRr:", paste(subset$seLogRr, collapse = ","), ", model: ", paste(model, collapse = ",")))
    ci <- EmpiricalCalibration::calibrateConfidenceInterval(
      logRr = subset$logRr,
      seLogRr = subset$seLogRr,
      ciWidth = subResult$ciWidth[j],
      model = model
    )
    ci$logLb95Rr[is.na(ci$logLb95Rr)] <- 0
    ci$logUb95Rr[is.na(ci$logUb95Rr)] <- 999
    below <- sum(subset$trueLogRr < ci$logLb95Rr)
    within <- sum(ci$logLb95Rr <= subset$trueLogRr & ci$logUb95Rr >= subset$trueLogRr)
    above <- sum(subset$trueLogRr > ci$logUb95Rr)

    return(c(below, within, above))
  }

  computeTheoreticalCoverage <- function(j, subResult, dataLeftOut) {
    subset <- dataLeftOut[dataLeftOut$strata == subResult$strata[j], ]
    ciWidth <- subResult$ciWidth[j]
    logLb95Rr <- subset$logRr + qnorm((1 - ciWidth) / 2) * subset$seLogRr
    logUb95Rr <- subset$logRr - qnorm((1 - ciWidth) / 2) * subset$seLogRr
    below <- sum(subset$trueLogRr < logLb95Rr)
    within <- sum(subset$trueLogRr >= logLb95Rr & subset$trueLogRr <= logUb95Rr)
    above <- sum(subset$trueLogRr > logUb95Rr)
    return(c(below, within, above))
  }

  computeLooCoverage <- function(leaveOutGroup, data, legacy) {
    dataLeaveOneOut <- data[data$crossValidationGroup != leaveOutGroup, ]
    dataLeftOut <- data[data$crossValidationGroup == leaveOutGroup, ]
    if (nrow(dataLeaveOneOut) == 0 || nrow(dataLeftOut) == 0) {
      return(data.frame())
    }

    model <- fitSystematicErrorModel(
      logRr = dataLeaveOneOut$logRr,
      seLogRr = dataLeaveOneOut$seLogRr,
      trueLogRr = dataLeaveOneOut$trueLogRr,
      estimateCovarianceMatrix = FALSE,
      legacy = legacy
    )
    strata <- unique(dataLeftOut$strata)
    ciWidth <- seq(0.01, 0.99, by = 0.01)
    subResult <- expand.grid(strata, ciWidth)
    names(subResult) <- c("strata", "ciWidth")
    coverage <- sapply(1:nrow(subResult), computeCoverage, subResult = subResult, dataLeftOut = dataLeftOut, model = model)
    subResult$below <- coverage[1, ]
    subResult$within <- coverage[2, ]
    subResult$above <- coverage[3, ]
    theoreticalCoverage <- sapply(1:nrow(subResult), computeTheoreticalCoverage, subResult = subResult, dataLeftOut = dataLeftOut)
    subResult$theoreticalBelow <- theoreticalCoverage[1, ]
    subResult$theoreticalWithin <- theoreticalCoverage[2, ]
    subResult$theoreticalAbove <- theoreticalCoverage[3, ]
    return(subResult)
  }
  writeLines("Fitting error models within leave-one-out cross-validation")
  coverages <- lapply(unique(data$crossValidationGroup), computeLooCoverage, data = data, legacy = legacy)
  coverage <- do.call("rbind", coverages)
  data$count <- 1
  counts <- aggregate(count ~ strata, data = data, sum)
  belowCali <- aggregate(below ~ strata + ciWidth, data = coverage, sum)
  belowCali <- merge(belowCali, counts, by = "strata")
  belowCali$coverage <- belowCali$below / belowCali$count
  belowCali$label <- "Below confidence interval"
  belowCali$type <- "Calibrated"
  withinCali <- aggregate(within ~ strata + ciWidth, data = coverage, sum)
  withinCali <- merge(withinCali, counts, by = "strata")
  withinCali$coverage <- withinCali$within / withinCali$count
  withinCali$label <- "Within confidence interval"
  withinCali$type <- "Calibrated"
  aboveCali <- aggregate(above ~ strata + ciWidth, data = coverage, sum)
  aboveCali <- merge(aboveCali, counts, by = "strata")
  aboveCali$coverage <- aboveCali$above / aboveCali$count
  aboveCali$label <- "Above confidence interval"
  aboveCali$type <- "Calibrated"
  belowUncali <- aggregate(theoreticalBelow ~ strata + ciWidth, data = coverage, sum)
  belowUncali <- merge(belowUncali, counts, by = "strata")
  belowUncali$coverage <- belowUncali$theoreticalBelow / belowUncali$count
  belowUncali$label <- "Below confidence interval"
  belowUncali$type <- "Uncalibrated"
  withinUncali <- aggregate(theoreticalWithin ~ strata + ciWidth, data = coverage, sum)
  withinUncali <- merge(withinUncali, counts, by = "strata")
  withinUncali$coverage <- withinUncali$theoreticalWithin / withinUncali$count
  withinUncali$label <- "Within confidence interval"
  withinUncali$type <- "Uncalibrated"
  aboveUncali <- aggregate(theoreticalAbove ~ strata + ciWidth, data = coverage, sum)
  aboveUncali <- merge(aboveUncali, counts, by = "strata")
  aboveUncali$coverage <- aboveUncali$theoreticalAbove / aboveUncali$count
  aboveUncali$label <- "Above confidence interval"
  aboveUncali$type <- "Uncalibrated"
  vizData <- rbind(
    belowCali[, c("strata", "label", "type", "ciWidth", "coverage")],
    withinCali[, c("strata", "label", "type", "ciWidth", "coverage")],
    aboveCali[, c("strata", "label", "type", "ciWidth", "coverage")],
    belowUncali[, c("strata", "label", "type", "ciWidth", "coverage")],
    withinUncali[, c("strata", "label", "type", "ciWidth", "coverage")],
    aboveUncali[, c("strata", "label", "type", "ciWidth", "coverage")]
  )
  names(vizData)[names(vizData) == "type"] <- "Confidence interval calculation"
  vizData$trueRr <- as.factor(exp(as.numeric(as.character(vizData$strata))))
  vizData$strata <- NULL
  return(vizData)
}
