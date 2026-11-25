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
#' `evaluateCiCalibration` performs a leave-one-out cross-validation to evaluate the calibration
#' confidence intervals.
#'
#' @details
#' The empirical calibration is performed using a leave-one-out design: The confidence interval of an
#' effect is computed by fitting a null using all other controls.
#'
#' @param estimate               A numeric vector of effect estimates. On the log relative risk scale
#'                               these are log relative risks.
#' @param se                     The standard error of the effect estimates. Hint: on the log relative
#'                               risk scale often the standard error = (log(<lower bound 95 percent
#'                               confidence interval>) - log(<effect estimate>))/qnorm(0.025).
#' @param trueEffect             The true effect (on the same scale as `estimate`).
#' @template Scale
#' @param strata                 Variable used to stratify the plot. Set `strata = NULL` for no
#'                               stratification.
#' @param crossValidationGroup   What should be the unit for the cross-validation? By default the unit
#'                               is a single control, but a different grouping can be provided, for
#'                               example linking a negative control to synthetic positive controls
#'                               derived from that negative control.
#' @param legacy                 If true, a legacy error model will be fitted, meaning standard
#'                               deviation is linear on the log scale. If false, standard deviation
#'                               is assumed to be simply linear. The legacy model is not supported for the
#'                               risk difference scale.
#' @param logRr                  DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr                DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#' @param trueLogRr              DEPRECATED. Use `trueEffect` with `scale = "log relative risk"` instead.
#'
#' @return
#' A data frame specifying the coverage per strata (usually true effect size) for a wide range of widths
#' of the confidence interval. The result also includes the fraction of estimates that was below and above
#' the confidence interval. The `trueEffect` column holds the true effect per stratum: the relative risk
#' (exponentiated true log relative risk) for `scale = "log relative risk"`, and the true risk difference for
#' `scale = "risk difference"`. For backwards compatibility the deprecated column `trueRr` holds the same values.
#'
#' @examples
#' \dontrun{
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' eval <- evaluateCiCalibration(data$estimate, data$se, data$trueEffect)
#' }
#' @export
evaluateCiCalibration <- function(estimate,
                                  se,
                                  trueEffect,
                                  scale = "log relative risk",
                                  strata = as.factor(trueEffect),
                                  crossValidationGroup = 1:length(estimate),
                                  legacy = FALSE,
                                  logRr = NULL,
                                  seLogRr = NULL,
                                  trueLogRr = NULL) {
  if (!is.null(logRr)) {
    warnDeprecatedArguments(c("logRr", "seLogRr", "trueLogRr"), c("estimate", "se", "trueEffect"))
    estimate <- logRr
    se <- seLogRr
    trueEffect <- trueLogRr
    scale <- LOG_RELATIVE_RISK
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimate, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(se, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(trueEffect, min.len = 1, add = errorMessages)
  checkmate::assertChoice(scale, choices = c(LOG_RELATIVE_RISK, RISK_DIFFERENCE), add = errorMessages)
  checkmate::assertFactor(strata, null.ok = TRUE, add = errorMessages)
  checkmate::assertVector(crossValidationGroup, len = length(estimate), add = errorMessages)
  checkmate::assertLogical(legacy, len = 1, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(trueLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  if (!is.null(strata) && !is.factor(strata)) {
    stop("Strata argument should be a factor (or null)")
  }
  if (is.null(strata)) {
    strata <- as.factor(-1)
  }
  data <- data.frame(
    estimate = estimate,
    se = se,
    trueEffect = trueEffect,
    strata = strata,
    crossValidationGroup = crossValidationGroup
  )
  if (any(is.infinite(data$se))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting error model")
    data <- data[!is.infinite(data$se), ]
  }
  if (any(is.infinite(data$estimate))) {
    warning("Estimate(s) with infinite estimate detected. Removing before fitting error model")
    data <- data[!is.infinite(data$estimate), ]
  }
  if (any(is.na(data$se))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting error model")
    data <- data[!is.na(data$se), ]
  }
  if (any(is.na(data$estimate))) {
    warning("Estimate(s) with NA estimate detected. Removing before fitting error model")
    data <- data[!is.na(data$estimate), ]
  }
  computeCoverage <- function(j, subResult, dataLeftOut, model) {
    subset <- dataLeftOut[dataLeftOut$strata == subResult$strata[j], ]
    if (nrow(subset) == 0) {
      return(0)
    }
    ci <- EmpiricalCalibration::calibrateConfidenceInterval(
      estimate = subset$estimate,
      se = subset$se,
      ciWidth = subResult$ciWidth[j],
      model = model
    )
    ci$lb[is.na(ci$lb)] <- 0
    ci$ub[is.na(ci$ub)] <- 999
    below <- sum(subset$trueEffect < ci$lb)
    within <- sum(ci$lb <= subset$trueEffect & ci$ub >= subset$trueEffect)
    above <- sum(subset$trueEffect > ci$ub)

    return(c(below, within, above))
  }

  computeTheoreticalCoverage <- function(j, subResult, dataLeftOut) {
    subset <- dataLeftOut[dataLeftOut$strata == subResult$strata[j], ]
    ciWidth <- subResult$ciWidth[j]
    lb <- subset$estimate + qnorm((1 - ciWidth) / 2) * subset$se
    ub <- subset$estimate - qnorm((1 - ciWidth) / 2) * subset$se
    below <- sum(subset$trueEffect < lb)
    within <- sum(subset$trueEffect >= lb & subset$trueEffect <= ub)
    above <- sum(subset$trueEffect > ub)
    return(c(below, within, above))
  }

  computeLooCoverage <- function(leaveOutGroup, data, legacy) {
    dataLeaveOneOut <- data[data$crossValidationGroup != leaveOutGroup, ]
    dataLeftOut <- data[data$crossValidationGroup == leaveOutGroup, ]
    if (nrow(dataLeaveOneOut) == 0 || nrow(dataLeftOut) == 0) {
      return(data.frame())
    }

    model <- fitSystematicErrorModel(
      estimate = dataLeaveOneOut$estimate,
      se = dataLeaveOneOut$se,
      trueEffect = dataLeaveOneOut$trueEffect,
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
  trueEffect <- as.numeric(as.character(vizData$strata))
  if (scale == LOG_RELATIVE_RISK) {
    trueEffect <- exp(trueEffect)
  }
  vizData$trueEffect <- as.factor(trueEffect)
  vizData$trueRr <- vizData$trueEffect
  vizData$strata <- NULL
  return(vizData)
}
