# @file Plots.R
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

#' Create a forest plot
#'
#' @description
#' `plotForest` creates a forest plot of effect size estimates.
#'
#' @details
#' Creates a forest plot of effect size estimates. For the log relative risk scale the estimates are shown as
#' ratios on a log axis, and estimates that are significantly different from 1 (alpha = 0.05) are marked in
#' orange, others are marked in blue. For the risk difference scale estimates are shown on a linear axis, and
#' estimates that are significantly different from 0 are marked in orange.
#'
#'
#' @param estimate   A numeric vector of effect estimates. On the log relative risk scale these are
#'                   log relative risks.
#' @param se         The standard error of the effect estimates. Hint: on the log relative risk scale
#'                   often the standard error = (log(<lower bound 95 percent confidence interval>) -
#'                   log(<effect estimate>))/qnorm(0.025)
#' @param names      A vector containing the names of the drugs or outcomes
#' @template Scale
#' @param xLabel     The label on the x-axis: the name of the effect estimate. Defaults to
#'                   "Relative risk" or "Risk difference", depending on the scale.
#' @param title      Optional: the main title for the plot
#' @param logRr      DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr    DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                   function `ggsave` in the ggplot2 package for supported file formats.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' plotForest(negatives$logRr, negatives$seLogRr, negatives$drugName)
#'
#' @export
plotForest <- function(estimate,
                       se,
                       names,
                       scale = "log relative risk",
                       xLabel = NULL,
                       title = NULL,
                       fileName = NULL,
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
  checkmate::assertVector(names, len = length(estimate), add = errorMessages)
  checkmate::assertChoice(scale, choices = c(LOG_RELATIVE_RISK, RISK_DIFFERENCE), add = errorMessages)
  checkmate::assertString(xLabel, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(xLabel)) {
    xLabel <- if (scale == LOG_RELATIVE_RISK) "Relative risk" else "Risk difference"
  }
  theme <- ggplot2::element_text(colour = "#000000", size = 6)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 5, hjust = 1)
  col <- c(rgb(0, 0, 0.8, alpha = 1), rgb(0.8, 0.4, 0, alpha = 1))
  colFill <- c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0.4, 0, alpha = 0.5))
  data <- data.frame(
    DRUG_NAME = as.factor(names),
    estimate = estimate,
    lb95 = estimate + qnorm(0.025) * se,
    ub95 = estimate + qnorm(0.975) * se
  )
  data$significant <- data$lb95 > 0 | data$ub95 < 0
  data$DRUG_NAME <- factor(data$DRUG_NAME, levels = rev(levels(data$DRUG_NAME)))
  if (scale == LOG_RELATIVE_RISK) {
    breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
    labels <- breaks
    nullValue <- 1
    data$estimate <- exp(data$estimate)
    data$lb95 <- exp(data$lb95)
    data$ub95 <- exp(data$ub95)
    yScale <- ggplot2::scale_y_continuous(xLabel, trans = "log10", breaks = breaks, labels = labels)
    limits <- c(0.25, 10)
  } else {
    nullValue <- 0
    yScale <- ggplot2::scale_y_continuous(xLabel)
    limits <- NULL
  }
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$DRUG_NAME,
      y = .data$estimate,
      ymin = .data$lb95,
      ymax = .data$ub95,
      colour = .data$significant,
      fill = .data$significant
    )
  )
  if (scale == LOG_RELATIVE_RISK) {
    plot <- plot + ggplot2::geom_hline(yintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2)
  }
  plot <- plot +
    ggplot2::geom_hline(yintercept = nullValue, size = 0.5) +
    ggplot2::geom_pointrange(shape = 23) +
    ggplot2::scale_colour_manual(values = col) +
    ggplot2::scale_fill_manual(values = colFill) +
    ggplot2::coord_flip(ylim = limits) +
    yScale +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = "#EEEEEE"),
      axis.ticks = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.key = ggplot2::element_blank(),
      strip.text.x = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 5, height = 2.5 + length(estimate) * 0.8, dpi = 400)
  }
  return(plot)
}

logRrtoSE <- function(logRr, alpha, mu, sigma) {
  phi <- (mu - logRr)^2 / qnorm(alpha / 2)^2 - sigma^2
  phi[phi < 0] <- 0
  se <- sqrt(phi)
  return(se)
}

#' Plot the effect of the calibration
#'
#' @description
#' `plotCalibrationEffect` creates a plot showing the effect of the calibration.
#'
#' @details
#' Creates a plot with the effect estimate on the x-axis and the standard error on the y-axis.
#' Negative controls are shown as blue dots, positive controls as yellow diamonds. The area below the
#' dashed line indicated estimates with p < 0.05. The orange area indicates estimates with calibrated
#' p < 0.05.
#'
#' @param estimateNegatives  A numeric vector of effect estimates of the negative controls.
#' @param seNegatives        The standard error of the effect estimates of the negative controls.
#' @param estimatePositives  Optional: A numeric vector of effect estimates of the positive controls.
#' @param sePositives        Optional: The standard error of the effect estimates of the positive controls.
#' @template Scale
#' @param null               An object representing the fitted null distribution as created by the
#'                           `fitNull` or `fitMcmcNull` functions. If not provided, a null
#'                           will be fitted before plotting.
#' @param alpha              The alpha for the hypothesis test.
#' @param xLabel             The label on the x-axis: the name of the effect estimate.
#' @param title              Optional: the main title for the plot
#' @param showCis            Show 95 percent credible intervals for the calibrated p = alpha boundary.
#' @param showExpectedAbsoluteSystematicError  Show the expected absolute systematic error. If `null` is of
#'                           type `mcmcNull` the 95 percent credible interval will also be shown.
#' @param fileName           Name of the file where the plot should be saved, for example 'plot.png'.
#'                           See the function `ggsave` in the ggplot2 package for supported file
#'                           formats.
#' @param xLimits            Vector of length 2 for limits of the plot x axis - defaults to 0.25, 10
#' @param yLimits            Vector of length 2 for size limits of the y axis - defaults to 0, 1.5
#' @param logRrNegatives     DEPRECATED. Use `estimateNegatives` with `scale = "log relative risk"` instead.
#' @param seLogRrNegatives   DEPRECATED. Use `seNegatives` with `scale = "log relative risk"` instead.
#' @param logRrPositives     DEPRECATED. Use `estimatePositives` with `scale = "log relative risk"` instead.
#' @param seLogRrPositives   DEPRECATED. Use `sePositives` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' positive <- sccs[sccs$groundTruth == 1, ]
#' plotCalibrationEffect(negatives$logRr, negatives$seLogRr, positive$logRr, positive$seLogRr)
#'
#' @export
plotCalibrationEffect <- function(estimateNegatives,
                                  seNegatives,
                                  estimatePositives = NULL,
                                  sePositives = NULL,
                                  null = NULL,
                                  scale = "log relative risk",
                                  alpha = 0.05,
                                  xLabel = "Effect estimate",
                                  title = NULL,
                                  showCis = FALSE,
                                  showExpectedAbsoluteSystematicError = FALSE,
                                  fileName = NULL,
                                  xLimits = NULL,
                                  yLimits = NULL,
                                  logRrNegatives = NULL,
                                  seLogRrNegatives = NULL,
                                  logRrPositives = NULL,
                                  seLogRrPositives = NULL) {
  if (!is.null(logRrNegatives)) {
    warnDeprecatedArguments(
      c("logRrNegatives", "seLogRrNegatives", "logRrPositives", "seLogRrPositives"),
      c("estimateNegatives", "seNegatives", "estimatePositives", "sePositives")
    )
    estimateNegatives <- logRrNegatives
    seNegatives <- seLogRrNegatives
    estimatePositives <- logRrPositives
    sePositives <- seLogRrPositives
    scale <- LOG_RELATIVE_RISK
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimateNegatives, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(seNegatives, min.len = 1, lower = 0, add = errorMessages)
  checkmate::assertNumeric(estimatePositives, min.len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(sePositives, min.len = 1, lower = 0, null.ok = TRUE, add = errorMessages)
  checkmate::assertChoice(scale, choices = c(LOG_RELATIVE_RISK, RISK_DIFFERENCE), add = errorMessages)
  checkmate::assertClass(null, classes = c("null"), null.ok = TRUE, add = errorMessages)
  checkmate::assertNumber(alpha, lower = 0, upper = 1, add = errorMessages)
  checkmate::assertString(xLabel, add = errorMessages)
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertLogical(showCis, len = 1, add = errorMessages)
  checkmate::assertLogical(showExpectedAbsoluteSystematicError, len = 1, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(xLimits, len = 2, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(yLimits, len = 2, null.ok = TRUE, lower = 0, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(null)) {
    if (showCis) {
      null <- fitMcmcNull(estimateNegatives, seNegatives)
    } else {
      null <- fitNull(estimateNegatives, seNegatives)
    }
  }
  if (showCis && !is(null, "mcmcNull")) {
    stop("Cannot show credible intervals when using asymptotic null. Please use 'fitMcmcNull' to fit the null")
  }
  if (scale == LOG_RELATIVE_RISK) {
    if (is.null(xLimits)) {
      xLimits <- c(log(0.1), log(10))
    }
    if (is.null(yLimits)) {
      yLimits <- c(0, 1.5)
    }
    breakLabels <- c(0.1, 0.25, 0.5, 1, 2 * (1:ceiling(exp(xLimits[2]) / 2.0)))
    breaks <- log(breakLabels)
    breakIdx <- breaks >= xLimits[1] & breaks <= xLimits[2]
    breakLabels <- breakLabels[breakIdx]
    breaks <- breaks[breakIdx]
  } else {
    # scale == RISK_DIFFERENCE
    if (is.null(xLimits)) {
      maxAbsEstimate <- quantile(abs(estimateNegatives), 0.9, na.rm = TRUE)
      if (!is.null(estimatePositives)) {
        maxAbsEstimate <- max(c(maxAbsEstimate, abs(estimatePositives)), na.rm = TRUE)
      }
      xLimits <- c(-maxAbsEstimate * 1.25, maxAbsEstimate * 1.25)
    }
    if (is.null(yLimits)) {
      maxSe <- quantile(seNegatives, 0.9, na.rm = TRUE)
      if (!is.null(sePositives)) {
        maxSe <- max(c(maxSe, sePositives), na.rm = TRUE)
      }
      yLimits <- c(0, maxSe * 1.10)
    }
    # Let ggplot figure out the breaks:
    breaks <- ggplot2::waiver()
    breakLabels <- ggplot2::waiver() 
  }
  
  
  x <- seq(xLimits[1], xLimits[2], length.out = 400)
  if (is(null, "mcmcNull")) {
    chain <- attr(null, "mcmc")$chain
    matrix <- apply(chain, 1, function(null) logRrtoSE(x, alpha, null[1], 1 / sqrt(null[2])))
    ys <- apply(matrix, 1, function(se) quantile(se, c(0.025, 0.50, 0.975), na.rm = TRUE))
    rm(matrix)
    y <- ys[2, ]
    yLb <- ys[1, ]
    yUb <- ys[3, ]
  } else {
    y <- logRrtoSE(x, alpha, null[1], null[2])
  } 
  seTheoretical <- sapply(x, FUN = function(x) {
    abs(x) / qnorm(1 - alpha / 2)
  })
  vizData <- data.frame(x = x, y = y, seTheoretical = seTheoretical)
  
  checkWithinLimits(yLimits, c(seNegatives, sePositives), "yLimits")
  checkWithinLimits(xLimits, c(estimateNegatives, estimatePositives), "xLimits")
  
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  plot <- ggplot2::ggplot(vizData, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_vline(xintercept = 0, size = 1) +
    ggplot2::geom_area(
      fill = rgb(1, 0.5, 0, alpha = 0.5),
      color = rgb(1, 0.5, 0),
      size = 1,
      alpha = 0.5
    )
  
  if (showCis) {
    plot <- plot +
      ggplot2::geom_ribbon(ggplot2::aes(
        ymin = yLb,
        ymax = yUb
      ), fill = rgb(0.8, 0.2, 0.2), alpha = 0.3) +
      ggplot2::geom_line(ggplot2::aes(y = yLb),
                         colour = rgb(0.8, 0.2, 0.2, alpha = 0.2),
                         size = 1
      ) +
      ggplot2::geom_line(ggplot2::aes(y = yUb),
                         colour = rgb(0.8, 0.2, 0.2, alpha = 0.2),
                         size = 1
      )
  }
  if (showExpectedAbsoluteSystematicError) {
    ease <- computeExpectedAbsoluteSystematicError(null)
    if (!is(null, "mcmcNull")) {
      label <- sprintf("Expected absolute systematic error = %0.2f", ease)
    } else {
      label <- sprintf(
        "Expected absolute systematic error = %0.2f (%0.2f - %0.2f)",
        ease$ease,
        ease$ciLb,
        ease$ciUb
      )
    }
    dummy <- data.frame(text = label)
    if (scale == LOG_RELATIVE_RISK) {
      labelX <- log10(0.26)
      labelY <- 1.49
    } else {
      labelX <- xLimits[1] + 0.02 * (xLimits[2] - xLimits[1])
      labelY <- yLimits[2] * 0.99
    }
    plot <- plot + ggplot2::geom_label(x = labelX, y = labelY, hjust = "left", vjust = "top", alpha = 0.9, ggplot2::aes(label = .data$text), data = dummy, size = 3.5)
  }
  
  plot <- plot +
    ggplot2::geom_area(ggplot2::aes(y = .data$seTheoretical),
                       fill = rgb(0, 0, 0),
                       colour = rgb(0, 0, 0, alpha = 0.1),
                       alpha = 0.1
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$seTheoretical),
                       colour = rgb(0, 0, 0),
                       linetype = "dashed",
                       size = 1,
                       alpha = 0.5
    ) +
    ggplot2::geom_point(
      shape = 16,
      ggplot2::aes(x, y),
      data = data.frame(x = estimateNegatives, y = seNegatives),
      size = 2,
      alpha = 0.5,
      color = rgb(0, 0, 0.8)
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_x_continuous(xLabel,
                                # trans = "log10",
                                limits = xLimits,
                                breaks = breaks,
                                labels = breakLabels
    ) +
    ggplot2::scale_y_continuous("Standard Error") +
    ggplot2::coord_cartesian(ylim = yLimits) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "#AAAAAA"),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      legend.key = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  if (!is.null(estimatePositives)) {
    plot <- plot + ggplot2::geom_point(
      shape = 23,
      ggplot2::aes(x = .data$x, y = .data$y),
      data = data.frame(
        x = estimatePositives,
        y = sePositives
      ),
      size = 4,
      fill = rgb(1, 1, 0),
      alpha = 0.8
    )
  }
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }
  return(plot)
}

checkWithinLimits <- function(limits, values, label) {
  values <- values[!is.na(values)]
  if (length(values) > 0) {
    if (limits[1] > min(values) || limits[2] < max(values)) {
      warning(sprintf("Values are outside plotted range. Consider adjusting %s parameter", label))
    }
  }
}

#' Create a calibration plot
#'
#' @description
#' `plotCalibration` creates a plot showing the calibration of our calibration procedure
#'
#' @details
#' Creates a calibration plot showing the number of effects with p < alpha for every level of alpha.
#' The empirical calibration is performed using a leave-one-out design: The p-value of an effect is
#' computed by fitting a null using all other negative controls. Ideally, the calibration line should
#' approximate the diagonal. The plot shows both theoretical (traditional) and empirically calibrated
#' p-values.
#'
#' @param estimate         A numeric vector of effect estimates. On the log relative risk scale these
#'                         are log relative risks.
#' @param se               The standard error of the effect estimates. Hint: on the log relative risk
#'                         scale often the standard error = (log(<lower bound 95 percent confidence
#'                         interval>) - log(<effect estimate>))/qnorm(0.025)
#' @template Scale
#' @param useMcmc          Use MCMC to estimate the calibrated P-value?
#' @param legendPosition   Where should the legend be positioned? ("none", "left", "right", "bottom",
#'                         "top")
#' @param title            Optional: the main title for the plot
#' @param fileName         Name of the file where the plot should be saved, for example 'plot.png'. See
#'                         the function `ggsave` in the ggplot2 package for supported file
#'                         formats.
#' @param logRr            DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr          DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' plotCalibration(negatives$logRr, negatives$seLogRr)
#'
#' @export
plotCalibration <- function(estimate,
                            se,
                            scale = "log relative risk",
                            useMcmc = FALSE,
                            legendPosition = "right",
                            title = NULL,
                            fileName = NULL,
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
  checkmate::assertChoice(scale, choices = c(LOG_RELATIVE_RISK, RISK_DIFFERENCE), add = errorMessages)
  checkmate::assertLogical(useMcmc, len = 1, add = errorMessages)
  checkmate::assertString(legendPosition, add = errorMessages)
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (any(is.infinite(se))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting null distribution")
    estimate <- estimate[!is.infinite(se)]
    se <- se[!is.infinite(se)]
  }
  if (any(is.infinite(estimate))) {
    warning("Estimate(s) with infinite estimate detected. Removing before fitting null distribution")
    se <- se[!is.infinite(estimate)]
    estimate <- estimate[!is.infinite(estimate)]
  }
  if (any(is.na(se))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting null distribution")
    estimate <- estimate[!is.na(se)]
    se <- se[!is.na(se)]
  }
  if (any(is.na(estimate))) {
    warning("Estimate(s) with NA estimate detected. Removing before fitting null distribution")
    se <- se[!is.na(estimate)]
    estimate <- estimate[!is.na(estimate)]
  }
  
  data <- data.frame(estimate = estimate, SE = se)
  data$Z <- data$estimate / data$SE
  data$P <- 2 * pmin(pnorm(data$Z), 1 - pnorm(data$Z)) # 2-sided p-value
  data$Y <- sapply(data$P, function(x) {
    sum(data$P < x) / nrow(data)
  })
  
  data$calibratedP <- vector(length = nrow(data))
  for (i in 1:nrow(data)) {
    dataLeaveOneOut <- data[seq(1, nrow(data)) != i, ]
    if (useMcmc) {
      null <- fitMcmcNull(dataLeaveOneOut$estimate, dataLeaveOneOut$SE)
    } else {
      null <- fitNull(dataLeaveOneOut$estimate, dataLeaveOneOut$SE)
    }
    data$calibratedP[i] <- calibrateP(null, data$estimate[i], data$SE[i], pValueOnly = TRUE)
  }
  data <- data[!is.na(data$calibratedP), ]
  data$AdjustedY <- sapply(data$calibratedP, function(x) {
    sum(data$calibratedP < x) / nrow(data)
  })
  
  catData <- data.frame(
    x = c(data$P, data$calibratedP),
    y = c(data$Y, data$AdjustedY),
    label = factor(c(
      rep("Theoretical", times = nrow(data)),
      rep("Empirical", times = nrow(data))
    ))
  )
  catData$label <- factor(catData$label, levels = c("Empirical", "Theoretical"))
  
  names(catData) <- c("x", "y", "P-value calculation")
  
  breaks <- c(0, 0.25, 0.5, 0.75, 1)
  theme <- ggplot2::element_text(colour = "#000000", size = 10)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 10, hjust = 1)
  plot <- ggplot2::ggplot(
    catData,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      colour = .data$`P-value calculation`,
      linetype = .data$`P-value calculation`
    )
  ) +
    ggplot2::geom_vline(
      xintercept = breaks,
      colour = "#AAAAAA",
      lty = 1,
      size = 0.3
    ) +
    ggplot2::geom_vline(xintercept = 0.05, colour = "#888888", linetype = "dashed", size = 1) +
    ggplot2::geom_hline(yintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.3) +
    ggplot2::geom_abline(colour = "#AAAAAA", lty = 1, size = 0.3) +
    ggplot2::geom_step(direction = "hv", size = 1) +
    ggplot2::scale_colour_manual(values = c(rgb(0, 0, 0), rgb(0, 0, 0), rgb(0.5, 0.5, 0.5))) +
    ggplot2::scale_linetype_manual(values = c("solid", "twodash")) +
    ggplot2::scale_x_continuous(expression(alpha), limits = c(0, 1), breaks = c(breaks, 0.05), labels = c("", ".25", ".50", ".75", "1", ".05")) +
    ggplot2::scale_y_continuous(expression(paste("Fraction with p < ", alpha)), limits = c(0, 1), breaks = breaks, labels = c("0", ".25", ".50", ".75", "1")) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      strip.text.x = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = legendPosition
    )
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }
  return(plot)
}

#' Create a confidence interval calibration plot
#'
#' @description
#' `plotCalibration` creates a plot showing the calibration of our confidence interval
#' calibration procedure
#'
#' @details
#' Creates a calibration plot showing the fraction of effects within the confidence interval. The
#' empirical calibration is performed using a leave-one-out design: The confidence interval of an
#' effect is computed by fitting a null using all other controls. Ideally, the calibration line should
#' approximate the diagonal. The plot shows the coverage for both theoretical (traditional) and
#' empirically calibrated confidence intervals.
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
#' @param legacy                 If true, a legacy error model will be fitted, meaning standard
#'                               deviation is linear on the log scale. If false, standard deviation
#'                               is assumed to be simply linear.
#' @param evaluation             A data frame as generated by the [evaluateCiCalibration()]
#'                               function. If provided, the estimate, se, trueEffect, scale, strata, and
#'                               legacy arguments will be ignored.
#' @param crossValidationGroup   What should be the unit for the cross-validation? By default the unit
#'                               is a single control, but a different grouping can be provided, for
#'                               example linking a negative control to synthetic positive controls
#'                               derived from that negative control.
#' @param legendPosition         Where should the legend be positioned? ("none", "left", "right",
#'                               "bottom", "top").
#' @param title                  Optional: the main title for the plot
#' @param fileName               Name of the file where the plot should be saved, for example
#'                               'plot.png'. See the function `ggsave` in the ggplot2 package for
#'                               supported file formats.
#' @param logRr                  DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr                DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#' @param trueLogRr              DEPRECATED. Use `trueEffect` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#'
#' @examples
#' \dontrun{
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' plotCiCalibration(data$estimate, data$se, data$trueEffect)
#' }
#' @export
plotCiCalibration <- function(estimate,
                              se,
                              trueEffect,
                              scale = "log relative risk",
                              strata = as.factor(trueEffect),
                              crossValidationGroup = 1:length(estimate),
                              legacy = FALSE,
                              evaluation = NULL,
                              legendPosition = "top",
                              title = NULL,
                              fileName = NULL,
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
  checkmate::assertDataFrame(evaluation, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(legendPosition, add = errorMessages)
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(trueLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(evaluation)) {
    evaluation <- evaluateCiCalibration(
      estimate = estimate,
      se = se,
      trueEffect = trueEffect,
      strata = strata,
      crossValidationGroup = crossValidationGroup,
      legacy = legacy
    )
  }
  vizData <- evaluation[evaluation$label == "Within confidence interval", ]
  breaks <- c(0, 0.25, 0.5, 0.75, 1)
  theme <- ggplot2::element_text(colour = "#000000", size = 10)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 10, hjust = 1)
  plot <- ggplot2::ggplot(
    vizData,
    ggplot2::aes(
      x = .data$ciWidth,
      y = .data$coverage,
      colour = .data$`Confidence interval calculation`,
      linetype = .data$`Confidence interval calculation`
    )
  ) +
    ggplot2::geom_vline(
      xintercept = breaks,
      colour = "#AAAAAA",
      lty = 1,
      size = 0.3
    ) +
    ggplot2::geom_vline(xintercept = 0.95, colour = "#888888", linetype = "dashed", size = 1) +
    ggplot2::geom_hline(yintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.3) +
    ggplot2::geom_abline(colour = "#AAAAAA", lty = 1, size = 0.3) +
    ggplot2::geom_line(size = 1) +
    ggplot2::scale_colour_manual(values = c(rgb(0, 0, 0), rgb(0, 0, 0), rgb(0.5, 0.5, 0.5))) +
    ggplot2::scale_linetype_manual(values = c("solid", "twodash")) +
    ggplot2::scale_x_continuous("Width of CI", limits = c(0, 1), breaks = c(breaks, 0.95), labels = c("0", ".25", ".50", ".75", "", ".95")) +
    ggplot2::scale_y_continuous("Coverage", limits = c(0, 1), breaks = breaks, labels = c("0", ".25", ".50", ".75", "1")) +
    ggplot2::facet_grid(. ~ trueEffect) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = legendPosition
    )
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    width <- 1 + 2 * length(levels(evaluation$trueEffect))
    ggplot2::ggsave(fileName, plot, width = width, height = 3.5, dpi = 400)
  }
  return(plot)
}

#' Plot true and observed values
#'
#' @description
#' Plot true and observed values, for example from a simulation study.
#'
#' @details
#' Creates a forest plot of effect size estimates. Estimates that are significantly different
#' from the true value (alpha = 0.05) are marked in orange, others are marked in blue. For the log relative
#' risk scale the estimates are shown as ratios on a log axis, for the risk difference scale on a linear axis.
#'
#' @param estimate    A numeric vector of effect estimates. On the log relative risk scale these are
#'                    log relative risks.
#' @param se          The standard error of the effect estimates. Hint: on the log relative risk scale
#'                    often the standard error = (log(<lower bound 95 percent confidence interval>) -
#'                    log(<effect estimate>))/qnorm(0.025).
#' @param trueEffect  A vector of the true effect sizes (on the same scale as `estimate`).
#' @template Scale
#' @param xLabel      The label on the x-axis: the name of the effect estimate. Defaults to
#'                    "Relative risk" or "Risk difference", depending on the scale.
#' @param title       Optional: the main title for the plot
#' @param fileName    Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                    function `ggsave` in the ggplot2 package for supported file formats.
#' @param logRr       DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr     DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#' @param trueLogRr   DEPRECATED. Use `trueEffect` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' plotTrueAndObserved(data$estimate, data$se, data$trueEffect)
#'
#' @export
plotTrueAndObserved <- function(estimate,
                                se,
                                trueEffect,
                                scale = "log relative risk",
                                xLabel = NULL,
                                title = NULL,
                                fileName = NULL,
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
  checkmate::assertString(xLabel, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(trueLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(xLabel)) {
    xLabel <- if (scale == LOG_RELATIVE_RISK) "Relative risk" else "Risk difference"
  }
  theme <- ggplot2::element_text(colour = "#000000", size = 6)
  col <- c(rgb(0, 0, 0.8, alpha = 1), rgb(0.8, 0.4, 0, alpha = 1))
  colFill <- c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0.4, 0, alpha = 0.5))
  data <- data.frame(
    estimate = estimate,
    lb95 = estimate + qnorm(0.025) * se,
    ub95 = estimate + qnorm(0.975) * se,
    trueEffect = trueEffect
  )
  data$significant <- data$lb95 > data$trueEffect | data$ub95 < data$trueEffect
  data <- data[order(data$trueEffect, data$estimate), ]
  data$order <- 1:nrow(data)
  if (scale == LOG_RELATIVE_RISK) {
    breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
    data$trueEffectLabel <- round(exp(data$trueEffect), 2)
    data$estimate <- exp(data$estimate)
    data$lb95 <- exp(data$lb95)
    data$ub95 <- exp(data$ub95)
    xScale <- ggplot2::scale_x_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks)
    xLimits <- ggplot2::coord_cartesian(xlim = c(0.25, 10))
  } else {
    data$trueEffectLabel <- signif(data$trueEffect, 3)
    xScale <- ggplot2::scale_x_continuous(xLabel)
    xLimits <- ggplot2::coord_cartesian()
  }
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$estimate,
      y = .data$order,
      xmin = .data$lb95,
      xmax = .data$ub95,
      colour = .data$significant,
      fill = .data$significant
    )
  )
  if (scale == LOG_RELATIVE_RISK) {
    plot <- plot + ggplot2::geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2)
  }
  plot <- plot +
    ggplot2::geom_errorbarh(ggplot2::aes(xmax = .data$trueEffectLabel, xmin = .data$trueEffectLabel), height = 1, color = rgb(0, 0, 0), size = 1) +
    ggplot2::geom_errorbarh(height = 0) +
    ggplot2::geom_point(shape = 21, size = 1.5) +
    ggplot2::scale_colour_manual(values = col) +
    ggplot2::scale_fill_manual(values = colFill) +
    xLimits +
    xScale +
    ggplot2::facet_grid(trueEffectLabel ~ ., scales = "free_y", space = "free") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = "#EEEEEE"),
      axis.ticks = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = theme, legend.key = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 5, height = 7, dpi = 400)
  }
  return(plot)
}

#' Plot the MCMC trace
#'
#' @details
#' Plot the trace of the MCMC for diagnostics purposes.
#'
#' @param mcmcNull   An object of type `mcmcNull` as generated using the `fitMcmcNull`
#'                   function.
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                   function `ggsave` in the ggplot2 package for supported file formats.
#'
#' @examples
#' \dontrun{
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitMcmcNull(negatives$logRr, negatives$seLogRr)
#' plotMcmcTrace(null)
#' }
#' @export
plotMcmcTrace <- function(mcmcNull, fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(mcmcNull, classes = "mcmcNull", add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  mcmc <- attr(mcmcNull, "mcmc")
  dataMean <- data.frame(x = 1:nrow(mcmc$chain), trace = as.numeric(ts(mcmc$chain[
    ,
    1
  ])), var = "Mean")
  dataPrecision <- data.frame(x = 1:nrow(mcmc$chain), trace = as.numeric(ts(mcmc$chain[
    ,
    2
  ])), var = "Precision")
  data <- rbind(dataMean, dataPrecision)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$trace)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::scale_x_continuous("Iterations") +
    ggplot2::facet_grid(var ~ ., scales = "free") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 5, height = 3.5, dpi = 400)
  }
  return(plot)
}


#' Create a confidence interval coverage plot
#'
#' @description
#' `plotCiCoverage` creates a plot showing the coverage before and after confidence interval
#' calibration at various widths of the confidence interval.
#'
#' @details
#' Creates a plot showing the fraction of effects above, within, and below the confidence interval. The
#' empirical calibration is performed using a leave-one-out design: The confidence interval of an
#' effect is computed by fitting a null using all other controls. The plot shows the coverage for
#' both theoretical (traditional) and empirically calibrated confidence intervals.
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
#' @param legacy                 If true, a legacy error model will be fitted, meaning standard
#'                               deviation is linear on the log scale. If false, standard deviation
#'                               is assumed to be simply linear.
#' @param evaluation             A data frame as generated by the [evaluateCiCalibration()]
#'                               function. If provided, the estimate, se, trueEffect, scale, strata, and
#'                               legacy arguments will be ignored.
#' @param crossValidationGroup   What should be the unit for the cross-validation? By default the unit
#'                               is a single control, but a different grouping can be provided, for
#'                               example linking a negative control to synthetic positive controls
#'                               derived from that negative control.
#' @param legendPosition         Where should the legend be positioned? ("none", "left", "right",
#'                               "bottom", "top").
#' @param title                  Optional: the main title for the plot
#' @param fileName               Name of the file where the plot should be saved, for example
#'                               'plot.png'. See the function `ggsave` in the ggplot2 package for
#'                               supported file formats.
#' @param logRr                  DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr                DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#' @param trueLogRr              DEPRECATED. Use `trueEffect` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#'
#' @examples
#' \dontrun{
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' plotCiCoverage(data$estimate, data$se, data$trueEffect)
#' }
#' @export
plotCiCoverage <- function(estimate,
                           se,
                           trueEffect,
                           scale = "log relative risk",
                           strata = as.factor(trueEffect),
                           crossValidationGroup = 1:length(estimate),
                           legacy = FALSE,
                           evaluation = NULL,
                           legendPosition = "top",
                           title = NULL,
                           fileName = NULL,
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
  checkmate::assertDataFrame(evaluation, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(legendPosition, add = errorMessages)
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(trueLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(evaluation)) {
    evaluation <- evaluateCiCalibration(
      estimate = estimate,
      se = se,
      trueEffect = trueEffect,
      strata = strata,
      crossValidationGroup = crossValidationGroup,
      legacy = legacy
    )
  }
  evaluation$label <- factor(evaluation$label, levels = rev(c("Below confidence interval", "Within confidence interval", "Above confidence interval")))
  evaluation$`Confidence interval calculation` <- factor(evaluation$`Confidence interval calculation`, levels = rev(c("Calibrated", "Uncalibrated")))
  evaluation$trueEffect <- paste0("True effect size = ", evaluation$trueEffect)
  breaks <- c(0, 0.25, 0.5, 0.75, 1)
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  plot <- ggplot2::ggplot(
    evaluation,
    ggplot2::aes(
      x = .data$ciWidth,
      y = .data$coverage,
      colour = .data$label,
      fill = .data$label,
      group = .data$label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = breaks,
      colour = "#AAAAAA",
      lty = 1,
      size = 0.3
    ) +
    ggplot2::geom_hline(yintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.3) +
    ggplot2::geom_area(position = "stack", alpha = 0.6) +
    ggplot2::geom_abline(colour = "#000000", size = 1, intercept = 0.5, slope = 0.5, linetype = "dashed") +
    ggplot2::geom_abline(colour = "#000000", size = 1, intercept = 0.5, slope = -0.5, linetype = "dashed") +
    ggplot2::scale_colour_manual(values = c(rgb(0.8, 0, 0), rgb(0, 0, 0), rgb(0, 0, 0.8))) +
    ggplot2::scale_fill_manual(values = c(rgb(0.8, 0, 0), rgb(0.75, 0.75, 0.75), rgb(0, 0, 0.8))) +
    ggplot2::scale_x_continuous("Width of the confidence interval", limits = c(0, 1), breaks = breaks, labels = c("0", ".25", ".50", ".75", "1")) +
    ggplot2::scale_y_continuous("Fraction", limits = c(0, 1), breaks = breaks, labels = c("0", ".25", ".50", ".75", "1")) +
    ggplot2::facet_grid(`Confidence interval calculation` ~ trueEffect) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = theme,
      strip.text.y = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = legendPosition,
      legend.title = ggplot2::element_blank(),
      legend.text = theme
    )
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    width <- 1 + 1.8 * length(levels(as.factor(evaluation$trueEffect)))
    ggplot2::ggsave(fileName, plot, width = width, height = 5, dpi = 400)
  }
  return(plot)
}

#' Plot the systematic error model
#'
#' @description
#' `plotErrorModel` creates a plot showing the systematic error model.
#'
#' @details
#' Creates a plot with the true effect size on the x-axis, and the mean plus and minus the standard
#' deviation shown on the y-axis. Also shown are simple error models fitted at each true effect
#' size in the input. For the log relative risk scale both axes show relative risks on a log scale, for the
#' risk difference scale both axes are linear.
#'
#' @param estimate               A numeric vector of effect estimates. On the log relative risk scale
#'                               these are log relative risks.
#' @param se                     The standard error of the effect estimates. Hint: on the log relative
#'                               risk scale often the standard error = (log(<lower bound 95 percent
#'                               confidence interval>) - log(<effect estimate>))/qnorm(0.025).
#' @param trueEffect             The true effect (on the same scale as `estimate`).
#' @template Scale
#' @param legacy                 If true, a legacy error model will be fitted, meaning standard
#'                               deviation is linear on the log scale. If false, standard deviation
#'                               is assumed to be simply linear. The legacy model is not supported for the
#'                               risk difference scale.
#' @param title              Optional: the main title for the plot
#' @param fileName           Name of the file where the plot should be saved, for example 'plot.png'.
#'                           See the function `ggsave` in the ggplot2 package for supported file
#'                           formats.
#' @param logRr              DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr            DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#' @param trueLogRr          DEPRECATED. Use `trueEffect` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' plotErrorModel(data$estimate, data$se, data$trueEffect)
#' @export
plotErrorModel <- function(estimate,
                           se,
                           trueEffect,
                           scale = "log relative risk",
                           title = NULL,
                           legacy = FALSE,
                           fileName = NULL,
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
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertLogical(legacy, len = 1, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(trueLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  model <- fitSystematicErrorModel(
    estimate = estimate,
    se = se,
    trueEffect = trueEffect,
    legacy = legacy
  )
  fitSimpleModel <- function(selectedTrueEffect) {
    idx <- trueEffect == selectedTrueEffect
    return(fitNull(estimate[idx], se[idx]))
  }
  simpleModels <- as.data.frame(t(sapply(unique(trueEffect), fitSimpleModel)))
  if (scale == LOG_RELATIVE_RISK) {
    simpleModels$trueEffect <- exp(unique(trueEffect))
    simpleModels$y <- exp(simpleModels$mean)
    simpleModels$ymin <- exp(simpleModels$mean - simpleModels$sd)
    simpleModels$ymax <- exp(simpleModels$mean + simpleModels$sd)
    
    breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
    x <- exp(seq(log(0.25), log(10), by = 0.01))
    y <- exp(model[1] + model[2] * log(x))
    if (legacy) {
      ymin <- exp(log(y) - exp(model[3] + model[4] * log(x)))
      ymax <- exp(log(y) + exp(model[3] + model[4] * log(x)))
    } else {
      ymin <- exp(log(y) - (model[3] + model[4] * abs(log(x))))
      ymax <- exp(log(y) + (model[3] + model[4] * abs(log(x))))
    }
    xScale <- ggplot2::scale_x_continuous("True effect size", trans = "log10", breaks = breaks, labels = breaks)
    yScale <- ggplot2::scale_y_continuous("Systematic error mean (plus and minus one SD)",
                                          trans = "log10",
                                          breaks = breaks,
                                          labels = breaks
    )
    limits <- ggplot2::coord_cartesian(xlim = c(0.25, 10), ylim = c(0.25, 10))
    nullValue <- 1
    errorBarWidth <- 0.1
  } else {
    simpleModels$trueEffect <- unique(trueEffect)
    simpleModels$y <- simpleModels$mean
    simpleModels$ymin <- simpleModels$mean - simpleModels$sd
    simpleModels$ymax <- simpleModels$mean + simpleModels$sd
    
    maxAbsEffect <- max(abs(trueEffect), na.rm = TRUE)
    if (maxAbsEffect == 0) {
      maxAbsEffect <- max(abs(estimate), na.rm = TRUE)
    }
    x <- seq(-maxAbsEffect * 1.25, maxAbsEffect * 1.25, length.out = 200)
    y <- model[1] + model[2] * x
    ymin <- y - (model[3] + model[4] * abs(x))
    ymax <- y + (model[3] + model[4] * abs(x))
    xScale <- ggplot2::scale_x_continuous("True effect size")
    yScale <- ggplot2::scale_y_continuous("Systematic error mean (plus and minus one SD)")
    limits <- ggplot2::coord_cartesian()
    nullValue <- 0
    errorBarWidth <- 0.05 * (max(x) - min(x))
  }
  data <- data.frame(x = x, y = y, ymin = ymin, ymax = ymax)
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      ymin = .data$ymin,
      ymax = .data$ymax
    )
  )
  if (scale == LOG_RELATIVE_RISK) {
    plot <- plot +
      ggplot2::geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
      ggplot2::geom_hline(yintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2)
  }
  plot <- plot +
    ggplot2::geom_vline(xintercept = nullValue, size = 1) +
    ggplot2::geom_hline(yintercept = nullValue, size = 1) +
    ggplot2::geom_ribbon(fill = rgb(0, 0, 0.8), alpha = 0.3) +
    ggplot2::geom_line(color = rgb(0, 0, 0.8)) +
    ggplot2::geom_errorbar(ggplot2::aes(x = .data$trueEffect), width = errorBarWidth, size = 1, color = rgb(0, 0, 0.8), data = simpleModels) +
    xScale +
    yScale +
    limits +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }
  return(plot)
}

plotIsobars <- function(null, alpha, xLabel = "Relative risk", sePositives = NULL, scale = LOG_RELATIVE_RISK, maxSe = 1) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(null, classes = "null", add = errorMessages)
  checkmate::assertNumber(alpha, add = errorMessages)
  checkmate::assertString(xLabel, add = errorMessages)
  checkmate::assertNumeric(sePositives, min.len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::assertChoice(scale, choices = c(LOG_RELATIVE_RISK, RISK_DIFFERENCE), add = errorMessages)
  checkmate::assertNumber(maxSe, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is(null, "mcmcNull")) {
    null <- c(null[1], 1 / sqrt(null[2]))
  }
  if (scale == LOG_RELATIVE_RISK) {
    x <- exp(seq(log(0.25), log(10), by = 0.01))
    y <- sapply(x, FUN = function(x) abs(log(x)) / qnorm(1 - alpha / 2))
    xOnEstimateScale <- log(x)
    nullOnPlotScale <- exp(null[1])
    breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
    breakLabels <- breaks
    xTrans <- "log10"
    xLimits <- c(0.25, 10)
    nullValue <- 1
  } else {
    xRange <- max(abs(null[1]) + 4 * null[2], maxSe * qnorm(1 - alpha / 2) * 1.25)
    x <- seq(-xRange, xRange, length.out = 400)
    y <- abs(x) / qnorm(1 - alpha / 2)
    xOnEstimateScale <- x
    nullOnPlotScale <- null[1]
    breaks <- ggplot2::waiver()
    breakLabels <- ggplot2::waiver()
    xTrans <- "identity"
    xLimits <- c(-xRange, xRange)
    nullValue <- 0
  }
  thresholds <- c(0.05, 0.25, 0.5, 0.75)
  isoBars <- lapply(thresholds, function(threshold) {
    data.frame(
      x = x,
      y = logRrtoSE(xOnEstimateScale, threshold, null[1], null[2]),
      ymax = y,
      threshold = threshold
    )
  })
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  plot <- ggplot2::ggplot(
    data.frame(x, y),
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_vline(xintercept = if (scale == LOG_RELATIVE_RISK) breaks else 0, colour = rgb(0, 0, 0), lty = 1, size = 0.5, alpha = 0.2) +
    ggplot2::geom_hline(yintercept = maxSe * 0:4 / 4, colour = rgb(0, 0, 0), lty = 1, size = 0.5, alpha = 0.2) +
    ggplot2::geom_vline(xintercept = nullValue, size = 1) +
    ggplot2::geom_line(
      colour = rgb(0, 0, 0),
      linetype = "dashed",
      size = 1,
      alpha = 0.5
    ) +
    ggplot2::scale_x_continuous(xLabel,
                                trans = xTrans,
                                breaks = breaks,
                                labels = breakLabels
    ) +
    ggplot2::scale_y_continuous("Standard Error") +
    ggplot2::coord_cartesian(xlim = xLimits, ylim = c(0, maxSe)) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      legend.key = ggplot2::element_blank(),
      strip.text.x = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  for (isoBar in isoBars) {
    isoBarLeft <- isoBar[isoBar$x < nullOnPlotScale, ]
    isoBarRight <- isoBar[isoBar$x > nullOnPlotScale, ]
    labelData <- data.frame(
      x = c(
        isoBarLeft$x[which.min(abs(isoBarLeft$y - 0.625 * maxSe))],
        isoBarRight$x[which.min(abs(isoBarRight$y - 0.625 * maxSe))]
      ),
      y = 0.625 * maxSe,
      label = paste0(100 * isoBar$threshold[1], "%")
    )
    plot <- plot + ggplot2::geom_ribbon(ggplot2::aes(
      ymin = .data$y,
      ymax = .data$ymax
    ),
    fill = rgb(1, 0.5, 0, alpha = 0.2),
    data = isoBar[isoBar$y < isoBar$ymax, ]
    ) +
      ggplot2::geom_line(
        color = rgb(1, 0.5, 0),
        size = 1,
        alpha = 0.5,
        data = isoBar
      ) +
      ggplot2::geom_label(ggplot2::aes(label = .data$label),
                          color = rgb(1, 0.5, 0),
                          data = labelData
      ) +
      ggplot2::geom_text(ggplot2::aes(label = .data$label),
                         data = labelData
      )
  }
  if (!is.null(sePositives)) {
    plot <- plot + ggplot2::geom_hline(yintercept = sePositives)
  }
  return(plot)
}

#' Plot the expected type 1 error as a function of standard error
#'
#' @description
#' `plotExpectedType1Error` creates a plot showing the expected type 1 error as a function of standard error.
#'
#' @details
#' Creates a plot with the standard error on the x-axis and the expected type 1 error on the y-axis. The
#' red line indicates the expected type 1 error  given the estimated empirical null distribution if no
#' calibration is performed. The dashed line indicated the nominal expected type 1 error rate, assuming
#' the theoretical null distribution.
#'
#' If standard errors are provided for non-negative estimates these will be plotted on the red line as
#' yellow diamonds.
#'
#' @param estimateNegatives  A numeric vector of effect estimates of the negative controls. On the log
#'                           relative risk scale these are log relative risks.
#' @param seNegatives        The standard error of the effect estimates of the negative controls.
#' @param sePositives        The standard error of the effect estimates of the positive controls.
#' @template Scale
#' @param alpha              The alpha (nominal type 1 error) to be used.
#' @param null               An object representing the fitted null distribution as created by the
#'                           `fitNull` function. If not provided, a null will be fitted before
#'                           plotting.
#' @param xLabel             If showing effect sizes, what label should be used for the effect size axis?
#'                           Defaults to "Relative risk" or "Risk difference", depending on the scale.
#' @param title              Optional: the main title for the plot
#' @param showCis            Show 95 percent credible intervals for the expected type 1 error.
#' @param showEffectSizes    Show the expected effect sizes alongside the expected type 1 error?
#' @param fileName           Name of the file where the plot should be saved, for example 'plot.png'.
#'                           See the function `ggsave` in the ggplot2 package for supported file
#'                           formats.
#' @param logRrNegatives     DEPRECATED. Use `estimateNegatives` with `scale = "log relative risk"` instead.
#' @param seLogRrNegatives   DEPRECATED. Use `seNegatives` with `scale = "log relative risk"` instead.
#' @param seLogRrPositives   DEPRECATED. Use `sePositives` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' positive <- sccs[sccs$groundTruth == 1, ]
#' plotExpectedType1Error(negatives$logRr, negatives$seLogRr, positive$seLogRr)
#'
#' @export
plotExpectedType1Error <- function(estimateNegatives,
                                   seNegatives,
                                   sePositives = NULL,
                                   scale = "log relative risk",
                                   alpha = 0.05,
                                   null = NULL,
                                   xLabel = NULL,
                                   title = NULL,
                                   showCis = FALSE,
                                   showEffectSizes = FALSE,
                                   fileName = NULL,
                                   logRrNegatives = NULL,
                                   seLogRrNegatives = NULL,
                                   seLogRrPositives = NULL) {
  if (!is.null(logRrNegatives)) {
    warnDeprecatedArguments(
      c("logRrNegatives", "seLogRrNegatives", "seLogRrPositives"),
      c("estimateNegatives", "seNegatives", "sePositives")
    )
    estimateNegatives <- logRrNegatives
    seNegatives <- seLogRrNegatives
    if (!is.null(seLogRrPositives)) {
      sePositives <- seLogRrPositives
    }
    scale <- LOG_RELATIVE_RISK
  }
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(estimateNegatives, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(seNegatives, min.len = 1, add = errorMessages)
  checkmate::assertNumeric(sePositives, min.len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::assertChoice(scale, choices = c(LOG_RELATIVE_RISK, RISK_DIFFERENCE), add = errorMessages)
  checkmate::assertNumber(alpha, add = errorMessages)
  checkmate::assertClass(null, classes = "null", null.ok = TRUE, add = errorMessages)
  checkmate::assertString(xLabel, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertLogical(showCis, len = 1, add = errorMessages)
  checkmate::assertLogical(showEffectSizes, len = 1, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRrNegatives, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRrNegatives, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRrPositives, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(xLabel)) {
    xLabel <- if (scale == LOG_RELATIVE_RISK) "Relative risk" else "Risk difference"
  }
  if (is.null(null)) {
    if (showCis) {
      null <- fitMcmcNull(estimateNegatives, seNegatives)
    } else {
      null <- fitNull(estimateNegatives, seNegatives)
    }
  }
  if (showCis && !is(null, "mcmcNull")) {
    stop("Cannot show credible intervals when using asymptotic null. Please use 'fitMcmcNull' to fit the null")
  }
  
  if (scale == LOG_RELATIVE_RISK) {
    maxSe <- 1
  } else {
    maxSe <- max(c(quantile(seNegatives, 0.9, na.rm = TRUE), if (!is.null(sePositives)) sePositives), na.rm = TRUE) * 1.25
  }
  se <- maxSe * (0:100) / 100
  if (is(null, "mcmcNull")) {
    computeExpected <- function(se, chain) {
      mean <- chain[, 1]
      sd <- sqrt(se^2 + (1 / sqrt(chain[, 2]))^2)
      type1Error <- pnorm(qnorm(1 - alpha / 2, 0, se), mean, sd, lower.tail = FALSE) +
        pnorm(qnorm(alpha / 2, 0, se), mean, sd, lower.tail = TRUE)
      return(quantile(type1Error, c(0.025, 0.5, 0.975)))
    }
    mcmc <- attr(null, "mcmc")
    estimates <- sapply(se, computeExpected, chain = mcmc$chain)
    type1Error <- estimates[2, ]
    lb <- estimates[1, ]
    ub <- estimates[3, ]
  } else {
    mean <- null[1]
    sd <- sqrt(se^2 + null[2]^2)
    type1Error <- pnorm(qnorm(1 - alpha / 2, 0, se), mean, sd, lower.tail = FALSE) +
      pnorm(qnorm(alpha / 2, 0, se), mean, sd, lower.tail = TRUE)
  } 
  breaks <- 0:4 / 4
  seBreaks <- maxSe * breaks
  seLabels <- if (scale == LOG_RELATIVE_RISK) breaks else signif(seBreaks, 3)
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  plot <- ggplot2::ggplot(
    data.frame(se, type1Error),
    ggplot2::aes(x = .data$se, y = .data$type1Error)
  ) +
    ggplot2::geom_vline(xintercept = seBreaks, colour = rgb(0, 0, 0), lty = 1, size = 0.5, alpha = 0.2) +
    ggplot2::geom_hline(yintercept = breaks, colour = rgb(0, 0, 0), lty = 1, size = 0.5, alpha = 0.2) +
    ggplot2::geom_hline(
      yintercept = alpha,
      colour = rgb(0, 0, 0),
      linetype = "dashed",
      size = 1,
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      color = rgb(0.8, 0, 0),
      size = 1,
      alpha = 0.5
    )
  
  if (showCis) {
    plot <- plot +
      ggplot2::geom_ribbon(ggplot2::aes(
        ymin = lb,
        ymax = ub
      ), fill = rgb(0.8, 0.2, 0.2), alpha = 0.3) +
      ggplot2::geom_line(ggplot2::aes(y = lb),
                         colour = rgb(0.8, 0.2, 0.2, alpha = 0.2),
                         size = 1
      ) +
      ggplot2::geom_line(ggplot2::aes(y = ub),
                         colour = rgb(0.8, 0.2, 0.2, alpha = 0.2),
                         size = 1
      )
  }
  plot <- plot + ggplot2::geom_label(
    label = paste("alpha ==", alpha),
    data = data.frame(
      se = maxSe * (0.875 + showEffectSizes * 0.125),
      type1Error = alpha + showEffectSizes * 0.04
    ),
    parse = TRUE
  ) +
    ggplot2::scale_x_continuous("Standard error",
                                limits = c(0, maxSe),
                                breaks = seBreaks,
                                labels = seLabels
    ) +
    ggplot2::scale_y_continuous("Expected type 1 error",
                                limits = c(0, 1),
                                breaks = breaks,
                                labels = breaks
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      legend.key = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  if (!is.null(sePositives)) {
    yPos <- sapply(sePositives, FUN = function(x) {
      type1Error[which.min(abs(x - se))]
    })
    posData <- data.frame(
      x = sePositives,
      y = yPos
    )
    plot <- plot + ggplot2::geom_vline(ggplot2::aes(xintercept = .data$x),
                                       data = posData
    ) +
      ggplot2::geom_point(
        shape = 23,
        ggplot2::aes(x = .data$x, y = .data$y),
        data = posData,
        size = 4,
        fill = rgb(1, 1, 0),
        alpha = 0.8
      )
  }
  if (showEffectSizes) {
    plot <- plot + ggplot2::coord_flip()
    plot <- plot + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
    plot2 <- if (is.null(sePositives)) {
      plotIsobars(null = null, alpha = alpha, xLabel = xLabel, scale = scale, maxSe = maxSe)
    } else {
      plotIsobars(null = null, alpha = alpha, xLabel = xLabel, sePositives = sePositives, scale = scale, maxSe = maxSe)
    }
    plots <- list(plot2, plot)
    grobs <- heights <- list()
    for (i in 1:length(plots)) {
      grobs[[i]] <- ggplot2::ggplotGrob(plots[[i]])
      heights[[i]] <- grobs[[i]]$heights[6:9]
    }
    maxHeight <- do.call(grid::unit.pmax, heights)
    for (i in 1:length(grobs)) {
      grobs[[i]]$heights[6:9] <- as.list(maxHeight)
    }
    plot <- gridExtra::grid.arrange(grobs[[1]], grobs[[2]], nrow = 1, ncol = 2, widths = c(2, 1), top = title)
    if (!is.null(fileName)) {
      ggplot2::ggsave(fileName, plot, width = 10, height = 5, dpi = 400)
    }
  } else {
    if (!is.null(title)) {
      plot <- plot + ggplot2::ggtitle(title)
    }
    if (!is.null(fileName)) {
      ggplot2::ggsave(fileName, plot, width = 5, height = 5, dpi = 400)
    }
  }
  return(plot)
}


#' Plot the effect of the CI calibration
#'
#' @description
#' Creates a plot with the effect estimate on the x-axis and the standard error on the y-axis. The plot
#' is trellised by true effect size. Negative and positive controls are shown as blue dots. The area below the
#' dashed line indicated estimates that are statistically significant different from the true effect size (p < 0.05).
#' The orange area indicates estimates with calibrated p < 0.05.
#'
#' @param estimate    A numeric vector of effect estimates. On the log relative risk scale these are
#'                    log relative risks.
#' @param se          The standard error of the effect estimates. Hint: on the log relative risk scale
#'                    often the standard error = (log(<lower bound 95 percent confidence interval>) -
#'                    log(<effect estimate>))/qnorm(0.025).
#' @param trueEffect  A vector of the true effect sizes (on the same scale as `estimate`).
#' @template Scale
#' @param legacy      If true, a legacy error model will be fitted, meaning standard
#'                    deviation is linear on the log scale. If false, standard deviation
#'                    is assumed to be simply linear. The legacy model is not supported for the
#'                    risk difference scale.
#' @param model       The fitted systematic error model. If not provided, it will be fitted on the
#'                    provided data.
#' @param xLabel      The label on the x-axis: the name of the effect estimate. Defaults to
#'                    "Relative risk" or "Risk difference", depending on the scale.
#' @param title       Optional: the main title for the plot
#' @param fileName    Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                    function `ggsave` in the ggplot2 package for supported file formats.
#' @param logRr       DEPRECATED. Use `estimate` with `scale = "log relative risk"` instead.
#' @param seLogRr     DEPRECATED. Use `se` with `scale = "log relative risk"` instead.
#' @param trueLogRr   DEPRECATED. Use `trueEffect` with `scale = "log relative risk"` instead.
#'
#' @return
#' A Ggplot object. Use the `ggsave` function to save to file.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueEffect = log(c(1, 2, 4)))
#' plotCiCalibrationEffect(data$estimate, data$se, data$trueEffect)
#'
#' @export
plotCiCalibrationEffect <- function(estimate,
                                    se,
                                    trueEffect,
                                    scale = "log relative risk",
                                    legacy = FALSE,
                                    model = NULL,
                                    xLabel = NULL,
                                    title = NULL,
                                    fileName = NULL,
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
  checkmate::assertLogical(legacy, len = 1, add = errorMessages)
  checkmate::assertClass(model, classes = "systematicErrorModel", null.ok = TRUE, add = errorMessages)
  checkmate::assertString(xLabel, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(title, null.ok = TRUE, add = errorMessages)
  checkmate::assertString(fileName, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(logRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(seLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(trueLogRr, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(xLabel)) {
    xLabel <- if (scale == LOG_RELATIVE_RISK) "Relative risk" else "Risk difference"
  }
  alpha <- 0.05
  if (is.null(model)) {
    model <- fitSystematicErrorModel(
      estimate = estimate,
      se = se,
      trueEffect = trueEffect,
      estimateCovarianceMatrix = FALSE,
      legacy = legacy
    )
  } else {
    legacy <- (names(model)[3] == "logSdIntercept")
  }
  d <- data.frame(
    estimate = estimate,
    se = se,
    trueEffect = trueEffect,
    ci95lb = estimate + qnorm(0.025) * se,
    ci95ub = estimate + qnorm(0.975) * se
  )
  # Label for the true effect on the scale used for display:
  if (scale == LOG_RELATIVE_RISK) {
    d$trueEffectLabel <- exp(trueEffect)
  } else {
    d$trueEffectLabel <- signif(trueEffect, 3)
  }
  d <- d[!is.na(d$estimate), ]
  d <- d[!is.na(d$se), ]
  if (nrow(d) == 0) {
    return(NULL)
  }
  d$Group <- as.factor(d$trueEffectLabel)
  d$Significant <- d$ci95lb > d$trueEffect | d$ci95ub < d$trueEffect
  
  temp1 <- aggregate(Significant ~ trueEffect + trueEffectLabel, data = d, length)
  temp2 <- aggregate(Significant ~ trueEffect + trueEffectLabel, data = d, mean)
  temp1$nLabel <- paste0(formatC(temp1$Significant, big.mark = ","), " estimates")
  temp1$Significant <- NULL
  temp2$meanLabel <- paste0(
    formatC(100 * (1 - temp2$Significant), digits = 1, format = "f"),
    "% of CIs includes ",
    temp2$trueEffectLabel
  )
  temp2$Significant <- NULL
  dd <- merge(temp1, temp2)
  
  theme <- ggplot2::element_text(colour = "#000000", size = 10)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 10, hjust = 1)
  
  d$Group <- paste("True", tolower(xLabel), "=", d$trueEffectLabel)
  dd$Group <- paste("True", tolower(xLabel), "=", dd$trueEffectLabel)
  
  if (scale == LOG_RELATIVE_RISK) {
    breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
    xLimits <- log(c(0.1, 10))
    xBreaks <- log(breaks)
    xLabels <- breaks
    maxSe <- 1
  } else {
    maxAbs <- max(abs(c(d$estimate, d$trueEffect)), na.rm = TRUE) * 1.25
    xLimits <- c(-maxAbs, maxAbs)
    xBreaks <- ggplot2::waiver()
    xLabels <- ggplot2::waiver()
    maxSe <- max(d$se, na.rm = TRUE) * 1.1
  }
  
  x <- seq(xLimits[1], xLimits[2], length.out = 500)
  calBounds <- data.frame()
  for (i in 1:nrow(dd)) {
    mu <- model[1] + model[2] * dd$trueEffect[i]
    if (legacy) {
      sigma <- exp(model[3] + model[4] * dd$trueEffect[i])
    } else {
      sigma <- model[3] + model[4] * abs(dd$trueEffect[i])
    }
    calBounds <- rbind(
      calBounds,
      data.frame(
        estimate = x,
        se = logRrtoSE(x, alpha, mu, sigma),
        Group = dd$Group[i]
      )
    )
  }
  plot <- ggplot2::ggplot(d, ggplot2::aes(x = .data$estimate, y = .data$se))
  if (scale == LOG_RELATIVE_RISK) {
    plot <- plot + ggplot2::geom_vline(xintercept = xBreaks, colour = "#AAAAAA", lty = 1, size = 0.5)
  }
  plot <- plot +
    ggplot2::geom_area(
      fill = rgb(1, 0.5, 0, alpha = 0.5),
      color = rgb(1, 0.5, 0),
      size = 1,
      alpha = 0.5, data = calBounds
    ) +
    ggplot2::geom_abline(ggplot2::aes(intercept = (-.data$trueEffect) / qnorm(0.025), slope = 1 / qnorm(0.025)), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5, data = dd) +
    ggplot2::geom_abline(ggplot2::aes(intercept = (-.data$trueEffect) / qnorm(0.975), slope = 1 / qnorm(0.975)), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5, data = dd) +
    ggplot2::geom_point(
      shape = 16,
      size = 2,
      alpha = 0.5,
      color = rgb(0, 0, 0.8)
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_label(x = xLimits[1] + 0.025 * (xLimits[2] - xLimits[1]), y = 0.95 * maxSe, alpha = 1, hjust = "left", ggplot2::aes(label = .data$nLabel), size = 3.5, data = dd) +
    ggplot2::geom_label(x = xLimits[1] + 0.025 * (xLimits[2] - xLimits[1]), y = 0.8 * maxSe, alpha = 1, hjust = "left", ggplot2::aes(label = .data$meanLabel), size = 3.5, data = dd) +
    ggplot2::scale_x_continuous(xLabel, limits = xLimits, breaks = xBreaks, labels = xLabels) +
    ggplot2::scale_y_continuous("Standard Error") +
    ggplot2::coord_cartesian(ylim = c(0, maxSe)) +
    ggplot2::facet_grid(. ~ Group) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      axis.title = theme,
      legend.key = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      strip.text.x = theme,
      strip.text.y = theme,
      strip.background = ggplot2::element_blank(),
      legend.position = "none"
    )
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 2 + 3.5 * nrow(dd), height = 2.8, dpi = 400)
  }
  return(plot)
}
