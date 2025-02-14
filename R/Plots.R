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
#' \code{plotForest} creates a forest plot of effect size estimates.
#'
#' @details
#' Creates a forest plot of effect size estimates (ratios). Estimates that are significantly different
#' from 1 (alpha = 0.05) are marked in orange, others are marked in blue.
#'
#'
#' @param logRr      A numeric vector of effect estimates on the log scale
#' @param seLogRr    The standard error of the log of the effect estimates. Hint: often the standard
#'                   error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                   estimate>))/qnorm(0.025)
#' @param names      A vector containing the names of the drugs or outcomes
#' @param xLabel     The label on the x-axis: the name of the effect estimate
#' @param title      Optional: the main title for the plot
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                   function \code{ggsave} in the ggplot2 package for supported file formats.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' plotForest(negatives$logRr, negatives$seLogRr, negatives$drugName)
#'
#' @export
plotForest <- function(logRr, seLogRr, names, xLabel = "Relative risk", title, fileName = NULL) {
  breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 6)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 5, hjust = 1)
  col <- c(rgb(0, 0, 0.8, alpha = 1), rgb(0.8, 0.4, 0, alpha = 1))
  colFill <- c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0.4, 0, alpha = 0.5))
  data <- data.frame(
    DRUG_NAME = as.factor(names),
    logRr = logRr,
    logLb95Rr = logRr + qnorm(0.025) *
      seLogRr, logUb95Rr = logRr + qnorm(0.975) * seLogRr
  )
  data$significant <- data$logLb95Rr > 0 | data$logUb95Rr < 0
  data$DRUG_NAME <- factor(data$DRUG_NAME, levels = rev(levels(data$DRUG_NAME)))
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$DRUG_NAME,
      y = exp(.data$logRr),
      ymin = exp(.data$logLb95Rr),
      ymax = exp(.data$logUb95Rr),
      colour = .data$significant,
      fill = .data$significant
    )
  ) +
    ggplot2::geom_hline(yintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
    ggplot2::geom_hline(yintercept = 1, size = 0.5) +
    ggplot2::geom_pointrange(shape = 23) +
    ggplot2::scale_colour_manual(values = col) +
    ggplot2::scale_fill_manual(values = colFill) +
    ggplot2::coord_flip(ylim = c(0.25, 10)) +
    ggplot2::scale_y_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) +
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
  if (!missing(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 5, height = 2.5 + length(logRr) * 0.8, dpi = 400)
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
#' \code{plotCalibrationEffect} creates a plot showing the effect of the calibration.
#'
#' @details
#' Creates a plot with the effect estimate on the x-axis and the standard error on the y-axis.
#' Negative controls are shown as blue dots, positive controls as yellow diamonds. The area below the
#' dashed line indicated estimates with p < 0.05. The orange area indicates estimates with calibrated
#' p < 0.05.
#'
#' @param logRrNegatives     A numeric vector of effect estimates of the negative controls on the log
#'                           scale.
#' @param seLogRrNegatives   The standard error of the log of the effect estimates of the negative
#'                           controls.
#' @param logRrPositives     Optional: A numeric vector of effect estimates of the positive controls on the log
#'                           scale.
#' @param seLogRrPositives   Optional: The standard error of the log of the effect estimates of the positive
#'                           controls.
#' @param null               An object representing the fitted null distribution as created by the
#'                           \code{fitNull} or \code{fitMcmcNull} functions. If not provided, a null
#'                           will be fitted before plotting.
#' @param alpha              The alpha for the hypothesis test.
#' @param xLabel             The label on the x-axis: the name of the effect estimate.
#' @param title              Optional: the main title for the plot
#' @param showCis            Show 95 percent credible intervals for the calibrated p = alpha boundary.
#' @param showExpectedAbsoluteSystematicError  Show the expected absolute systematic error. If \code{null} is of
#'                           type \code{mcmcNull} the 95 percent credible interval will also be shown.
#' @param fileName           Name of the file where the plot should be saved, for example 'plot.png'.
#'                           See the function \code{ggsave} in the ggplot2 package for supported file
#'                           formats.
#' @param xLimits            Vector of length 2 for limits of the plot x axis - defaults to 0.25, 10
#' @param yLimits            Vector of length 2 for size limits of the y axis - defaults to 0, 1.5
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' positive <- sccs[sccs$groundTruth == 1, ]
#' plotCalibrationEffect(negatives$logRr, negatives$seLogRr, positive$logRr, positive$seLogRr)
#'
#' @export
plotCalibrationEffect <- function(logRrNegatives,
                                  seLogRrNegatives,
                                  logRrPositives = NULL,
                                  seLogRrPositives = NULL,
                                  null = NULL,
                                  alpha = 0.05,
                                  xLabel = "Relative risk",
                                  title,
                                  showCis = FALSE,
                                  showExpectedAbsoluteSystematicError = FALSE,
                                  fileName = NULL,
                                  xLimits = c(0.25, 10),
                                  yLimits = c(0, 1.5)) {
  if (is.null(null)) {
    if (showCis) {
      null <- fitMcmcNull(logRrNegatives, seLogRrNegatives)
    } else {
      null <- fitNull(logRrNegatives, seLogRrNegatives)
    }
  }
  if (showCis && is(null, "null")) {
    stop("Cannot show credible intervals when using asymptotic null. Please use 'fitMcmcNull' to fit the null")
  }

  if (!is.vector(xLimits) | !length(xLimits) >= 2) {
    stop("xLimits must be a vector of length 2")
  }
  if (!is.vector(yLimits) | !length(yLimits) >= 2) {
    stop("yLimits must be a vector of length 2")
  }

  x <- exp(seq(log(xLimits[1]), log(xLimits[2]), by = 0.01))
  if (is(null, "null")) {
    y <- logRrtoSE(log(x), alpha, null[1], null[2])
  } else {
    chain <- attr(null, "mcmc")$chain
    matrix <- apply(chain, 1, function(null) logRrtoSE(log(x), alpha, null[1], 1 / sqrt(null[2])))
    ys <- apply(matrix, 1, function(se) quantile(se, c(0.025, 0.50, 0.975), na.rm = TRUE))
    rm(matrix)
    y <- ys[2, ]
    yLb <- ys[1, ]
    yUb <- ys[3, ]
  }
  seTheoretical <- sapply(x, FUN = function(x) {
    abs(log(x)) / qnorm(1 - alpha / 2)
  })

  breaks <- c(0.25, 0.5, 1, 2 * (1:ceiling(xLimits[2] / 2.0)))
  if (xLimits[1] < 1) {
    breaks <- c(0.125, breaks)
  }

  checkWithinLimits(yLimits, c(seLogRrNegatives, seLogRrPositives), "yLimits")
  checkWithinLimits(log(xLimits), c(logRrNegatives, logRrPositives), "xLimits")

  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  plot <- ggplot2::ggplot(
    data.frame(x, y, seTheoretical),
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.5) +
    ggplot2::geom_vline(xintercept = 1, size = 1) +
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
    if (is(null, "null")) {
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
    plot <- plot + ggplot2::geom_label(x = log10(0.26), y = 1.49, hjust = "left", vjust = "top", alpha = 0.9, ggplot2::aes(label = .data$text), data = dummy, size = 3.5)
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
      data = data.frame(x = exp(logRrNegatives), y = seLogRrNegatives),
      size = 2,
      alpha = 0.5,
      color = rgb(0, 0, 0.8)
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_x_continuous(xLabel,
      trans = "log10",
      limits = xLimits,
      breaks = breaks,
      labels = breaks
    ) +
    ggplot2::scale_y_continuous("Standard Error") +
    ggplot2::coord_cartesian(ylim = yLimits) +
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
  if (!missing(logRrPositives)) {
    plot <- plot + ggplot2::geom_point(
      shape = 23,
      ggplot2::aes(x = .data$x, y = .data$y),
      data = data.frame(
        x = exp(logRrPositives),
        y = seLogRrPositives
      ),
      size = 4,
      fill = rgb(1, 1, 0),
      alpha = 0.8
    )
  }
  if (!missing(title)) {
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
#' \code{plotCalibration} creates a plot showing the calibration of our calibration procedure
#'
#' @details
#' Creates a calibration plot showing the number of effects with p < alpha for every level of alpha.
#' The empirical calibration is performed using a leave-one-out design: The p-value of an effect is
#' computed by fitting a null using all other negative controls. Ideally, the calibration line should
#' approximate the diagonal. The plot shows both theoretical (traditional) and empirically calibrated
#' p-values.
#'
#' @param logRr            A numeric vector of effect estimates on the log scale
#' @param seLogRr          The standard error of the log of the effect estimates. Hint: often the
#'                         standard error = (log(<lower bound 95 percent confidence interval>) -
#'                         log(<effect estimate>))/qnorm(0.025)
#' @param useMcmc          Use MCMC to estimate the calibrated P-value?
#' @param legendPosition   Where should the legend be positioned? ("none", "left", "right", "bottom",
#'                         "top")
#' @param title            Optional: the main title for the plot
#' @param fileName         Name of the file where the plot should be saved, for example 'plot.png'. See
#'                         the function \code{ggsave} in the ggplot2 package for supported file
#'                         formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' plotCalibration(negatives$logRr, negatives$seLogRr)
#'
#' @export
plotCalibration <- function(logRr,
                            seLogRr,
                            useMcmc = FALSE,
                            legendPosition = "right",
                            title,
                            fileName = NULL) {
  if (any(is.infinite(seLogRr))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting null distribution")
    logRr <- logRr[!is.infinite(seLogRr)]
    seLogRr <- seLogRr[!is.infinite(seLogRr)]
  }
  if (any(is.infinite(logRr))) {
    warning("Estimate(s) with infinite logRr detected. Removing before fitting null distribution")
    seLogRr <- seLogRr[!is.infinite(logRr)]
    logRr <- logRr[!is.infinite(logRr)]
  }
  if (any(is.na(seLogRr))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting null distribution")
    logRr <- logRr[!is.na(seLogRr)]
    seLogRr <- seLogRr[!is.na(seLogRr)]
  }
  if (any(is.na(logRr))) {
    warning("Estimate(s) with NA logRr detected. Removing before fitting null distribution")
    seLogRr <- seLogRr[!is.na(logRr)]
    logRr <- logRr[!is.na(logRr)]
  }

  data <- data.frame(logRr = logRr, SE = seLogRr)
  data$Z <- data$logRr / data$SE
  data$P <- 2 * pmin(pnorm(data$Z), 1 - pnorm(data$Z)) # 2-sided p-value
  data$Y <- sapply(data$P, function(x) {
    sum(data$P < x) / nrow(data)
  })

  data$calibratedP <- vector(length = nrow(data))
  for (i in 1:nrow(data)) {
    dataLeaveOneOut <- data[seq(1, nrow(data)) != i, ]
    if (useMcmc) {
      null <- fitMcmcNull(dataLeaveOneOut$logRr, dataLeaveOneOut$SE)
    } else {
      null <- fitNull(dataLeaveOneOut$logRr, dataLeaveOneOut$SE)
    }
    data$calibratedP[i] <- calibrateP(null, data$logRr[i], data$SE[i], pValueOnly = TRUE)
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
  if (!missing(title)) {
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
#' \code{plotCalibration} creates a plot showing the calibration of our confidence interval
#' calibration procedure
#'
#' @details
#' Creates a calibration plot showing the fraction of effects within the confidence interval. The
#' empirical calibration is performed using a leave-one-out design: The confidence interval of an
#' effect is computed by fitting a null using all other controls. Ideally, the calibration line should
#' approximate the diagonal. The plot shows the coverage for both theoretical (traditional) and
#' empirically calibrated confidence intervals.
#'
#' @param logRr                  A numeric vector of effect estimates on the log scale.
#' @param seLogRr                The standard error of the log of the effect estimates. Hint: often the
#'                               standard error = (log(<lower bound 95 percent confidence interval>) -
#'                               log(<effect estimate>))/qnorm(0.025).
#' @param trueLogRr              The true log relative risk.
#' @param strata                 Variable used to stratify the plot. Set \code{strata = NULL} for no
#'                               stratification.
#' @param legacy                 If true, a legacy error model will be fitted, meaning standard
#'                               deviation is linear on the log scale. If false, standard deviation
#'                               is assumed to be simply linear.
#' @param evaluation             A data frame as generated by the \code{\link{evaluateCiCalibration}}
#'                               function. If provided, the logRr, seLogRr, trueLogRr, strata, and legacy
#'                               arguments will be ignored.
#' @param crossValidationGroup   What should be the unit for the cross-validation? By default the unit
#'                               is a single control, but a different grouping can be provided, for
#'                               example linking a negative control to synthetic positive controls
#'                               derived from that negative control.
#' @param legendPosition         Where should the legend be positioned? ("none", "left", "right",
#'                               "bottom", "top").
#' @param title                  Optional: the main title for the plot
#' @param fileName               Name of the file where the plot should be saved, for example
#'                               'plot.png'. See the function \code{ggsave} in the ggplot2 package for
#'                               supported file formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' \dontrun{
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' plotCiCalibration(data$logRr, data$seLogRr, data$trueLogRr)
#' }
#' @export
plotCiCalibration <- function(logRr,
                              seLogRr,
                              trueLogRr,
                              strata = as.factor(trueLogRr),
                              crossValidationGroup = 1:length(logRr),
                              legacy = FALSE,
                              evaluation,
                              legendPosition = "top",
                              title,
                              fileName = NULL) {
  if (missing(evaluation) || is.null(evaluation)) {
    evaluation <- evaluateCiCalibration(
      logRr = logRr,
      seLogRr = seLogRr,
      trueLogRr = trueLogRr,
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
    ggplot2::facet_grid(. ~ trueRr) +
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
  if (!missing(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    width <- 1 + 2 * length(levels(evaluation$trueRr))
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
#' Creates a forest plot of effect size estimates (ratios). Estimates that are significantly different
#' from the true value (alpha = 0.05) are marked in orange, others are marked in blue.
#'
#' @param logRr       A numeric vector of effect estimates on the log scale.
#' @param seLogRr     The standard error of the log of the effect estimates. Hint: often the standard
#'                    error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                    estimate>))/qnorm(0.025).
#' @param trueLogRr   A vector of the true effect sizes.
#' @param xLabel      The label on the x-axis: the name of the effect estimate.
#' @param title       Optional: the main title for the plot
#' @param fileName    Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                    function \code{ggsave} in the ggplot2 package for supported file formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' plotTrueAndObserved(data$logRr, data$seLogRr, data$trueLogRr)
#'
#' @export
plotTrueAndObserved <- function(logRr,
                                seLogRr,
                                trueLogRr,
                                xLabel = "Relative risk",
                                title,
                                fileName = NULL) {
  breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 6)
  col <- c(rgb(0, 0, 0.8, alpha = 1), rgb(0.8, 0.4, 0, alpha = 1))
  colFill <- c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0.4, 0, alpha = 0.5))
  data <- data.frame(
    logRr = logRr,
    logLb95Rr = logRr + qnorm(0.025) * seLogRr,
    logUb95Rr = logRr + qnorm(0.975) * seLogRr,
    trueLogRr = trueLogRr,
    trueRr = round(exp(trueLogRr), 2)
  )
  data$significant <- data$logLb95Rr > data$trueLogRr | data$logUb95Rr < data$trueLogRr
  data <- data[order(data$trueLogRr, data$logRr), ]
  data$order <- 1:nrow(data)
  coverage <- aggregate(!significant ~ trueRr, data = data, mean)
  names(coverage)[2] <- "coverage"
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = exp(.data$logRr),
      y = .data$order,
      xmin = exp(.data$logLb95Rr),
      xmax = exp(.data$logUb95Rr),
      colour = .data$significant,
      fill = .data$significant
    )
  ) +
    ggplot2::geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmax = .data$trueRr, xmin = .data$trueRr), height = 1, color = rgb(0, 0, 0), size = 1) +
    ggplot2::geom_errorbarh(height = 0) +
    ggplot2::geom_point(shape = 21, size = 1.5) +
    ggplot2::scale_colour_manual(values = col) +
    ggplot2::scale_fill_manual(values = colFill) +
    ggplot2::coord_cartesian(xlim = c(0.25, 10)) +
    ggplot2::scale_x_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) +
    ggplot2::facet_grid(trueRr ~ ., scales = "free_y", space = "free") +
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
  if (!missing(title)) {
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
#' @param mcmcNull   An object of type \code{mcmcNull} as generated using the \code{fitMcmcNull}
#'                   function.
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                   function \code{ggsave} in the ggplot2 package for supported file formats.
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
#' \code{plotCiCoverage} creates a plot showing the coverage before and after confidence interval
#' calibration at various widths of the confidence interval.
#'
#' @details
#' Creates a plot showing the fraction of effects above, within, and below the confidence interval. The
#' empirical calibration is performed using a leave-one-out design: The confidence interval of an
#' effect is computed by fitting a null using all other controls. The plot shows the coverage for
#' both theoretical (traditional) and empirically calibrated confidence intervals.
#'
#' @param logRr                  A numeric vector of effect estimates on the log scale.
#' @param seLogRr                The standard error of the log of the effect estimates. Hint: often the
#'                               standard error = (log(<lower bound 95 percent confidence interval>) -
#'                               log(<effect estimate>))/qnorm(0.025).
#' @param trueLogRr              The true log relative risk.
#' @param strata                 Variable used to stratify the plot. Set \code{strata = NULL} for no
#'                               stratification.
#' @param legacy                 If true, a legacy error model will be fitted, meaning standard
#'                               deviation is linear on the log scale. If false, standard deviation
#'                               is assumed to be simply linear.
#' @param evaluation             A data frame as generated by the \code{\link{evaluateCiCalibration}}
#'                               function. If provided, the logRr, seLogRr, trueLogRr, strata, and legacy
#'                               arguments will be ignored.
#' @param crossValidationGroup   What should be the unit for the cross-validation? By default the unit
#'                               is a single control, but a different grouping can be provided, for
#'                               example linking a negative control to synthetic positive controls
#'                               derived from that negative control.
#' @param legendPosition         Where should the legend be positioned? ("none", "left", "right",
#'                               "bottom", "top").
#' @param title                  Optional: the main title for the plot
#' @param fileName               Name of the file where the plot should be saved, for example
#'                               'plot.png'. See the function \code{ggsave} in the ggplot2 package for
#'                               supported file formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' \dontrun{
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' plotCiCoverage(data$logRr, data$seLogRr, data$trueLogRr)
#' }
#' @export
plotCiCoverage <- function(logRr,
                           seLogRr,
                           trueLogRr,
                           strata = as.factor(trueLogRr),
                           crossValidationGroup = 1:length(logRr),
                           legacy = FALSE,
                           evaluation,
                           legendPosition = "top",
                           title,
                           fileName = NULL) {
  if (missing(evaluation) || is.null(evaluation)) {
    evaluation <- evaluateCiCalibration(
      logRr = logRr,
      seLogRr = seLogRr,
      trueLogRr = trueLogRr,
      strata = strata,
      crossValidationGroup = crossValidationGroup,
      legacy = legacy
    )
  }
  evaluation$label <- factor(evaluation$label, levels = rev(c("Below confidence interval", "Within confidence interval", "Above confidence interval")))
  evaluation$`Confidence interval calculation` <- factor(evaluation$`Confidence interval calculation`, levels = rev(c("Calibrated", "Uncalibrated")))
  evaluation$trueRr <- paste0("True effect size = ", evaluation$trueRr)
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
    ggplot2::facet_grid(`Confidence interval calculation` ~ trueRr) +
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
  if (!missing(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    width <- 1 + 1.8 * length(levels(as.factor(evaluation$trueRr)))
    ggplot2::ggsave(fileName, plot, width = width, height = 5, dpi = 400)
  }
  return(plot)
}

#' Plot the systematic error model
#'
#' @description
#' \code{plotErrorModel} creates a plot showing the systematic error model.
#'
#' @details
#' Creates a plot with the true effect size on the x-axis, and the mean plus and minus the standard
#' deviation shown on the y-axis. Also shown are simple error models fitted at each true relative
#' risk in the input.
#'
#' @param logRr                  A numeric vector of effect estimates on the log scale.
#' @param seLogRr                The standard error of the log of the effect estimates. Hint: often the
#'                               standard error = (log(<lower bound 95 percent confidence interval>) -
#'                               log(<effect estimate>))/qnorm(0.025).
#' @param trueLogRr              The true log relative risk.
#' @param legacy                 If true, a legacy error model will be fitted, meaning standard
#'                               deviation is linear on the log scale. If false, standard deviation
#'                               is assumed to be simply linear.
#' @param title              Optional: the main title for the plot
#' @param fileName           Name of the file where the plot should be saved, for example 'plot.png'.
#'                           See the function \code{ggsave} in the ggplot2 package for supported file
#'                           formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' plotErrorModel(data$logRr, data$seLogRr, data$trueLogRr)
#' @export
plotErrorModel <- function(logRr, seLogRr, trueLogRr, title, legacy = FALSE, fileName = NULL) {
  model <- fitSystematicErrorModel(
    logRr = logRr,
    seLogRr = seLogRr,
    trueLogRr = trueLogRr,
    legacy = legacy
  )
  fitSimpleModel <- function(selectedTrueLogRr) {
    idx <- trueLogRr == selectedTrueLogRr
    return(fitNull(logRr[idx], seLogRr[idx]))
  }
  simpleModels <- as.data.frame(t(sapply(unique(trueLogRr), fitSimpleModel)))
  simpleModels$trueRr <- exp(unique(trueLogRr))
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
  data <- data.frame(x = x, y = y, ymin = ymin, ymax = ymax)
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      ymin = .data$ymin,
      ymax = .data$ymax
    )
  ) +
    ggplot2::geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
    ggplot2::geom_hline(yintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
    ggplot2::geom_vline(xintercept = 1, size = 1) +
    ggplot2::geom_hline(yintercept = 1, size = 1) +
    ggplot2::geom_ribbon(fill = rgb(0, 0, 0.8), alpha = 0.3) +
    ggplot2::geom_line(color = rgb(0, 0, 0.8)) +
    ggplot2::geom_errorbar(ggplot2::aes(x = .data$trueRr), width = 0.1, size = 1, color = rgb(0, 0, 0.8), data = simpleModels) +
    ggplot2::scale_x_continuous("True effect size",
      trans = "log10",
      breaks = breaks,
      labels = breaks
    ) +
    ggplot2::scale_y_continuous("Systematic error mean (plus and minus one SD)",
      trans = "log10",
      breaks = breaks,
      labels = breaks
    ) +
    ggplot2::coord_cartesian(xlim = c(0.25, 10), ylim = c(0.25, 10)) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA),
      panel.grid.major = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
  if (!missing(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }
  return(plot)
}

plotIsobars <- function(null, alpha, xLabel = "Relative risk", seLogRrPositives) {
  if (is(null, "mcmcNull")) {
    null <- c(null[1], 1 / sqrt(null[2]))
  }
  x <- exp(seq(log(0.25), log(10), by = 0.01))
  y <- sapply(x, FUN = function(x) abs(log(x)) / qnorm(1 - alpha / 2))
  thresholds <- c(0.05, 0.25, 0.5, 0.75)
  isoBars <- lapply(thresholds, function(threshold) {
    data.frame(
      x = x,
      y = logRrtoSE(log(x), threshold, null[1], null[2]),
      ymax = y,
      threshold = threshold
    )
  })
  breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  plot <- ggplot2::ggplot(
    data.frame(x, y),
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_vline(xintercept = breaks, colour = rgb(0, 0, 0), lty = 1, size = 0.5, alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 0:4 / 4, colour = rgb(0, 0, 0), lty = 1, size = 0.5, alpha = 0.2) +
    ggplot2::geom_vline(xintercept = 1, size = 1) +
    ggplot2::geom_line(
      colour = rgb(0, 0, 0),
      linetype = "dashed",
      size = 1,
      alpha = 0.5
    ) +
    ggplot2::scale_x_continuous(xLabel,
      trans = "log10",
      breaks = breaks,
      labels = breaks
    ) +
    ggplot2::scale_y_continuous("Standard Error") +
    ggplot2::coord_cartesian(xlim = c(0.25, 10), ylim = c(0, 1)) +
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
    isoBarLeft <- isoBar[isoBar$x < exp(null[1]), ]
    isoBarRight <- isoBar[isoBar$x > exp(null[1]), ]
    labelData <- data.frame(
      x = c(
        isoBarLeft$x[which.min(abs(isoBarLeft$y - 0.625))],
        isoBarRight$x[which.min(abs(isoBarRight$y - 0.625))]
      ),
      y = 0.625,
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
  if (!missing(seLogRrPositives)) {
    plot <- plot + ggplot2::geom_hline(yintercept = seLogRrPositives)
  }
  return(plot)
}

#' Plot the expected type 1 error as a function of standard error
#'
#' @description
#' \code{plotExpectedType1Error} creates a plot showing the expected type 1 error as a function of standard error.
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
#' @param logRrNegatives     A numeric vector of effect estimates of the negative controls on the log
#'                           scale.
#' @param seLogRrNegatives   The standard error of the log of the effect estimates of the negative
#'                           controls.
#' @param seLogRrPositives   The standard error of the log of the effect estimates of the positive
#'                           controls.
#' @param alpha              The alpha (nominal type 1 error) to be used.
#' @param null               An object representing the fitted null distribution as created by the
#'                           \code{fitNull} function. If not provided, a null will be fitted before
#'                           plotting.
#' @param xLabel             If showing effect sizes, what label should be used for the effect size axis?
#' @param title              Optional: the main title for the plot
#' @param showCis            Show 95 percent credible intervals for the expected type 1 error.
#' @param showEffectSizes    Show the expected effect sizes alongside the expected type 1 error?
#' @param fileName           Name of the file where the plot should be saved, for example 'plot.png'.
#'                           See the function \code{ggsave} in the ggplot2 package for supported file
#'                           formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' positive <- sccs[sccs$groundTruth == 1, ]
#' plotExpectedType1Error(negatives$logRr, negatives$seLogRr, positive$seLogRr)
#'
#' @export
plotExpectedType1Error <- function(logRrNegatives,
                                   seLogRrNegatives,
                                   seLogRrPositives,
                                   alpha = 0.05,
                                   null = NULL,
                                   xLabel = "Relative risk",
                                   title,
                                   showCis = FALSE,
                                   showEffectSizes = FALSE,
                                   fileName = NULL) {
  if (is.null(null)) {
    if (showCis) {
      null <- fitMcmcNull(logRrNegatives, seLogRrNegatives)
    } else {
      null <- fitNull(logRrNegatives, seLogRrNegatives)
    }
  }
  if (showCis && is(null, "null")) {
    stop("Cannot show credible intervals when using asymptotic null. Please use 'fitMcmcNull' to fit the null")
  }

  se <- (0:100) / 100
  if (is(null, "null")) {
    mean <- null[1]
    sd <- sqrt(se^2 + null[2]^2)
    type1Error <- pnorm(qnorm(1 - alpha / 2, 0, se), mean, sd, lower.tail = FALSE) +
      pnorm(qnorm(alpha / 2, 0, se), mean, sd, lower.tail = TRUE)
  } else {
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
  }
  breaks <- 0:4 / 4
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 12, hjust = 1)
  plot <- ggplot2::ggplot(
    data.frame(se, type1Error),
    ggplot2::aes(x = .data$se, y = .data$type1Error)
  ) +
    ggplot2::geom_vline(xintercept = breaks, colour = rgb(0, 0, 0), lty = 1, size = 0.5, alpha = 0.2) +
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
      se = 0.875 + showEffectSizes * 0.125,
      type1Error = alpha + showEffectSizes * 0.04
    ),
    parse = TRUE
  ) +
    ggplot2::scale_x_continuous("Standard error",
      limits = c(0, 1),
      breaks = breaks,
      labels = breaks
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
  if (!missing(seLogRrPositives)) {
    yPos <- sapply(seLogRrPositives, FUN = function(x) {
      type1Error[which.min(abs(x - se))]
    })
    posData <- data.frame(
      x = seLogRrPositives,
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
    plot2 <- plotIsobars(null = null, alpha = alpha, xLabel = xLabel, seLogRrPositives = seLogRrPositives)
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
    if (missing(title)) {
      title <- NULL
    }
    plot <- gridExtra::grid.arrange(grobs[[1]], grobs[[2]], nrow = 1, ncol = 2, widths = c(2, 1), top = title)
    if (!is.null(fileName)) {
      ggplot2::ggsave(fileName, plot, width = 10, height = 5, dpi = 400)
    }
  } else {
    if (!missing(title)) {
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
#' @param logRr       A numeric vector of effect estimates on the log scale.
#' @param seLogRr     The standard error of the log of the effect estimates. Hint: often the standard
#'                    error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                    estimate>))/qnorm(0.025).
#' @param trueLogRr   A vector of the true effect sizes.
#' @param legacy      If true, a legacy error model will be fitted, meaning standard
#'                    deviation is linear on the log scale. If false, standard deviation
#'                    is assumed to be simply linear.
#' @param model       The fitted systematic error model. If not provided, it will be fitted on the
#'                    provided data.
#' @param xLabel      The label on the x-axis: the name of the effect estimate.
#' @param title       Optional: the main title for the plot
#' @param fileName    Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                    function \code{ggsave} in the ggplot2 package for supported file formats.
#'
#' @return
#' A Ggplot object. Use the \code{ggsave} function to save to file.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' plotCiCalibrationEffect(data$logRr, data$seLogRr, data$trueLogRr)
#'
#' @export
plotCiCalibrationEffect <- function(logRr,
                                    seLogRr,
                                    trueLogRr,
                                    legacy = FALSE,
                                    model = NULL,
                                    xLabel = "Relative risk",
                                    title,
                                    fileName = NULL) {
  alpha <- 0.05
  if (is.null(model)) {
    model <- fitSystematicErrorModel(
      logRr = logRr,
      seLogRr = seLogRr,
      trueLogRr = trueLogRr,
      estimateCovarianceMatrix = FALSE,
      legacy = legacy
    )
  } else {
    legacy <- (names(model)[3] == "logSdIntercept")
  }
  d <- data.frame(
    logRr = logRr,
    seLogRr = seLogRr,
    trueLogRr = trueLogRr,
    trueRr = exp(trueLogRr),
    logCi95lb = logRr + qnorm(0.025) * seLogRr,
    logCi95ub = logRr + qnorm(0.975) * seLogRr
  )
  d <- d[!is.na(d$logRr), ]
  d <- d[!is.na(d$seLogRr), ]
  if (nrow(d) == 0) {
    return(NULL)
  }
  d$Group <- as.factor(d$trueRr)
  d$Significant <- d$logCi95lb > d$trueLogRr | d$logCi95ub < d$trueLogRr

  temp1 <- aggregate(Significant ~ trueRr, data = d, length)
  temp2 <- aggregate(Significant ~ trueRr, data = d, mean)
  temp1$nLabel <- paste0(formatC(temp1$Significant, big.mark = ","), " estimates")
  temp1$Significant <- NULL
  temp2$meanLabel <- paste0(
    formatC(100 * (1 - temp2$Significant), digits = 1, format = "f"),
    "% of CIs includes ",
    temp2$trueRr
  )
  temp2$Significant <- NULL
  dd <- merge(temp1, temp2)

  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  theme <- ggplot2::element_text(colour = "#000000", size = 10)
  themeRA <- ggplot2::element_text(colour = "#000000", size = 10, hjust = 1)

  d$Group <- paste("True", tolower(xLabel), "=", d$trueRr)
  dd$Group <- paste("True", tolower(xLabel), "=", dd$trueRr)

  x <- seq(log(0.1), log(10), by = 0.01)
  calBounds <- data.frame()
  for (i in 1:nrow(dd)) {
    mu <- model[1] + model[2] * log(dd$trueRr[i])
    if (legacy) {
      sigma <- exp(model[3] + model[4] * log(dd$trueRr[i]))
    } else {
      sigma <- model[3] + model[4] * abs(log(dd$trueRr[i]))
    }
    calBounds <- rbind(
      calBounds,
      data.frame(
        logRr = x,
        seLogRr = logRrtoSE(x, alpha, mu, sigma),
        Group = dd$Group[i]
      )
    )
  }
  plot <- ggplot2::ggplot(d, ggplot2::aes(x = .data$logRr, y = .data$seLogRr)) +
    ggplot2::geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.5) +
    ggplot2::geom_area(
      fill = rgb(1, 0.5, 0, alpha = 0.5),
      color = rgb(1, 0.5, 0),
      size = 1,
      alpha = 0.5, data = calBounds
    ) +
    ggplot2::geom_abline(ggplot2::aes(intercept = (-log(.data$trueRr)) / qnorm(0.025), slope = 1 / qnorm(0.025)), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5, data = dd) +
    ggplot2::geom_abline(ggplot2::aes(intercept = (-log(.data$trueRr)) / qnorm(0.975), slope = 1 / qnorm(0.975)), colour = rgb(0, 0, 0), linetype = "dashed", size = 1, alpha = 0.5, data = dd) +
    ggplot2::geom_point(
      shape = 16,
      size = 2,
      alpha = 0.5,
      color = rgb(0, 0, 0.8)
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_label(x = log(0.15), y = 0.95, alpha = 1, hjust = "left", ggplot2::aes(label = .data$nLabel), size = 3.5, data = dd) +
    ggplot2::geom_label(x = log(0.15), y = 0.8, alpha = 1, hjust = "left", ggplot2::aes(label = .data$meanLabel), size = 3.5, data = dd) +
    ggplot2::scale_x_continuous(xLabel, limits = log(c(0.1, 10)), breaks = log(breaks), labels = breaks) +
    ggplot2::scale_y_continuous("Standard Error") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
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
  if (!missing(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 2 + 3.5 * nrow(dd), height = 2.8, dpi = 400)
  }
  return(plot)
}
