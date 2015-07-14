# @file CiUsingMcmc.R
#
# Copyright 2015 Observational Health Data Sciences and Informatics
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

#' Compute p-value confidence intervals using MCMC
#'
#' @details
#' This is an experimental function for computing the 95 percent confidence interval of a calibrated
#' p-value using Markov-Chain Monte Carlo (MCMC). This should give better estimates than the default
#' function when the standard deviation of the error distribution is close to zero.
#'
#' @param logRrNegatives     A numeric vector of effect estimates of the negative controls on the log
#'                           scale.
#' @param seLogRrNegatives   The standard error of the log of the effect estimates of the negative
#'                           controls.
#' @param logRrPositives     A numeric vector of effect estimates of the positive controls on the log
#'                           scale.
#' @param seLogRrPositives   The standard error of the log of the effect estimates of the positive
#'                           controls.
#' @param scale              A vector of two numbers representing the scale of the likelihood space
#'                           around the mean and standard deviation of the error distribution,
#'                           respecitively.
#' @param iter               Number of iterations of the MCMC.
#'
#' @examples
#' controls <- simulateControls()
#' p <- calibratePWithCiUsingMcmc(controls$logRr,
#'                                controls$seLogRr,
#'                                controls$logRr[1],
#'                                controls$seLogRr[1],
#'                                scale = c(0.05, 25),
#'                                iter = 10000)
#' mcmc <- attr(p, "mcmc")
#' mean(mcmc$acc)  # Acceptance rate
#' plot(ts(mcmc$chain[, 1]))  # Trace for the mean
#' plot(ts(mcmc$chain[, 2]))  # Trace for the precision
#' mean(mcmc$chain[, 1])
#' mean(mcmc$chain[, 2])
#'
#'
#' @export
calibratePWithCiUsingMcmc <- function(logRrNegatives,
                                      seLogRrNegatives,
                                      logRrPositives,
                                      seLogRrPositives,
                                      scale = c(0.1, 200),
                                      iter = 10000) {
  gaussianProduct <- function(mu1, mu2, sd1, sd2) {
    (2 * pi)^(-1/2) * (sd1^2 + sd2^2)^(-1/2) * exp(-(mu1 - mu2)^2/(2 * (sd1^2 + sd2^2)))
  }

  LL <- function(theta, estimate, se) {
    result <- 0
    for (i in 1:length(estimate)) result <- result + log(gaussianProduct(estimate[i],
                                                                         theta[1],
                                                                         se[i],
                                                                         1/sqrt(theta[2])))
    if (is.infinite(result))
      result <- -99999
    result <- result + dgamma(theta[2], shape = 0.001, rate = 0.001, log = TRUE)
    return(-result)
  }

  fit <- optim(c(0, 0.1), LL, estimate = logRrNegatives, se = seLogRrNegatives)
  mcmc <- runMetropolisMcmc(fit$par, LL, iterations = iter, scale, logRrNegatives, seLogRrNegatives)

  adjustedP <- data.frame(p = rep(1, length(logRrPositives)), lb95ci = 0, ub95ci = 0)
  for (i in 1:length(logRrPositives)) {
    P_upper_bound <- pnorm((mcmc$chain[,
                            1] - logRrPositives[i])/sqrt((1/sqrt(mcmc$chain[,
                                                                 2]))^2 + seLogRrPositives[i]^2))
    P_lower_bound <- pnorm((logRrPositives[i] - mcmc$chain[,
                            1])/sqrt((1/sqrt(mcmc$chain[, 2]))^2 + seLogRrPositives[i]^2))
    p <- P_upper_bound
    p[P_lower_bound < p] <- P_lower_bound[P_lower_bound < p]
    p <- p * 2
    adjustedP$p[i] <- quantile(p, 0.5)
    adjustedP$lb95ci[i] <- quantile(p, 0.025)
    adjustedP$ub95ci[i] <- quantile(p, 0.975)
  }
  attr(adjustedP, "mcmc") <- mcmc
  return(adjustedP)
}
