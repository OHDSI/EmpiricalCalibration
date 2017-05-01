# @file Simulation.R
#
# Copyright 2017 Observational Health Data Sciences and Informatics
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
