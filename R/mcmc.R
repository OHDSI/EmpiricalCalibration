# @file mcmc.R
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

proposalFunction <- function(param, scale = c(0.1, 0.1), fixed = NULL) {
  dim <- length(param)
  draw <- rnorm(dim, mean = param, sd = scale)
  if (!is.null(fixed)) {
    draw[fixed] <- param[fixed]
  }
  draw[2] <- abs(draw[2])
  return(draw)
}

runMetropolisMcmc <- function(startValue, ll, iterations, scale, logRr, seLogRr, fixed = NULL) {
  dim <- length(startValue)
  chain <- array(dim = c(iterations + 1, dim))
  logLik <- array(dim = c(iterations + 1, 1))
  acc <- array(dim = c(iterations + 1, 1))

  logLik[1] <- -ll(startValue, logRr, seLogRr)
  chain[1, ] <- c(startValue)
  acc[1] <- 1

  for (i in 1:iterations) {
    # print(paste('itr =', i))
    proposal <- proposalFunction(chain[i, ], scale = scale, fixed = fixed)

    newLogLik <- tryCatch(-ll(proposal, logRr, seLogRr), error = function(e) {
      -1e+10
    })

    prob <- exp(newLogLik - logLik[i])
    if (runif(1) < prob) {
      chain[i + 1, ] <- proposal
      logLik[i + 1] <- newLogLik
      acc[i + 1] <- 1
    } else {
      chain[i + 1, ] <- chain[i, ]
      logLik[i + 1] <- logLik[i]
      acc[i + 1] <- 0
    }
  }
  result <- list(logLik = logLik, chain = chain, acc = acc)
  return(result)
}
