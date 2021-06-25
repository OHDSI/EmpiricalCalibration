library(dplyr)
library(EmpiricalCalibration)
# Simulate from Poisson with unbalanced exposure groups to achieve non-normal likelihood functions

nControls <- 100
timeExposed <- 1
timeUnexposed <- 10
backgroundRate <- 5
mean <- 0.2
sd <- 0.2

# Simulate data ------------------------
set.seed(123)
ncsNormal <- tibble(logRr = rep(NA, nControls),
                    seLogRr = rep(NA, nControls))
profiles <- list()
for (i in 1:nControls) {
  expTheta <- exp(rnorm(1, mean, sd))
  
  # Adjust rate in unexposed to keep overall power constant:
  # rateUnexposed <- backgroundRate * ((timeExposed + timeUnexposed) / timeExposed) / (expTheta + (timeUnexposed / timeExposed))
  # rateExposed <- expTheta * rateUnexposed
  # outcome <- rpois(2, c(rateExposed * timeExposed, rateUnexposed * timeUnexposed))
  
  outcome <- rpois(2, c(expTheta * backgroundRate * timeExposed, backgroundRate * timeUnexposed))
  exposure <- c(1, 0)
  time <- c(timeExposed, timeUnexposed)
  cyclopsData <- Cyclops::createCyclopsData(outcome ~ exposure + offset(log(time)), modelType = "pr") 
  fit <- Cyclops::fitCyclopsModel(cyclopsData)
  
  profiles[[i]] <- Cyclops::getCyclopsProfileLogLikelihood(fit, "exposure", bounds = c(log(0.1), log(10)))
  # plot(profile$point, profile$value)
  ncsNormal$logRr[i] <- coef(fit)["exposure"]
  ci <- confint(fit, "exposure")
  ncsNormal$seLogRr[i] <- (ci[3] - ci[2]) / (2 * qnorm(0.975))
}

# Fit empirical null ----------------------------------------

# Using normal approximation
ncsNormal <- ncsNormal[!is.na(ncsNormal$seLogRr), ]
null <- fitNull(ncsNormal$logRr, ncsNormal$seLogRr)
null

null <- fitNullNonNormalLl(ncsNormal)
null

# Using grid approximation

plotCalibrationEffect(ncsNormal$logRr, ncsNormal$seLogRr, null = null)

null <- fitNullNonNormalLl(profiles)
null

plotCalibrationEffect(ncsNormal$logRr, ncsNormal$seLogRr, null = null)
