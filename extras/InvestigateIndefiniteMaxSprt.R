library(EmpiricalCalibration)
library(Cyclops)
library(survival)

alpha <- 0.05
beta <- 0.20

nSimulations <- 100
maxT <- 1000
irr <- 1.5

cvSignal <- log((1 - beta) / alpha)
cvSafe <- log(beta / (1 - alpha))





cvSignal <- 3.09
cvSafe <- -0.65
# llrs <- data.frame(sim = rep(1:nSimulations, each = maxT),
#                    t = rep(1:maxT, times = nSimulations),
#                    llr)
signals <- 0
safe <- 0
notDone <- 0
for (i in 1:nSimulations) {
  observed <- 0
  expected <- 0
  for (t in 1:maxT) {
    groupSize <- round(runif(1, 1, 10))
    expected <- expected + groupSize
    observed <- observed + rpois(1, groupSize * irr)
    
    # writeLines(sprintf("%0.2f at t = %d", observed / expected, t))
    
    # SPRT using H1: IRR = 2
    # llr <- dpois(observed, expected * 2, log = TRUE) - dpois(observed, expected, log = TRUE)
    # if (llr > cvSignal) {
    #   signals <- signals + 1
    #   break
    # } else if (llr < cvSafe) {
    #   safe <- safe + 1
    #   break
    # }
    
    # MaxSPRT    
    # llr <- dpois(observed, observed, log = TRUE) - dpois(observed, expected, log = TRUE)
    # if (observed < expected) {
    #   llr <- -llr
    # }
    # if (observed >= expected) {
    #   if (llr > cvSignal) {
    #     signals <- signals + 1
    #     break
    #   }
    # } else {
    #   if (llr < cvSafe) {
    #     safe <- safe + 1
    #     break
    #   }
    # }
    # llrs$llr[(i - 1) * maxT + t] <- llr
    
    # Confidence intervals
    if (qpois(0.975, observed) / expected < 1.5) {
      safe <- safe + 1
      break
    } else if (qpois(0.025, observed) / expected > 1){
      signals <- signals + 1
      break
    }
  }
  if (t == maxT) {
    notDone <- notDone + 1
  }
}
signals / nSimulations
safe / nSimulations
notDone / nSimulations

library(ggplot2)
ggplot(llrs, aes(x = t, y = llr, group = sim)) +
  geom_line(alpha = 0.2)


data <- simulateMaxSprtData(looks = 10, 
                            nullMu = 0, 
                            nullSigma = 0,
                            numberOfNegativeControls = 100,
                            numberOfPositiveControls = 0)

signals <- 0
for (outcomeId in unique(data$outcomeId)) {
  writeLines(sprintf("Outcome %s", outcomeId))
  signal <- FALSE
  for (time in unique(data$lookTime)) {
    subset <- data[data$outcomeId == outcomeId & data$lookTime == time, ]
    cyclopsData <- createCyclopsData(Surv(time, outcome) ~ exposure , 
                                     modelType = "cox", 
                                     data = subset)
    fit <- fitCyclopsModel(cyclopsData)
    if (fit$return_flag == "SUCCESS") {
      llNull <- getCyclopsProfileLogLikelihood(object = fit,
                                               parm = "exposureTRUE",
                                               x = 0)$value
      llr <- fit$log_likelihood - llNull
      
      if (llr > cvSignal) {
        signal <- TRUE
      }
    }
  }
  if (signal) {
    signals <- signals + 1
  }
}

# Fit null distribution:
nullLookTt <- fitNullNonNormalLl(negativeControlProfilesLookTt)
