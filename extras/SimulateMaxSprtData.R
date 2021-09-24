# This code is used to precompute some of the results in the MaxSPRT vignette

set.seed(123)
maxSprtSimulationData <- simulateMaxSprtData()

cvs <- c()
for (i in 1:50) {
  dataOutcomeI <- maxSprtSimulationData[maxSprtSimulationData$outcomeId == i, ]
  outcomesPerLook <- aggregate(outcome ~ lookTime, dataOutcomeI, sum)
  # Need incremental outcomes per look:
  outcomesPerLook <- outcomesPerLook$outcome[order(outcomesPerLook$lookTime)]
  outcomesPerLook[2:10] <- outcomesPerLook[2:10] - outcomesPerLook[1:9]
  
  exposedTime <- sum(dataOutcomeI$time[dataOutcomeI$exposure == TRUE & 
                                         dataOutcomeI$lookTime == 100])
  unexposedTime <- sum(dataOutcomeI$time[dataOutcomeI$exposure == FALSE & 
                                            dataOutcomeI$lookTime == 100])
  cv <- computeCvBinomial(groupSizes = outcomesPerLook,
                          z = unexposedTime / exposedTime,
                          minimumEvents = 1,
                          alpha = 0.05)
  cvs[i] <- cv
}
saveRDS(cvs, "vignettes/cvs.rds")


allLlrs <- data.frame()
for (t in unique(maxSprtSimulationData$lookTime)) {
  
  # Compute likelihood profiles for all negative controls:
  negativeControlProfilesLookTt <- list()
  dataLookTt <- maxSprtSimulationData[maxSprtSimulationData$lookTime == t, ]
  for (i in 1:50) {
    dataOutcomeIlookTt <- dataLookTt[dataLookTt$outcomeId == i, ]
    cyclopsData <- createCyclopsData(Surv(time, outcome) ~ exposure , 
                                     modelType = "cox", 
                                     data = dataOutcomeIlookTt)
    fit <- fitCyclopsModel(cyclopsData)
    llProfile <- getCyclopsProfileLogLikelihood(object = fit,
                                                parm = "exposureTRUE",
                                                bounds = log(c(0.1, 10)))
    negativeControlProfilesLookTt[[i]] <- llProfile
  }
  
  # Fit null distribution:
  nullLookTt <- fitNullNonNormalLl(negativeControlProfilesLookTt)

  # Compute calibrated and uncalibrated LLRs for all negative controls:
  llrsLookTt <- c()
  calibratedLlrsLookTt <- c()
  for (i in 1:50) {
    dataOutcomeIlookTt <- dataLookTt[dataLookTt$outcomeId == i, ]
    cyclopsData <- createCyclopsData(Surv(time, outcome) ~ exposure , 
                                     modelType = "cox", 
                                     data = dataOutcomeIlookTt)
    fit <- fitCyclopsModel(cyclopsData)
    llProfile <- getCyclopsProfileLogLikelihood(object = fit,
                                                parm = "exposureTRUE",
                                                bounds = log(c(0.1, 10)))
    
    # Calibrated LLR:
    calibrateLlr <- calibrateLlr(null = nullLookTt, 
                 likelihoodApproximation = llProfile)
    calibratedLlrsLookTt[i] <- calibrateLlr
    
    # Uncalibrated LLR:
    llNull <- getCyclopsProfileLogLikelihood(object = fit,
                                             parm = "exposureTRUE",
                                             x = 0)$value
    if (fit$return_flag == "ILLCONDITIONED" || coef(fit) < 0) {
      llr <- 0
    } else {
      llr <- fit$log_likelihood - llNull
    }
    llrsLookTt[i] <- llr
  }
  
  # Store in a data frame:
  allLlrs <- rbind(allLlrs,
                   data.frame(t = t,
                              outcomeId = 1:50,
                              llr = llrsLookTt,
                              calibrateLlr = calibratedLlrsLookTt))
  
}
saveRDS(allLlrs, "vignettes/allLlrs.rds")

