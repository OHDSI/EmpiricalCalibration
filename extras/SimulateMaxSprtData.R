# This code is used to precompute some of the results in the MaxSPRT vignette

set.seed(123)
maxSprtSimulationData <- simulateMaxSprtData(
  n = 10000,
  pExposure = 0.5,
  backgroundHazard = 0.001,
  tar = 10,
  nullMu = 0.2,
  nullSigma = 0.2,
  maxT = 100,
  looks = 10,
  numberOfNegativeControls = 50,
  numberOfPositiveControls = 1,
  positiveControlEffectSize = 4
)

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
  cv <- computeCvBinomial(
    groupSizes = outcomesPerLook,
    z = unexposedTime / exposedTime,
    minimumEvents = 1,
    alpha = 0.05
  )
  cvs[i] <- cv
}
saveRDS(cvs, "vignettes/cvs.rds")


allCvsAndLlrs <- data.frame()
for (t in unique(maxSprtSimulationData$lookTime)) {

  # Compute likelihood profiles and LLR for all negative controls:
  negativeControlProfilesLookTt <- list()
  llrsLookTt <- c()
  dataLookTt <- maxSprtSimulationData[maxSprtSimulationData$lookTime == t, ]
  for (i in 1:50) {
    dataOutcomeIlookTt <- dataLookTt[dataLookTt$outcomeId == i, ]
    cyclopsData <- createCyclopsData(Surv(time, outcome) ~ exposure,
      modelType = "cox",
      data = dataOutcomeIlookTt
    )
    fit <- fitCyclopsModel(cyclopsData)
    
    # likelihood profile:
    llProfile <- getCyclopsProfileLogLikelihood(
      object = fit,
      parm = "exposureTRUE",
      bounds = log(c(0.1, 10))
    )
    negativeControlProfilesLookTt[[i]] <- llProfile
    
    # LLR:
    llNull <- getCyclopsProfileLogLikelihood(
      object = fit,
      parm = "exposureTRUE",
      x = 0
    )$value
    if (fit$return_flag == "ILLCONDITIONED" || coef(fit) < 0) {
      llr <- 0
    } else {
      llr <- fit$log_likelihood - llNull
    }
    llrsLookTt[i] <- llr
  }

  # Fit null distribution:
  nullLookTt <- fitNullNonNormalLl(negativeControlProfilesLookTt)

  # Compute calibrated and uncalibrated CV for all negative controls:
  cvs <- c()
  calibratedCvsLookT <- c()
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
    
    # Note: uncalibrated CV will be same for every t, but computing in loop 
    # over t for clarity of code:
    cv <- computeCvBinomial(
      groupSizes = outcomesPerLook,
      z = unexposedTime / exposedTime,
      minimumEvents = 1,
      alpha = 0.05
    )
    cvs[i] <- cv
    
    calibratedCv <- computeCvBinomial(
      groupSizes = outcomesPerLook,
      z = unexposedTime / exposedTime,
      minimumEvents = 1,
      alpha = 0.05,
      nullMean = nullLookTt[1],
      nullSd = nullLookTt[2]
    )
    calibratedCvsLookT[i] <- calibratedCv
  }

  # Store in a data frame:
  allCvsAndLlrs <- rbind(
    allCvsAndLlrs,
    data.frame(
      t = t,
      outcomeId = 1:50,
      llr = llrsLookTt,
      cv = cvs,
      calibrateCv = calibratedCvsLookT
    )
  )
}
saveRDS(allCvsAndLlrs, "vignettes/allCvsAndLlrs.rds")
