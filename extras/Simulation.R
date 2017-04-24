negControls <- simulateControls()
null <- fitNull(negControls$logRr, negControls$seLogRr)
plotCalibrationEffect(negControls$logRr, negControls$seLogRr, null = null)

data <- simulateControls(n = 50 * 3, trueLogRr = log(c(1, 2, 4)))

plotTrueAndObserved(data$logRr, data$seLogRr, data$trueLogRr)
plotCoverage(data$logRr, data$seLogRr, data$trueLogRr)

model <- fitSystematicErrorModel(data$logRr, data$seLogRr, data$trueLogRr)

cal <- calibrateConfidenceInterval(data$logRr, data$seLogRr, model)


plotTrueAndObserved(cal$logRr, cal$seLogRr, data$trueLogRr)



logRr <- data$logRr
seLogRr <- data$seLogRr
trueLogRr <- data$trueLogRr

data <- simulateControls()
data("sccs")
data <- sccs
p <- calibratePWithCiUsingMcmc(data$logRr,
                               data$seLogRr,
                               data$logRr[1],
                               data$seLogRr[1],
                               scale = c(0.05, 25),
                               iter = 10000)
mcmc <- attr(p, "mcmc")
mean(mcmc$acc)  # Acceptance rate
plot(ts(mcmc$chain[, 1]))  # Trace for the mean
plot(ts(mcmc$chain[, 2]))  # Trace for the precision (= 1/sqr(sd) )
mean(mcmc$chain[, 1])
mean(mcmc$chain[, 2])
1/(sqrt(quantile(mcmc$chain[, 2], c(0.025, 0.5, 0.975))))
p
null <- fitNull(data$logRr, data$seLogRr)
calibrateP(data$logRr[1], data$seLogRr[1], null, pValueConfidenceInterval = TRUE)
