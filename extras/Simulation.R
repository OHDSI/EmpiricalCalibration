set.seed(123)
negControls <- simulateControls(n = 50, mean = 0.01, sd = 0, seLogRr = runif(50, 0.05, 0.5))

# plotCalibration(negControls$logRr, negControls$seLogRr, useMcmc = FALSE)
null <- fitNull(negControls$logRr, negControls$seLogRr)
null

plotCalibrationEffect(negControls$logRr, negControls$seLogRr, null = null, showCis = FALSE, showExpectedSystematicError = TRUE, fileName = "c:/temp/plot1.png")
calibrateP(null, negControls$logRr[1], negControls$seLogRr[1])

# plotCalibration(negControls$logRr, negControls$seLogRr, useMcmc = TRUE)
null <- fitMcmcNull(negControls$logRr, negControls$seLogRr)
null
plotCalibrationEffect(negControls$logRr, negControls$seLogRr, null = null, showCis = TRUE, showExpectedSystematicError = TRUE)

calibrateP(null, negControls$logRr[1], negControls$seLogRr[1])


data <- simulateControls(n = 50 * 3, trueLogRr = log(c(1, 2, 4)))

plotTrueAndObserved(data$logRr, data$seLogRr, data$trueLogRr)


model <- fitSystematicErrorModel(data$logRr, data$seLogRr, data$trueLogRr)

cal <- calibrateConfidenceInterval(data$logRr, data$seLogRr, model)


plotTrueAndObserved(cal$logRr, cal$seLogRr, data$trueLogRr)

eval <- evaluateCiCalibration(data$logRr, data$seLogRr, data$trueLogRr)
plotCiCalibration(evaluation = eval)
plotCiCoverage(evaluation = eval, fileName = "c:/temp/plot.png")
plotErrorModel(data$logRr, data$seLogRr, data$trueLogRr, fileName = "c:/temp/plot.png")

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


