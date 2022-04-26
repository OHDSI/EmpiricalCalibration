# Need develop version of EmpiricalCalibration to compute one-sided calibrated p-values
# remotes::install_github("ohdsi/EmpiricalCalibration", ref = "develop")
library(EmpiricalCalibration)

computePFromLlr <- function(llr, nullRr = 0) {
  ifelse(logRr < nullRr, 1 - (1 - pchisq(2 * llr, df = 1))/2, (1 - pchisq(2 * llr, df = 1))/2)
}

computeLlrFromP <- function(p) {
  if (p > 0.5) {
    return(0)
  } else {
    return(qchisq(1 - 2 * p, df = 1) / 2)
  }
}

# Example assuming complete asymptotics -------------------------------------------------------
logRr <- log(1.5)
seLogRr <- 0.5
likelihood <- function(x) {
  dnorm(x, mean = logRr, sd = seLogRr)
}

# Uncalibrated 1-sided p-value using seLogRr
1 - pnorm(logRr/seLogRr)
# [1] 0.2087029

# Uncalibrated 1-sided p-value using LLR
llr <- log(likelihood(logRr) / likelihood(0))
computePFromLlr(llr)
# [1] 0.2087029


# Calibration when sigma = 0
null <- c(0.1, 0)
names(null) <- c("mean", "sd")
class(null) <- "null"

# Calibrated 1-sided p-value using seLogRr
calibrateP(null = null, logRr = logRr, seLogRr = seLogRr, twoSided = FALSE, upper = TRUE)
# 0.2706229 # Gold standard

# Calibrated 1-sided p-value using LLR
llrCalibrated <- log(likelihood(logRr) / likelihood(null[1]))
computePFromLlr(llrCalibrated)
# 0.2706229 # Same as gold standard


# Calibration when sigma != 0
null <- c(0.25, 0.25)
names(null) <- c("mean", "sd")
class(null) <- "null"

calibrateP(null = null, logRr = logRr, seLogRr = seLogRr, twoSided = FALSE, upper = TRUE)
# [1] 0.3904661 # Gold standard

f <- function(x) {
  llr <- log(likelihood(logRr) / likelihood(x))
  p <- computePFromLlr(llr, nullRr = x)
  p * dnorm(x, null[1], null[2])
}
calibratedP <- integrate(f = f, lower = null[1] - 10 * null[2], upper = null[1] + 10 * null[2])$value
calibratedP
# [1] 0.3904661 # Same as gold standard

calibratedLlr <- computeLlrFromP(calibratedP)

computePFromLlr(calibratedLlr)
# 0.3904661 # Same as gold standard


# Simulate data where asymptotics do not hold -------------------------------------------
library(Cyclops)
set.seed(1)
n <- 400
background <- 0.05
rr <- 2
treatment <- runif(n) < 0.05
pOutcome <- ifelse(treatment, background * rr, background)
outcome <- runif(n) < pOutcome
sum(outcome)
cyclopsData <- createCyclopsData(outcome ~ treatment, modelType = "lr")
fit <- fitCyclopsModel(cyclopsData)
logRr <- coef(fit)[2]
ci <- confint(fit, 2)
seLogRr <- (ci[3] - ci[2])/(2 * qnorm(0.975))
likelihood <- function(x) {
  exp(Cyclops::getCyclopsProfileLogLikelihood(object = fit,
                                              parm = 2,
                                              x = x,
                                              includePenalty = FALSE)$value)
}

x <- seq(-1, 2, by = 0.1)
plot(x, likelihood(x))

# Uncalibrated 1-sided p-value using seLogRr
1 - pnorm(logRr/seLogRr)
# 0.1659778   

# Uncalibrated 1-sided p-value using LLR
llr <- log(likelihood(logRr) / likelihood(0))
computePFromLlr(llr)
# 0.1729788 # Not the same because asymptotics don't hold 

# Calibration when sigma != 0
null <- c(0.25, 0.25)
names(null) <- c("mean", "sd")
class(null) <- "null"

calibrateP(null = null, logRr = logRr, seLogRr = seLogRr, twoSided = FALSE, upper = TRUE)
# 0.2264623

f <- function(x) {
  # print(x)
  llr <- log(likelihood(logRr) / likelihood(x))
  p <- computePFromLlr(llr, nullRr = x)
  p * dnorm(x, null[1], null[2])
}
calibratedP <- integrate(f = f, subdivisions = 10000, lower = null[1] - 5 * null[2], upper = null[1] + 5 * null[2])$value
calibratedP
# 0.2274269  # Not the same because asymptotics don't hold

# x <- seq(null[1] - 5 * null[2], null[1] + 5 * null[2], length.out = 1000)
# plot(x, f(x))

calibratedLlr <- computeLlrFromP(calibratedP)
calibratedLlr
# [1] 0.2792607

values <- EmpiricalCalibration:::sampleBinomialMaxLrr(groupSizes = sum(outcome),
                                                      p = mean(treatment),
                                                      minimumEvents = 1,
                                                      sampleSize = 1e6,
                                                      nullMean = 0,
                                                      nullSd = 0)
values <- values[order(-values)]
computeNewPFromLlr <- function(llr, nullRr = 0) {
  ifelse(logRr < nullRr,
         sapply(llr, function(x) mean(x > values)),
         1-sapply(llr, function(x) mean(x > values)))
}
computeNewLlrFromP <- function(p) {
    values[p*length(values)]
}
f <- function(x) {
  # print(x)
  llr <- log(likelihood(logRr) / likelihood(x))
  p <- computeNewPFromLlr(llr, nullRr = x)
  p * dnorm(x, null[1], null[2])
}
calibratedP <- integrate(f = f, subdivisions = 10000, lower = null[1] - 5 * null[2], upper = null[1] + 5 * null[2])$value
computeNewLlrFromP(calibratedP)
# [1] 0.4106933

# Simulate example using grid likelihood approximation with MCMC ----------------------------
library(Cyclops)
library(EmpiricalCalibration)
set.seed(1)
n <- 400
background <- 0.05
rr <- 2
treatment <- runif(n) < 0.05
pOutcome <- ifelse(treatment, background * rr, background)
outcome <- runif(n) < pOutcome
sum(outcome)
cyclopsData <- createCyclopsData(outcome ~ treatment, modelType = "lr")
fit <- fitCyclopsModel(cyclopsData)
likelihoodApproximation <- EvidenceSynthesis::approximateLikelihood(fit, approximation = "grid")

data <- simulateControls(n = 50, mean = 0, sd = 0, trueLogRr = 0)
null <- fitMcmcNull(data$logRr, data$seLogRr)

system.time(
  calibrateLlr(null, likelihoodApproximation) 
)
# user  system elapsed 
# 10.25    0.00   10.25

# user  system elapsed 
# 5.83    0.00    5.85 
