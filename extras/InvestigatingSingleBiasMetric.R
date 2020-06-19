mu <- 0.1
sigma <- 0.2

# Numeric solution ---------------------------------------
fun <- function(x, mu, sigma) {
   return(abs(x) * dnorm(x, mu, sigma)) 
}

integrate(fun, lower = -Inf, upper = Inf, mu = mu, sigma = sigma)

system.time(
  for (i in 1:10000)
    integrate(fun, lower = -Inf, upper = Inf, mu = mu, sigma = sigma)
)

# Closed form solution -----------------------------------
closedFormIntegral <- function(x, mu ,sigma) {
  mu * pnorm(x, mu, sigma) - 1 - sigma^2 * dnorm(x, mu, sigma)
}

closedFormIntegeralAbsolute <- function(mu, sigma) {
  closedFormIntegral(Inf, mu = mu, sigma = sigma) - 2*closedFormIntegral(0, mu = mu, sigma = sigma) + closedFormIntegral(-Inf, mu = mu, sigma = sigma)
}

closedFormIntegeralAbsolute(mu, sigma)

system.time(
  for (i in 1:10000)
    closedFormIntegeralAbsolute(mu, sigma)
)

# Coverage --------------------------------------------

grid <- expand.grid(mean = seq(0,0.5, 0.1), sd = seq(0,0.5, 0.1))



evaluateGridPoint <- function(row) {
  singleEvaluation <- function(i, row) {
    data <- simulateControls(n = 50, mean = row$mean, sd = row$sd, runif(50, min = 0.1, max = 1)) 
    null <- fitMcmcNull(data$logRr, data$seLogRr)
    return(computeExpectedSystematicError(null))
  }
  
  eval <- lapply(1:100, singleEvaluation, row = row)
  eval <- dplyr::bind_rows(eval)
  eval$mean <- row$mean
  eval$sd <- row$sd
  null <- c(mean = row$mean, sd = row$sd)
  class(null) <- "null"
  eval$trueEsr <- computeExpectedSystematicError(null)
  return(eval)
}

cluster <- ParallelLogger::makeCluster(3)
ParallelLogger::clusterRequire(cluster, "EmpiricalCalibration")
eval <- ParallelLogger::clusterApply(cluster, split(grid, 1:nrow(grid)), evaluateGridPoint)
ParallelLogger::stopCluster(cluster)

eval <- dplyr::bind_rows(eval)
saveRDS(eval, "c:/temp/esrEval.rds")
eval$coverage <- eval$trueEsr >= eval$lb95ci & eval$trueEsr < eval$lb95ub
agg <- aggregate(coverage ~ mean + sd + trueEsr, data = eval, mean)

row <- data.frame(mean = 0.05, sd = 0.05)
pointEval <- evaluateGridPoint(row)
pointEval$coverage <- pointEval$trueEsr >= pointEval$lb95ci & pointEval$trueEsr < pointEval$lb95ub
mean(pointEval$coverage)

data <- simulateControls(n = 50, mean = row$mean, sd = row$sd, runif(50, min = 0.1, max = 0.5)) 
plotCalibrationEffect(data$logRr, data$seLogRr, null = null, showCis = TRUE)
null <- fitMcmcNull(data$logRr, data$seLogRr)
null

null <- fitNull(data$logRr, data$seLogRr)
null

# Constructing a null for absolute systematic error ---------------------
data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]
null <- fitMcmcNull(negatives$logRr, negatives$seLogRr)
plotCalibrationEffect(negatives$logRr, negatives$seLogRr, showCis = TRUE, showExpectedSystematicError = TRUE)


data <- simulateControls(n = 50, mean = 0, sd = 0, runif(50, min = 0.1, max = 1))
plotCalibrationEffect(data$logRr, data$seLogRr, showCis = F, showExpectedSystematicError = TRUE)
plotCalibrationEffect(data$logRr, data$seLogRr, showCis = T, showExpectedSystematicError = TRUE)


negatives <- simulateControls(n = 25, mean = 0, sd = 0, seLogRr = runif(25, 0.1, 0.2))
null <- fitMcmcNull(negatives$logRr, negatives$seLogRr)


fitNull(negatives$logRr, negatives$seLogRr)



computeOne <- function(i) {
  data <- simulateControls(n = nrow(negatives), mean = 0, sd = 0, seLogRr = negatives$seLogRr)
  null <- fitNull(data$logRr, data$seLogRr)
  return(computeExpectedSystematicError(null))
}
system.time(
dist <- sapply(1:1000, computeOne)
)

# user  system elapsed 
4.25    0.07    4.36 

hist(dist)

null <- fitNull(negatives$logRr, negatives$seLogRr)


e <- computeExpectedSystematicError(null)
mean(dist > e)

plotCalibrationEffect(data$logRr, data$seLogRr, showCis = T, showExpectedSystematicError = TRUE)

hist(attr(null, "mcmc")$chain[, 1])
hist(attr(null, "mcmc")$chain[, 2])

mean(attr(null, "mcmc")$chain[, 2] == 0)

hist(1/sqrt(attr(null, "mcmc")$chain[, 2]))



logLikelihoodMuSigma <- function(logRr, seLogRr) {
  result <- 0
  for (i in 1:length(logRr)) {
    result <- result + dnorm(logRr[i], seLogRr[i], log = TRUE)
  }
  if (is.infinite(result)) {
    return(-99999)
  } else {
    return(result)
  }
}

negatives <- simulateControls(n = 25, mean = 00, sd = 0, seLogRr = runif(25, 0.1, 0.5))

# Cumputing absolute error likelihood
llAbsError <- function(theta, estimate, se) {
  absError <- rep(1/sqrt(theta), length(estimate))
  result <- -sum(log(dnorm(absError, estimate, se) +
    dnorm(-absError, estimate, se)))
  print(paste(theta, result))
  if (length(result) == 0 || is.infinite(result))
    result <- 99999
  result
}

fit <- optimize(llAbsError, interval = c(-10, 10), estimate = negatives$logRr, se = negatives$seLogRr)
exp(fit$minimum)

fit <- optim(0, llAbsError, estimate = negatives$logRr, se = negatives$seLogRr, hessian = T)
exp(fit$par)

x <- seq(fit$par-2, fit$par + 2, by = 0.1)
y <- sapply(x, llAbsError, estimate = negatives$logRr, se = negatives$seLogRr)
plot(x,y)

negatives <- simulateControls(n = 25, mean = 0, sd = 0, seLogRr = runif(25, 0.3, 1))
fitNull(negatives$logRr, negatives$seLogRr)
meta <- meta::metagen(negatives$logRr, negatives$seLogRr, studlab = 1:nrow(negatives), sm = "RR")
summary(meta)
meta$tau
meta$TE.random
