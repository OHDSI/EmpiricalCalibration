library(EmpiricalCalibration)

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

row <- data.frame(mean = 0.1, sd = 0.1)
pointEval <- evaluateGridPoint(row)
pointEval$coverage <- pointEval$trueEsr >= pointEval$lb95ci & pointEval$trueEsr < pointEval$lb95ub
mean(pointEval$coverage)

data <- simulateControls(n = 50, mean = row$mean, sd = row$sd, runif(50, min = 0.1, max = 0.5)) 
plotCalibrationEffect(data$logRr, data$seLogRr, null = null, showCis = TRUE)
null <- fitMcmcNull(data$logRr, data$seLogRr)
null

null <- fitNull(data$logRr, data$seLogRr)
null

# Likelihood curves -----------------------------------------
library(ggplot2)
closedFormIntegeralAbsolute <- EmpiricalCalibration:::closedFormIntegeralAbsolute
truth <- c(0.0, 0.1)
data <- simulateControls(n = 50, mean = truth[1], sd = truth[2], runif(50, min = 0.1, max = 1)) 
null <- fitMcmcNull(data$logRr, data$seLogRr, iter = 1000000)
chain <- attr(null, "mcmc")$chain
dist <- apply(chain, 1, function(x) closedFormIntegeralAbsolute(x[1], 1 / sqrt(x[2])))
# hist(dist, breaks = 100)

trueEsr <- closedFormIntegeralAbsolute(truth[1], truth[2])

plot <- ggplot(data.frame(dist = log(dist)), aes(x = dist)) +
  geom_histogram(bins = 100, alpha = 0.7) +
  geom_vline(xintercept = trueEsr) +
  scale_x_continuous("Absolute systematic error")
ggsave(filename = sprintf("c:/temp/Esr_%s.png", trueEsr), plot = plot)


null <- fitNull(data$logRr, data$seLogRr)
closedFormIntegeralAbsolute(null[1], null[2])


dist2 <- abs(rnorm(length(dist), trueEsr, 0.07))
ggplot(data.frame(dist = dist2), aes(x = dist)) +
  geom_histogram(bins = 100, alpha = 0.7) +
  geom_vline(xintercept = trueEsr)

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

# Computing test for correlated esr -------------------------
method1Ncs <- simulateControls(n = 50, mean = 0.1, sd = 0.1, seLogRr = runif(50, 0.1, 1))
# second method slightly more biased:
method2Ncs <- method1Ncs
method2Ncs$logRr <- method2Ncs$logRr + rnorm(nrow(method2Ncs), mean = 0.05, sd = 0.05)

null1 <- fitMcmcNull(method1Ncs$logRr, method1Ncs$seLogRr)
null2 <- fitMcmcNull(method2Ncs$logRr, method2Ncs$seLogRr)

computeExpectedSystematicError(null1)
computeExpectedSystematicError(null2)

deltaLogRr <- method2Ncs$logRr - method1Ncs$logRr
mean(deltaLogRr)
sd(deltaLogRr)

# Direct estimation (no mu or sigma) --------------------------------------
mu <- 0.5
sigma <- 0.5
ncs <- simulateControls(n = 60, mean = mu, sd = sigma, seLogRr = runif(20, 0.1, 0.1))
closedFormIntegeralAbsolute(mu, sigma)

null <- fitMcmcNull(ncs$logRr, ncs$seLogRr)
computeExpectedSystematicError(null)

computeCi(ncs)

# x <- seq(0, 1, length.out = 100)
# plot(x, -fun(x, ncs = ncs))

computeCi <- function(ncs, alpha = 0.05) {
  fun <- function(x, ncs) {
    p <- sapply(1:nrow(ncs), function(i) log(dnorm(x, ncs$logRr[i], ncs$seLogRr[i]) + dnorm(-x, ncs$logRr[i], ncs$seLogRr[i])))
    if (is.null(dim(p))) {
      p <- sum(p)
    } else {
      p <- apply(p, 1, sum)
    }
    # print(paste(x, x + p))
    return(-p)
  }
  
  fit <- nlm(fun, 0.6, ncs = ncs)
  ese <- abs(fit$estimate)
  
  x <- seq(0,4, length.out = 1000)
  y <- exp(-fun(x, ncs = ncs))
  # plot(x,y)
  d <- list(x = x, y = y)
  class(d) <- "density"
  
  ci <- HDInterval::hdi(d, credMass = 1 - alpha)
  
  result <- data.frame(ese = ese,
                       lb = ci[1],
                       ub = ci[2])
  return(result)
  # 
  # 
  # threshold <- -fit$minimum - qchisq(1 - alpha, df = 1)/2
  # 
  # precision <- 1e-07
  # 
  # # Binary search for upper bound
  # L <- ese
  # H <- 2
  # ub <- Inf
  # while (H >= L) {
  #   M <- L + (H - L)/2
  #   llM <- -fun(M, ncs = ncs)
  #   metric <- threshold - llM
  #   if (metric > precision) {
  #     H <- M
  #   } else if (-metric > precision) {
  #     L <- M
  #   } else {
  #     ub <- M
  #     break
  #   }
  #   if (M == ese) {
  #     warning("Error finding upper bound")
  #     break
  #   } else if (M == 10) {
  #     warning("Confidence interval upper bound out of range")
  #     break
  #   }
  # }
  # 
  # # Binary search for lower bound
  # if (threshold < -fun(0, ncs = ncs)) {
  #   lb <- 0
  # } else {
  #   L <- 0
  #   H <- ese
  #   lb <- -Inf
  #   while (H >= L) {
  #     M <- L + (H - L)/2
  #     llM <- -fun(M, ncs = ncs)
  #     metric <- threshold - llM
  #     if (metric > precision) {
  #       L <- M
  #     } else if (-metric > precision) {
  #       H <- M
  #     } else {
  #       lb <- M
  #       break
  #     }
  #     if (M == ese) {
  #       warning("Error finding lower bound")
  #       break
  #     } else if (M == 0) {
  #       warning("Confidence interval lower bound out of range")
  #       break
  #     }
  #   }
  # }
  # result <- data.frame(ese = ese,
  #                      lb = lb,
  #                      ub = ub)
  # return(result)
}



# Coverage of direct approach --------------------------------------------

grid <- expand.grid(mean = seq(0,0.5, 0.1), sd = seq(0,0.5, 0.1))



evaluateGridPoint <- function(row) {
  fun <- function(x, ncs) {
    
    p <- sapply(1:nrow(ncs), function(i) log(dnorm(x, ncs$logRr[i], ncs$seLogRr[i]) + dnorm(-x, ncs$logRr[i], ncs$seLogRr[i])))
    if (is.null(dim(p))) {
      p <- sum(p)
    } else {
      p <- apply(p, 1, sum)
    }
    # print(paste(x, x + p))
    return(-p)
  }
  
  computeCi <- function(ncs, alpha = 0.05) {
    fit <- nlm(fun, 0.6, ncs = ncs)
    ese <- abs(fit$estimate)
    x <- seq(0,4, length.out = 1000)
    y <- exp(-fun(x, ncs = ncs))
    d <- list(x = x, y = y)
    class(d) <- "density"
    
    ci <- HDInterval::hdi(d, credMass = 1 - alpha)
    
    result <- data.frame(ese = ese,
                         lb = ci[1],
                         ub = ci[2])
    return(result)
  }
  
  singleEvaluation <- function(i, row) {
    data <- simulateControls(n = 50, mean = row$mean, sd = row$sd, runif(50, min = 0.1, max = 1)) 
    return(computeCi(data))
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
saveRDS(eval, "c:/temp/esrDirectEval.rds")
eval$coverage <- eval$trueEsr >= eval$lb & eval$trueEsr < eval$ub
agg <- aggregate(coverage ~ mean + sd + trueEsr, data = eval, mean)



# Mean Squared Systematic Error -----------------------------------------------
# This is not a good idea, as the score will go up when sample size goes down.
mu <- 0.3
sigma <- 0.3
# Numeric solution 
fun2 <- function(x, mu, sigma) {
  return(x^2 * dnorm(x, mu, sigma)) 
}

integrate(fun2, lower = -Inf, upper = Inf, mu = mu, sigma = sigma)

# Direct estimation
ncs <- simulateControls(n = 50, mean = mu, sd = sigma, seLogRr = runif(50, 0.05, 0.05))

fun <- function(x, ncs) {
  
  p <- sapply(1:nrow(ncs), function(i) log(dnorm(x, ncs$logRr[i], ncs$seLogRr[i]) + dnorm(-x, ncs$logRr[i], ncs$seLogRr[i])))
  if (is.null(dim(p))) {
    p <- sum(p)
  } else {
    p <- apply(p, 1, sum)
  }
  # print(paste(x, x + p))
  return(-x^2 - p)
}


x <- seq(0, 1, length.out = 100)
plot(x, -fun(x, ncs = ncs))
nlm(fun, 0.6, ncs = ncs)$estimate
