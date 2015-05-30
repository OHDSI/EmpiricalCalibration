.simulation <- function(){
  
  n = 5000
  mu = 0
  sigma2 = 0.5
  y_hat <- 0
  theta <- rnorm(n,mean=mu,sd=sigma2)
  se <- runif(n,min=0.000001,max=0.000001)
  y <- rnorm(n,mean=y_hat + theta,sd=se)
  data <- data.frame(LOGRR = y, SE = se)
  data$GROUND_TRUTH <- 0
  
  null <- fitNull(data$LOGRR, data$SE)
  calibrateP(1,0.25,null,TRUE)
  
  rand <- MASS::mvrnorm(1000000,c(null[1], log(null[2])),attr(null,"CovarianceMatrix"))
  P_upper_bound = pnorm((rand[,1]-positive$logRr)/sqrt(exp(rand[,2])^2+positive$seLogRr^2)) 
  P_lower_bound = pnorm((positive$logRr-rand[,1])/sqrt(exp(rand[,2])^2+positive$seLogRr^2)) 
  p <- P_upper_bound
  p[P_lower_bound < p] <- P_lower_bound[P_lower_bound < p]
  p <- p * 2
  hist(log(p),breaks=100)
}


.simulationForSusan <- function(){
  
  n = 50
  theta <- runif(n,min = log(0.4), max = log(0.8))
  se <- runif(n,min=0.001,max=1)
  y <- rnorm(n,mean=0 + theta,sd=se)
  data <- data.frame(LOGRR = y, SE = se)
  data$GROUND_TRUTH <- 0
  plotForest(data$LOGRR, data$SE,1:nrow(data))
  plotCalibration(data$LOGRR, data$SE)
}

.simulationSusan <- function(){
  getNegCtrl <- function(n, niter = 1, a = 0, b = 0) {
    res <- matrix(NA, ncol =3, nrow = niter)
    colnames(res) = c("logOR", "SE", "p-value")
    for (i in 1:niter){
      W1 <- rnorm(n, 1, 1)
      U <- rnorm(n)
      aprioriKnowledge <- runif(n, a, b)
      A <- rbinom(n, 1, plogis(.1  - .1* W1 - .2*U + 0.4 * aprioriKnowledge))
      Y <- rbinom(n, 1, plogis(-1 + .2*W1 + .1*U + 0.4 * aprioriKnowledge))
      m <- glm(Y~A + W1, family = "binomial")  # unmeasured confounding - same for all
      res[i,] <- c(logOR = coef(m)[2], SE = coef(summary(m))[2,2], pvalue = coef(summary(m))[2,4])
    }
    return(res)
  }
  
  n <- 10000
  set.seed(10)
  knownNeg.mild.train <- getNegCtrl(n, 50, 1,2)
  knownNeg.mild.test  <- getNegCtrl(n, 1000, 1, 2)
  
  #forestPlot(knownNeg.mild.train[,"logOR"], knownNeg.mild.train[,"SE"],1:nrow(knownNeg.mild.train))
  null.mild <- fitNull(knownNeg.mild.train[,"logOR"], knownNeg.mild.train[,"SE"])
  plotCalibration(knownNeg.mild.train[,"logOR"], knownNeg.mild.train[,"SE"])
  
  p.cal.mild.known <- calibrateP(knownNeg.mild.test[,"logOR"], knownNeg.mild.test[,"SE"],null.mild,pValueConfidenceInterval = TRUE)
  mean(p.cal.mild.known <= 0.05)
  mean( knownNeg.mild.test[,"p-value"] <= 0.05)
  
  
  null.mild <- fitNull(knownNeg.mild.test[,"logOR"], knownNeg.mild.test[,"SE"])
  p.cal.mild.known <- calibrateP(knownNeg.mild.train[,"logOR"], knownNeg.mild.train[,"SE"],null.mild,pValueConfidenceInterval = TRUE)
  
}
