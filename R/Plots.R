# @file Plots.R
#
# Copyright 2014 Observational Health Data Sciences and Informatics
#
# This file is part of EmpiricalCalibration.R
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Create a forest plot
#'
#' @description
#' \code{forestPlot} creates a forest plot of effect size estimates.
#'
#' @details
#' Creates a forest plot of effect size estimates (ratios). Estimates that are significantly
#' different from 1 (alpha = 0.05) are marked in orange, others are marked in blue.
#' 
#' 
#' @param logRr     A numeric vector of effect estimates on the log scale
#' @param seLogRr  	The standard error of the log of the effect estimates. Hint: often the standard 
#' error = (log(<lower bound 95 percent confidence interval>) - log(<effect estimate>))/qnorm(0.025) 
#' @param names		  A vector containing the names of the drugs or outcomes
#' @param xLabel  	The label on the x-axis: the name of the effect estimate
#' 
#' @return A Ggplot object. Use the \code{ggsave} function to save to file.
#' 
#' @examples 
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0,]
#' forestPlot(negatives$logRr,negatives$seLogRr, negatives$drugName)
#' 
#' @export
forestPlot <- function(logRr,seLogRr,names,xLabel="Relative risk"){
  breaks <- c(0.25,0.5,1,2,4,6,8,10) 
  theme <- ggplot2::element_text(colour="#000000", size=6) 
  themeRA <- ggplot2::element_text(colour="#000000", size=5,hjust=1) 
  themeLA <- ggplot2::element_text(colour="#000000", size=10,hjust=0) 
  col <- c(rgb(0,0,0.8,alpha=1),rgb(0.8,0.4,0,alpha=1))
  colFill <- c(rgb(0,0,1,alpha=0.5),rgb(1,0.4,0,alpha=0.5))
  data <- data.frame(DRUG_NAME = as.factor(names), 
                     LOGRR = logRr, 
                     LOGLB95RR = logRr + qnorm(0.025)*seLogRr,
                     LOGUB95RR = logRr + qnorm(0.975)*seLogRr)
  data$SIGNIFICANT <- data$LOGLB95RR > 0 | data$LOGUB95RR < 0 
  data$DRUG_NAME<-factor(data$DRUG_NAME, levels=rev(levels(data$DRUG_NAME))) 
  
  ggplot2::ggplot(data, aes(x= DRUG_NAME , y=exp(LOGRR), ymin=exp(LOGLB95RR), ymax=exp(LOGUB95RR),colour=SIGNIFICANT, fill=SIGNIFICANT), environment=environment()) +
    geom_hline(yintercept=breaks, colour ="#AAAAAA", lty=1, lw=0.2) +
    geom_hline(yintercept=1, lw=0.5) + 
    geom_pointrange(shape=23) +
    scale_colour_manual(values=col) +
    scale_fill_manual(values=colFill) +
    coord_flip(ylim=c(0.25,10)) + 
    scale_y_continuous(xLabel,trans="log10", breaks=breaks, labels = breaks) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background= element_rect(fill="#FAFAFA", colour = NA),
      panel.grid.major = element_line(colour = "#EEEEEE"),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(), 
      axis.title.x = element_blank(), 
      axis.text.y = themeRA,
      axis.text.x = theme,
      legend.key= element_blank(),
      strip.text.x = theme,
      strip.background = element_blank(),
      legend.position = "none"
    ) 	
}

logRRtoSE <- function(logRR,p,null) {
  sapply(logRR, FUN = function(logRR){
    precision <- 0.001
    if (calibrateP(logRR,precision,null) > p)
      return(0)
    L <- 0
    H <- 100
    while (H>=L) { 
      M=L+(H-L) / 2;
      if	(calibrateP(logRR,M,null) - p > precision)
        H<-M 
      else if (p - calibrateP(logRR,M,null)  > precision)
        L<-M 
      else 
        return(M)
    }
    return(L-1)
  })
} 

#' Plot the effect of the calibration
#'
#' @description
#' \code{plotCalibrationEffect} creates a plot showing the effect of the calibration.
#'
#' @details
#' Creates a plot with the effect estimate on the x-axis and the standard error on the y-axis. Negative
#' controls are shown as blue dots, positive controls as yellow diamonds. The area below the dashed line
#' indicated estimates with p < 0.05. The orange area indicates estimates with calibrated p < 0.05.
#' 
#' @param logRrNegatives     A numeric vector of effect estimates of the negative controls on the log scale
#' @param seLogRrNegatives   The standard error of the log of the effect estimates of the negative controls
#' @param logRrPositives     A numeric vector of effect estimates of the positive controls on the log scale
#' @param seLogRrPositives   The standard error of the log of the effect estimates of the positive controls
#' @param xLabel    The label on the x-axis: the name of the effect estimate
#' 
#' @return A Ggplot object. Use the \code{ggsave} function to save to file.
#' 
#' @examples 
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0,]
#' positive <- sccs[sccs$groundTruth == 1,]
#' null <- fitNull(negatives$logRr,negatives$seLogRr)
#' plotCalibrationEffect(negatives$logRr,negatives$seLogRr,positive$logRr,positive$seLogRr,null)
#' 
#' @export
plotCalibrationEffect <- function(logRrNegatives, seLogRrNegatives, logRrPositives, seLogRrPositives,xLabel="Relative risk"){
  x <- exp(seq(log(0.25),log(10),by=0.01))
  y <- logRRtoSE(log(x),0.05,null)
  seTheoretical <- sapply(x,FUN=function(x){abs(log(x))/qnorm(0.975)})
  breaks <- c(0.25,0.5,1,2,4,6,8,10) 
  theme <- ggplot2::element_text(colour="#000000", size=12) 
  themeRA <- ggplot2::element_text(colour="#000000", size=12,hjust=1) 
  themeLA <- ggplot2::element_text(colour="#000000", size=12,hjust=0) 
  ggplot2::ggplot(data.frame(x,y,seTheoretical),aes(x=x,y=y), environment=environment())+
    geom_vline(xintercept=breaks, colour ="#AAAAAA", lty=1, lw=0.5) +
    geom_vline(xintercept=1, lw=1) + 			
    geom_area(fill=rgb(1,0.5,0,alpha = 0.5),color=rgb(1,0.5,0),size=1,alpha=0.5) +
    geom_area(aes(y=seTheoretical),fill=rgb(0,0,0),colour=rgb(0,0,0,alpha=0.1),alpha = 0.1)+
    geom_line(aes(y=seTheoretical),colour=rgb(0,0,0), linetype="dashed", size=1,alpha=0.5)+
    geom_point(shape=23,aes(x,y),data=data.frame(x=exp(logRrPositives), y=seLogRrPositives), size=4,fill=rgb(1,1,0),alpha=0.8) +
    geom_point(shape=21,aes(x,y),data=data.frame(x=exp(logRrNegatives), y=seLogRrNegatives), size=2,fill=rgb(0,0,1,alpha=0.5),colour=rgb(0,0,0.8)) +
    geom_hline(yintercept=0) +
    scale_x_continuous(xLabel,trans="log10",limits = c(0.25,10), breaks=breaks,labels=breaks) + 
    scale_y_continuous("Standard Error",limits = c(0,1.5)) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background= element_rect(fill="#FAFAFA", colour = NA),
      panel.grid.major= element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      legend.key= element_blank(),
      strip.text.x = theme,
      strip.background = element_blank(),
      legend.position = "none"
    ) 
}

#' Create a calibration plot
#'
#' @description
#' \code{plotCalibration} creates a plot showing the calibration of our calibration procedure
#'
#' @details
#' Creates a calibration plot showing the number of effects with p < alpha for every level of alpha. 
#' The empirical calibration is performed using a leave-one-out design: The p-value of an effect is computed
#' by fitting a null using all other negative controls. Ideally, the calibration line should approximate the 
#' diagonal. The plot shows both theoretical (traditional) and empirically calibrated p-values.
#' 
#' @param logRr     A numeric vector of effect estimates on the log scale
#' @param seLogRr    The standard error of the log of the effect estimates. Hint: often the standard 
#' error = (log(<lower bound 95 percent confidence interval>) - log(<effect estimate>))/qnorm(0.025) 
#' 
#' @return A Ggplot object. Use the \code{ggsave} function to save to file.
#' 
#' @examples 
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0,]
#' plotCalibration(negatives$logRr,negatives$seLogRr)
#' 
#' @export
plotCalibration <- function(logRr,seLogRr){
  data <- data.frame(LOGRR = logRr, SE = seLogRr)
  data$Z <- data$LOGRR / data$SE
  data$P <- 2*pmin(pnorm(data$Z),1-pnorm(data$Z)) # 2-sided p-value   
  data$Y <- sapply(data$P,FUN <- function(x){sum(data$P < x)/nrow(data)})
  
  data$calibratedP <- vector(length=nrow(data))
  for (i in 1:nrow(data)){
    dataLeaveOneOut <- data[seq(1,nrow(data))!=i,]
    null <- fitNull(dataLeaveOneOut$LOGRR, dataLeaveOneOut$SE)
    data$calibratedP[i] = calibrateP(data$LOGRR[i],data$SE[i],null)
  }
  
  data$AdjustedY <- sapply(data$calibratedP,FUN <- function(x){sum(data$calibratedP < x)/nrow(data)})
  
  catData <- data.frame(x = c(data$P,data$calibratedP),
                        y = c(data$Y,data$AdjustedY), 
                        label = factor(c(rep("Theoretical",times=nrow(data)),rep("Empirical",times=nrow(data)))))
  catData$label <- factor(catData$label, levels = c("Empirical", "Theoretical"))
  
  names(catData) <- c("x","y","P-value calculation")
  
  breaks <- c(0,0.25,0.5,0.75,1) 
  theme <- ggplot2::element_text(colour="#000000", size=10) 
  themeRA <- ggplot2::element_text(colour="#000000", size=10,hjust=1) 
  themeLA <- ggplot2::element_text(colour="#000000", size=10,hjust=0) 
  ggplot2::ggplot(catData, aes(x=x,y=y,colour=`P-value calculation`,linetype=`P-value calculation`), environment=environment()) +
    geom_vline(xintercept=breaks, colour ="#AAAAAA", lty=1, lw=0.3) +
    geom_vline(xintercept=0.05, colour ="#888888", linetype="dashed", lw=1) +
    geom_hline(yintercept=breaks, colour ="#AAAAAA", lty=1, lw=0.3) +
    geom_abline(colour ="#AAAAAA", lty=1, lw=0.3) + 
    geom_step(direction="hv", size=1) +
    scale_colour_manual(values=c(rgb(0,0,0),rgb(0,0,0),rgb(0.5,0.5,0.5))) +
    scale_linetype_manual(values=c("solid","twodash")) +
    scale_x_continuous("Alpha",limits = c(0,1), breaks=c(breaks,0.05),labels=c("",".25",".50",".75","1",".05")) + 
    scale_y_continuous("Fraction with p < alpha",limits = c(0,1), breaks=breaks,labels=c("0",".25",".50",".75","1")) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background= element_rect(fill="#FAFAFA", colour = NA),
      panel.grid.major= element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = themeRA,
      axis.text.x = theme,
      strip.text.x = theme,
      strip.background = element_blank(),
      legend.position = "right"
    )
  
}


