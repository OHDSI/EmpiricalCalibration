library(EmpiricalCalibration)
library(ggplot2)

options('fftempdir' = 's:/fftemp')

analysisSum <- readRDS("s:/temp/cohortMethodVignette2/analysisSummary.rds")


analysisSum$hoi <- analysisSum$outcomeId == 192671
a1 <- analysisSum[analysisSum$analysisId == 1, ]
a4 <- analysisSum[analysisSum$analysisId == 4, ]
a4$logRrAdj <- a4$logRr
a4$seLogRrAdj <- a4$seLogRr
a1and4 <- merge(a1[, c("outcomeId", "logRr", "seLogRr")], a4[, c("outcomeId", "logRrAdj", "seLogRrAdj")])

x <- exp(seq(log(0.25), log(10), by = 0.01))
seTheoretical <- sapply(x, FUN = function(x) {
  abs(log(x))/qnorm(0.975)
})
breaks <- c(0.25, 0.5, 1, 2, 4, 6, 8, 10)
theme <- element_text(colour = "#000000", size = 12)
themeRA <- element_text(colour = "#000000", size = 12, hjust = 1)
ggplot(a1, aes(x = exp(logRr), y = seLogRr), environment = environment()) +
  geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.5) +
  geom_vline(xintercept = 1, size = 1) +
  geom_area(aes(x = x, y = seTheoretical),
            fill = rgb(0, 0, 0),
            data = data.frame(x = x, seTheoretical = seTheoretical),
            colour = rgb(0, 0, 0, alpha = 0.1),
            alpha = 0.1) +
  geom_line(aes(x = x, y = seTheoretical),
            colour = rgb(0, 0, 0),
            linetype = "dashed",
            data = data.frame(x = x, seTheoretical = seTheoretical),
            size = 1,
            alpha = 0.5) +
  geom_point(shape = 21,
             size = 2,
             data = a1[!a1$hoi, ],
             fill = rgb(0, 0, 0, alpha = 0.5),
             colour = rgb(0, 0, 0)) +
  geom_point(shape = 21,
             size = 2,
             data = a4[!a4$hoi, ],
             fill = rgb(0, 0, 1, alpha = 0.5),
             colour = rgb(0, 0, 0.8)) +  
  geom_point(shape = 23,
             size = 4,
             data = a1[a1$hoi, ],
             fill = rgb(0, 0, 0, alpha = 0.5)) +   
  geom_point(shape = 23,
             size = 4,
             data = a4[a4$hoi, ],
             fill = rgb(1, 1, 0, alpha = 0.5)) +  
  geom_segment(data = a1and4,
               aes(x = exp(logRr),
                   y = seLogRr,
                   xend = exp(logRrAdj),
                   yend = seLogRrAdj),
               arrow = arrow(length = unit(0.01, "npc"), type='closed'),
               alpha = 0.5) +
  geom_hline(yintercept = 0) +
  scale_x_continuous("Hazard ratio",
                     trans = "log10",
                     limits = c(0.25, 10),
                     breaks = breaks,
                     labels = breaks) +
  scale_y_continuous("Standard Error", limits = c(0, 1.5)) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#FAFAFA", colour = NA),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = themeRA,
        axis.text.x = theme,
        legend.key = element_blank(),
        strip.text.x = theme,
        strip.background = element_blank(),
        legend.position = "none")
ggsave("s:/temp/biasdist.png", width = 6, height = 4.5, dpi = 400)


negCons <- analysisSum[analysisSum$analysisId == 1 & analysisSum$outcomeId != 192671, ]
hoi <-  analysisSum[analysisSum$analysisId == 1 & analysisSum$outcomeId == 192671, ]
negConsAdj <- analysisSum[analysisSum$analysisId == 4 & analysisSum$outcomeId != 192671, ]
hoiAdj <-  analysisSum[analysisSum$analysisId == 4 & analysisSum$outcomeId == 192671, ]

negConsAdj$logRrAdj <- negConsAdj$logRr
negConsAdj$seLogRrAdj <- negConsAdj$seLogRr
ncMerge <- merge(negCons[, c("outcomeId", "logRr", "seLogRr")], negConsAdj[, c("outcomeId", "logRrAdj", "seLogRrAdj")])
