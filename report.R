#====================================================================
# EJ (20251118)
# MSE
# Licensed under creative commons CC BY-SA 4.0
# (https://creativecommons.org/licenses/by-sa/4.0/)
#====================================================================

#====================================================================
# load libraries and data and setup run
#====================================================================
library(FLa4a)
library(FLasher)
library(FLBRP)
library(ggplotFL)

# data
load("results/mp_varsel.rda")
png("report/mp_varsel.png", 1000, 1200)
plot(FLStocks(perceivedstock.lst)) + labs(title = "Tracking of OM by MP with time varying selection pattern estimator/stock assessment model")
dev.off()

load("results/mp_sepsel.rda")
png("report/mp_sepsel.png", 1000, 1200)
plot(FLStocks(perceivedstock.lst)) + labs(title = "Tracking of OM by MP with F separable estimator/stock assessment model")
dev.off()


