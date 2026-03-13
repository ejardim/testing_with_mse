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
data(ple4)
data(ple4.index)
stk <- window(ple4, start=1980)
idx <- FLIndices(window(ple4.index, start=1980))
fit0 <- sca(stk, idx)

# set up
its <- 250
seed <- 123

#====================================================================
# OM conditioning based on stock assessment
# Uncertainty is introduced by simulating from multivariate normal
# Uncertainty is propagated through the SR model into reference points
#====================================================================

# fit loaded from model_a4a.R
fit <- simulate(fit0, its, seed=seed)
pred <- predict(fit)
om <- stk + fit
# repeating these to have sr and brp distributions
om.sr <- fmle(as.FLSR(om, model="bevholt"), control = list(trace = 0))
om.brp <- brp(FLBRP(om, sr=om.sr))
plot(om)

save(om, om.sr, om.brp, its, seed, fit, pred, idx, stk, file="model/om.rda")

