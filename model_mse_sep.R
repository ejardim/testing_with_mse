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
load("model/om.rda")

# set up (iterations were setup in om file)
tos <- mean(range(idx[[1]])[c("startf","endf")]) # time of survey
iar <- ac(range(idx[[1]])["min"]:range(idx[[1]])["max"]) # index age range
rng <- range(stk)
ny <- 10
yrs <- rng["maxyear"]+0:(ny-1)
af <- 1 # advice frequency
dl <- 1 # data lag
al <- 1 # advice lag
frefpt <- "f0.1"
estimator.lst <- split(yrs, yrs)
perceivedstock.lst <- split(yrs, yrs)
its <- 250
seed <- 123
tracking <- FLQuant(dimnames=list(quant=c("om.f", "mp.f", paste("mp",frefpt, sep="."), "mp.catch_advice"), year=yrs, iter=1:its))

#====================================================================
# Observation error model based on stock assessment estimates
# Uncertainty is introduced by adding observation error to catches
# and survey
# This process creates a stock object with estimation and observation
# errors up to maxyear and the observation error estimates for the
# future. All objects have range up to the end of the analyis.
#====================================================================

fit.oe <- simulate(fit, its, seed=seed, obserror=TRUE)
pred.oe <- predict(fit.oe)

# index
idx.oe <- idx + fit.oe
idx.oe <- window(idx.oe, end=rng["maxyear"]+ny)

flq0 <- index(idx.oe[[1]])
flq0[] <- pred$qmodel[[1]][,1] # q constant over time, no oe
index.q(idx.oe[[1]]) <- flq0

flq0[] <- pred$vmodel[[1]][,1] # q oe constant over time
q.oe <- flq0

# catch
stk.oe <- stk + fit.oe
stk.oe <- fwdWindow(stk.oe, end=rng["maxyear"]+ny)
c.oe <- catch.n(stk.oe)
c.oe[] <- pred$vmodel$catch[,1] # catch oe constant over time

#====================================================================
# Management procedure
#====================================================================

# estimator setup
# EWG 23-12 Final selected sub models:
# qmod  <- list(~ factor(replace(age, age>5, 5)))
fmod  <- ~ s(age, k=4) + factor(year)
# srmod <- ~ geomean(CV=0.6)
# n1mod <- defaultN1mod(stk.oe)
# vmod <- defaultVmod(stk.oe, idx.oe)

qmod  <- formula(qmodel(fit))
n1mod <- formula(n1model(fit))
vmod <- formula(vmodel(fit))

#--------------------------------------------------------------------
# Year 1, providing advice for year 2, data up to previous year (1-1)
#--------------------------------------------------------------------

ay <- ac(rng["maxyear"])
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# OEM
# Observation model generates data up to year 1-1
# Already sorted out in conditioning
idx <- window(idx.oe, end=dy)
stk <- window(stk.oe, end=dy)

# MP
# Estimator is sca with FLa4a
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk + estimator
# beverton and holt stock recruitment is fitted every year
sr0 <- fmle(as.FLSR(window(stk0, start=2000), model="geomean"), control = list(trace = 0))
# HCR is fmax as target, also estimated every year
ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]

# Projection
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
# catch advice is based on the projected catch if the target is applied
catch_advice <- catch(stk_fwd)[,py]

# IEM
# No implementation model

# Update tracking matrix and lists
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 2, providing advice for year 3, data up to previous year (2-1)
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
# Already sorted out in conditioning
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=dy)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk + estimator
#sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,py]

# IEM

# Update tracking matrix
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 3, providing advice for year 4, data up to previous year (3-1)
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=dy)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk3 <- stk + estimator
#sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 4, providing advice for year 5, data up to previous year (4-1)
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=dy)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk3 <- stk + estimator
#sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 5, providing advice for year 6, data up to previous year (5-1)
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=dy)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk3 <- stk + estimator
#sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 6, providing advice for year 7, data up to previous year (6-1)
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=dy)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk3 <- stk + estimator
#sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 7, providing advice for year 8, data up to previous year (7-1)
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=dy)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk3 <- stk + estimator
#sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 8, providing advice for year 9, data up to previous year (8-1)
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=dy)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk3 <- stk + estimator
#sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 9, providing advice for year 10, data up to previous year (9-1)
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=dy)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fit="MP", qmodel=qmod, vmodel=vmod, n1model=n1mod, fmodel=fmod)
stk0 <- stk3 <- stk + estimator
#sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c((as.numeric(ay)-1):as.numeric(py))
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ay] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ay] <- ftrg
tracking["mp.catch_advice",ay] <- catch_advice
estimator.lst[[ay]] <- fit
perceivedstock.lst[[ay]] <- stk0

#--------------------------------------------------------------------
# Year 10
#--------------------------------------------------------------------

ay <- ac(as.numeric(ay) + 1)
dy <- ac(as.numeric(ay) - dl)
py <- ac(as.numeric(ay) + al)
# store om value for later analysis
tracking["om.f",ay] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

save(estimator.lst, perceivedstock.lst, tracking, om, om.brp, om.sr, file="model/ara_mp_bc.rda")

#====================================================================
# Analysis of results
#====================================================================

#--------------------------------------------------------------------
# kobe is your friend
#--------------------------------------------------------------------
# ssb and F relative to MSY reference points
ssb_ssbmsy <- ssb(om)/refpts(om.brp)["msy", "ssb"]
f_fmsy <- fbar(om)/refpts(om.brp)["msy", "harvest"]

ssb_ssbmsy <- window(iterMedians(ssb(om)/refpts(om.brp)["msy", "ssb"]), start=2020)
f_fmsy <- window(iterMedians(fbar(om)/refpts(om.brp)["msy", "harvest"]), start=2020)

# code to compute kobe plot
ggplot(mapping=aes(y=c(f_fmsy), x=c(ssb_ssbmsy))) +
  # Add quadrants
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = 0, ymax = 1), fill = "green", alpha = 0.5) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "yellow", alpha = 0.5) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 1, ymax = Inf), fill = "red", alpha = 0.5) +
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = 1, ymax = Inf), fill = "yellow", alpha = 0.5) +
  # Reference lines
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # Points
  geom_point(size = 2) +
  # lines
  geom_path(arrow = arrow(type = "open", length = unit(0.15, "inches")), linewidth = 0.5) +
  # Labels and theme
  labs(
    x = expression(B / B[MSY]),
    y = expression(F / F[MSY]),
  ) +
  theme_minimal()

#--------------------------------------------------------------------
# Did the assessment track the OM?
#--------------------------------------------------------------------

plot(log(tracking["mp.f"])~log(tracking["om.f"]), ylab="f estimated in the mp", xlab="F in the om", pch=19)

perceivedstock.lst[[10]] <- om
names(perceivedstock.lst) <- c(1:9, "om")
plot(FLStocks(perceivedstock.lst))

#--------------------------------------------------------------------
# Were the reference points stable?
#--------------------------------------------------------------------

bwplot(data~year, data=tracking[paste("mp",frefpt, sep=".")], ylab=frefpt, xlab="Year")

#--------------------------------------------------------------------
# How was catch advice in relation to the status of the stock?
# (note there's no evaluation of status in the HCR)
#--------------------------------------------------------------------

plot(tracking["mp.catch_advice"]~I(tracking["mp.f"]/tracking[paste("mp",frefpt, sep=".")]), ylab="catch advice", xlab="F status in the mp", pch=19)

#--------------------------------------------------------------------
# Economic performance: C/F This ratio essentially represents the
# Efficiency of Extraction. If C/F is declining, it means the industry is
# working harder (higher F) to get the same or less output (C), signaling # a decline in the underlying natural capital (biomass) and likely a
# squeeze on profit margins. It can be seen as Return on Natural Capital
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# kobe-like C/F ~ F
#--------------------------------------------------------------------
# ssb and F relative to MSY reference points
# ssb and F relative to MSY reference points
f_fmsy <- iterMedians(fbar(om)/refpts(om.brp)["msy", "harvest"])
ssb_ssbmsy <- iterMedians(ssb(om)/refpts(om.brp)["msy", "ssb"])
c_f_msy <- refpts(om.brp)["msy", "yield"]/refpts(om.brp)["msy", "harvest"]
c_f <- catch(om)/fbar(om)
c_f_msy <- iterMedians(c_f/c_f_msy)

# code to compute kobe plot
ggplot(mapping=aes(y=c(c_f_msy), x=c(f_fmsy))) +
  # Add quadrants
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = 0, ymax = 1), fill = "red", alpha = 0.5) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "yellow", alpha = 0.5) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 1, ymax = Inf), fill = "green", alpha = 0.5) +
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = 1, ymax = Inf), fill = "yellow", alpha = 0.5) +
  # Reference lines
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # Points
  geom_point(size = 2) +
  # lines
  geom_path(arrow = arrow(type = "open", length = unit(0.15, "inches")), linewidth = 0.5) +
  # Labels and theme
  labs(
    x = expression(F / F[MSY]),
    y = c("Return on Natural Capital (C/F standardized to MSY level)"),
  ) +
  theme_minimal()

save(om, perceivedstock.lst, tracking, q.oe, c.oe, file="results/mp_sepsel.rda")


#====================================================================
# Extensions
#====================================================================

#- Use a different stock assessment model
#- Use a different assessment frequency
#- Don't reestimate reference points every year
#- Use a different catchability for the survey
#- Introduce uncertainty in catchability, or recruitment
#- ...

