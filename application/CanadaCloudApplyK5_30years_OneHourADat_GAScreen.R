library(VGAM)
library(mvtnorm)
library(foreach)
library(doMC)
library(Rcpp)
library(RcppArmadillo)
library(changepointGA)

source("tools/rfunc.R")
sourceCpp("src/QuasiNewton.cpp")

################################################################################
################################################################################

TargetHour = 15
mytitles = paste0("Hour ", TargetHour)

################################################################################
################################################################################

datdir = "~/OneDrive - University of Louisiana Lafayette/data/"
dat = read.table(file=paste0(datdir, "CloudCover/FortStJohnAirport/hourly_day.txt"))
dat$Date = ISOdate(year = dat$V1, month = dat$V2, day = dat$V3, hour = dat$V4)

## filter1: 30 years data from 1965 to 1994
dat1 = dat[dat$V1>=1965 & dat$V1 <=1994,]
dat1$V5 = as.numeric(as.character(dat1$V5))

## filter2: to exclude Feb 29 for exact period T=365
premain = which(!(dat1$V2 == 2 & dat1$V3 == 29))
dat1 = dat1[premain, ]

# Computation setting
ss = 365
nCoreuse = 40
mytol    = 1e-05
stepsize = c(1,1)

dat2 = dat1[dat1$V4==TargetHour, ]
dat2$V6 = catg_agg(dat2$V5)
mm = dat2$V6
  
Ti    = length(mm)
K     = max(mm) + 1
y_log = matrix(0, nrow=Ti, ncol=K)
for (tt in 1:Ti) {y_log[tt,mm[tt]+1] = 1}

#---- Candidate changepoint set
lc     = floor(0.05*Ti)
uc     = floor(0.95*Ti)
tmp    = lc:uc
tauClc = lc:(lc+(floor(length(tmp)/nCoreuse) + 1)*nCoreuse-1)

nTasks = length(tauClc)
  
#---- Design matrix
AValue     = rep(1,Ti)
TrendValue = 1:Ti/Ti
BValue     = cos(2*pi*(1:Ti)/ss)
DValue     = sin(2*pi*(1:Ti)/ss)
DesignXH0  = cbind(AValue, BValue, DValue) ## Model 1 (No Trend)
# DesignXH0  = cbind(AValue, BValue, DValue) ## Model 2 (With Trend)


#---- Model fitting under H0
parinitialH0  = ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXH0)
resEstH0 = SeqNewtonCpp(parinitialH0, y_log, mm, DesignXH0, stepsize, 1000, mytol)
parEstH0 = resEstHa$parEst
logLH0 = LoglikCalcCpp(parEstH0, DesignXH0, y_log, mm, conv = mytol, stepNR=stepsize[2])
logLH0
## [1] -15934.64

#---- Model fitting under Ha via GA
IslandGA_param = list(
  subpopsize   = 4,
  Islandsize   = 10,
  Pcrossover   = 0.95,
  Pmutation    = 0.15,
  Pchangepoint = 0.1,
  minDist      = 2,
  mmax         = Ti/2 - 1,
  lmax         = 2 + Ti/2 - 1,
  maxMig       = 1000,
  maxgen       = 50,
  maxconv      = 100,
  option       = "cp",
  monitoring   = TRUE,
  parallel     = TRUE, ###
  nCore        = 10,
  tol          = 1e-5,
  seed         = 1234
)

IslandGA_operators = list(population = "AMOCpopulation",
                          selection = "AMOCselection",
                          crossover = "AMOCcrossover",
                          mutation = "AMOCmutation")

IslandGA.res = IslandGA(ObjFunc=GAfit, 
                        N=Ti,
                        IslandGA_param=IslandGA_param,
                        IslandGA_operators=IslandGA_operators,
                        DesignXH0=DesignXH0, 
                        y_log=y_log, 
                        mm=mm, 
                        stepsize=stepsize, 
                        mytol=mytol,
                        logLH0=logLH0)
IslandGA.res$overbestfit
IslandGA.res$overbestchrom


