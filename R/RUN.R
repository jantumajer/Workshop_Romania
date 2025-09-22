library(dplR)
library(ggplot2)
library(reshape2)


################################
################################
#### Prepare input datasets ####
################################
################################

setwd("C:/Users/Jan/Desktop/VSLite_R_Workshop")

# 1] Load rwl file, detrend, build chronology
spruceCZ.series <- read.rwl("data/spruceCZ.rwl")
pineES.series <- read.rwl("data/pineES.rwl")
firCA.series <- read.rwl("data/firCA.rwl")

## Averaging cores from the same tree
spruceCZ.ids <- read.ids(spruceCZ.series, stc = c(5, 2, 1)) 
spruceCZ.trees <- treeMean(spruceCZ.series, spruceCZ.ids, na.rm = T)
#
pineES.ids <- read.ids(pineES.series, stc = c(3, 2, 1)) 
pineES.trees <- treeMean(pineES.series, pineES.ids, na.rm = T)
#
firCA.ids <- read.ids(firCA.series, stc = c(3, 2, 1)) 
firCA.trees <- treeMean(firCA.series, firCA.ids, na.rm = T)

## Detrending and chronology building
spruceCZ.det.series <- dplR::detrend(spruceCZ.trees, method = "Spline", nyrs = 60)
pineES.det.series <- dplR::detrend(pineES.trees, method = "Spline", nyrs = 60)
firCA.det.series <- dplR::detrend(firCA.trees, method = "Spline", nyrs = 60) 

spruceCZ.chronology <- chron.ars(spruceCZ.det.series, biweight = T) 
pineES.chronology <- chron.ars(pineES.det.series, biweight = T) 
firCA.chronology <- chron.ars(firCA.det.series, biweight = T) 

# 2] Load available climatic data
Temp.CZ <- read.delim("data/spruceCZ_T.txt", header = F, sep = "") # 1960-2019
Prec.CZ <- read.delim("data/spruceCZ_P.txt", header = F, sep = "")

Temp.ES <- read.delim("data/pineES_T.txt", header = F, sep = "") # 1950-2020
Prec.ES <- read.delim("data/pineES_P.txt", header = F, sep = "")

Temp.CA <- read.delim("data/firCA_T.txt", header = F, sep = "") # 1940-2016
Prec.CA <- read.delim("data/firCA_P.txt", header = F, sep = "")

###############################
###############################
#### Running VS-Lite model ####
###############################
###############################

###########################################
#### 0] Loading all required functions ####
###########################################

## 0a] Model sub-algorithms ##
# Calculation of partial growth rates to photoperiod
source("R/compute.gE.R") 
source("R/daylength.factor.from.lat.R") 

# Soil moisture model
source("R/leakybucket.monthly_annualTm.R") 
source("R/leakybucket.submonthly.R") 

# Ramp functions
source("R/std.ramp.R") # Original non-declining ramp functions
source("R/mod.ramp.R") # Modified increasing-stable-decreasing ramp functions

# Integration functions
source("R/integrate.orig.R") # Integration based on MINIMUM of growth rates (following Liebig's law)
source("R/integrate.multiplic.R") # Integration based on PRODUCT of growth rates (following initial TRACH model)

# Model definitions
source("R/VSLite.R")

## 0b] Functions to calibrate the model against local site chronology ##
source("R/randomization.R")
source("R/VSLite.iterative.R")

## 0c] Graphical functions
source("R/charts.R")

## 0d] OPTIONAL - detrending climatic data
# Implementation of approaches according to Ols et al. (2023): Dendrochronologia: 126094 and Brian and Nolan (2023): Frontiers in Climatology

source("R/climate.detrend.R")
# Prec.CZ.detrend <- climate.detrend(Prec.CZ, var = "prec", spline = 60, resolution = "month")
# Temp.CZ.detrend <- climate.detrend(Temp.CZ, var = "temp", spline = 60, resolution = "month")
# plot(rowMeans(Temp.CZ[,c(2:13)]), type = "l", col = "blue"); lines(rowMeans(Temp.CZ.detrend[,c(2:13)]), col = "red")
# plot(rowMeans(Prec.CZ[,c(2:13)]), type = "l", col = "blue"); lines(rowMeans(Prec.CZ.detrend[,c(2:13)]), col = "red")

################################
### 1] Calibrating the model ###
################################
# Site latitudes: CZ = 49.75, ES = 40.75, CA = 51.6
# Starting years (start of climatic data): CZ = 1960, ES = 1950, CA = 1940
# Ending years (last year of chronology or climatic data): CZ = 2016, ES = 2001, CA = 2016

random.par <- randomization(iter = 4500)
tuning <- VSLite.iterative(parameters = random.par, obs.chron = spruceCZ.chronology, chron.variant = "res",
                 syear = 1960, eyear = 2016,
                 phi = 49.75, Pinput = Prec.CZ, Tinput = Temp.CZ,
                 ramp = "modif", integration = "orig")

############################################################################################
### 2] Running the model for the set of parameters that previously produced the best fit ###
############################################################################################

simulation <- VSLite(phi = 49.75, Pinput = Prec.CZ, Tinput = Temp.CZ,
                     syear = 1960, eyear = 2016,
                     T1 = tuning$best$T1, T2 = tuning$best$T2, T3 = tuning$best$T3, T4 = tuning$best$T4,
                     M1 = tuning$best$M1, M2 = tuning$best$M2, M3 = tuning$best$M3, M4 = tuning$best$M4,
                     Acor = tuning$best$Acor, I_0 = tuning$best$I_0,
                     ramp  = tuning$ramp, integration = tuning$integration)


# Alternatively, you might skip calibration procedure from the step 1] and supply the model with fixed values of parameter
# (e.g., taken from literature, physiological measurements, ...)

# simulation <- VSLite(phi = 49.75, P = t(Prec), T = t(Temp),
#                     syear = min(as.numeric(rownames(chronology.sub))), eyear = max(as.numeric(rownames(chronology.sub))),
#                     T1 = 1, T2 = 11, T3 = 20, T4 = 30,
#                     M1 = 0.01, M2 = 0.3, M3 = 0.6, M4 = 0.85,
#                     Acor = 1, I_0 = -12,
#                     ramp  = "modif", integration = "orig")


###############################
### 3] Plotting the results ###
###############################
growth.matrix(simulation)

obsmod.chronologies(simulation, tuning)

growth.rates(simulation)

growth.rates.trends(simulation)

growth.rates.cumul(simulation)

growth.deficit(simulation)

phenology(simulation)

#### Other interesting variables calculated but not returned (for this exercise) by the VS-Lite model:
## Soil moisture dynamics calculated by functions leakybucket.monthly.R and leakybucket.submonthly.R
# Potential evapotranspiration = potEv
# Evapotranspiration = Etrans
# Groundwater drainage = G
# Surface runoff = R
# Model residuals (might store interesting environmental information - Kirdyanov et al. 2020: Ecology Letters)
observed <- tuning$obs.trw[as.numeric(rownames(tuning$obs.trw)) %in% c(simulation$syear : simulation$eyear), "res"]
plot((observed - mean(observed))/sd(observed) - t(simulation$mod.trw) ~ c(simulation$syear : simulation$eyear), type = "l")

