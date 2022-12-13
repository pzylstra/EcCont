# QUANTIFY CONTROLS

library(dplyr)
library(frame)
library(berryFunctions)
library(data.table)
library(extraDistr)
source("tinglePar.r")

# Import data  
datAll <- read.csv("tingleSurvey.csv") # Survey data
default.species.params <- read.csv("TraitsTingleFin.csv") # Plant traits
weather <- read.csv("7_10March2001.csv") # Weather conditions
Dynamics <- read.csv("Dynamics.csv") # Forest dynamics

REPLICATES <- 5 # Replicates for each day of weather
control <- read.csv("runsSept.csv") # Control data from script 3
disturbance <- 52 #Period of disturbance with all controls in place, from script 3

# SET CONTROLS
# Set to FALSE for the control(s) being tested
thin <- TRUE      # Self-thinning
growth <- TRUE    # Growth
prune <-  TRUE    # Self-pruning
traits <- TRUE   # Plant traits
sLit <- TRUE      # Surface litter
susLit <- FALSE    # Suspended litter

#______________________________________________________________________
if (traits == TRUE) {
  plant.traits <- default.species.params
} else {
  plant.traits <- ctrlDiversity(default.species.params)
}
pointRich <- rich(datAll, thres = 5, pnts = 10)

# Create pseudo-transects
dat <- data.frame()
for (Age in seq(from = 1, to = 100, by = 1)) {
  out <- pseudoTransect(Dynamics, pointRich, plant.traits, pnts = 10, Age = Age, 
                        growth = growth, thin = thin, prune = prune)
  dat<- rbind(dat,out)
}
# Format for modelling
tabs <- frameSurvey(dat, plant.traits, pN ="Point", spName ="Species", base = "base", top = "top", he = "he", ht = "ht",
                    wid = "width", rec = "Site", sN = "SiteName", negEx = 2, max = 54.22, rate = 0.026, a = 3.35, b = 0.832, age = "Age", 
                    surf = 10, cover = 0.8, aQ = -0.0113, bQ = 0.9731, cQ = -1.5268, maxNS = NA, rateNS = NA, 
                    thin = thin, sLit  = sLit, dec = susLit)

# Model behaviour for each age
Results <- behaviourDynamics(tabs, default.species.params, weather, REPLICATES = 5, l = 0.1,
                             Ms = 0.01, Pm = 1, Mr = 1.001, maxH = 65.4, freeCores = 1)

###Analysis stats
runs <- Results[[1]]

IP <- Results[[2]]

# ANALYSIS
# Control data
conDist <- control$fh[control$Age < disturbance]
conMat <- control$fh[control$Age >= disturbance]
CdistMean <- mean(conDist)
CdistSD <- sd(conDist)
CmatMean <- mean(conMat)
CMatSD <- sd(conMat)
cFS <- round(CdistMean/CmatMean,1)

# Treatment data
x <- runs$Age
y <- runs$fh
fit <- lm(y ~ bs(x, degree=3))
pred <- as.vector(predict(fit))
fhTr <- data.frame('Mod' = pred) 
fhTr$rec <- seq(1:nrow(fhTr))
fhTr$Age <- ceiling(fhTr$rec/(REPLICATES*4))
trend <- fhTr %>%
  group_by(Age) %>%
  summarise_all(mean) %>%
  select(Age, Mod)

meanFh <- as.numeric(mean(trend$Mod))
high <- as.numeric(which(trend$Mod >= meanFh))
Mature <- max(high)+1
distStats <- runs$fh[runs$Age < Mature]
matStats <- runs$fh[runs$Age >= Mature]
distMean <- mean(distStats)
matMean <- mean(matStats)
Test <- t.test(distStats,matStats)
FS <- round(distMean/matMean,1)
p <- Test$p.value

# Glass' Delta effect size
distES <- round((distMean-CdistMean)/CdistSD,2)
matES <- round((matMean-CmatMean)/CMatSD,2)

cat("Disturbance period", max(high), "years", "\n")
cat("Disturbance fh", round(distMean,1), "m", "\n")
cat("Mature fh", round(matMean,1), "m", "\n")
cat("Feedback Strength", FS, "\n")
cat("Significance", p, "\n")
cat("Disturbed Effect Size", distES, "\n")
cat("Mature Effect Size", matES, "\n")
cat("Total Effect Size", matES-distES, "\n")
