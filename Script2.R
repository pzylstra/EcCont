# SCRIPT 2
# Validate FRaME and Cheney et al

rm(list=ls())
gc()
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(frame)))
suppressWarnings(suppressMessages(library(berryFunctions)))
suppressWarnings(suppressMessages(library(extraDistr)))

# Import data  
datAll <- read.csv("tingleSurvey.csv") # Survey data
default.species.params <- read.csv("TraitsTingleFin.csv") # Plant traits
weather <- read.csv("7_10March2001.csv") # Weather conditions
Dynamics <- read.csv("Dynamics.csv") # Forest dynamics

nsH <- 0.5 # Height of suspended litter
density <- 300 # Mean density of woody suspended litter
REPLICATES <- 5 # Replicates for each day of weather

# Model species richness
pointRich <- rich(datAll, thres = 5, pnts = 10)

# Create pseudo-transects and format for modelling
dat <- pseudoTransect(Dynamics, pointRich, default.species.params, pnts = 100, Age = 64, 
                      growth = TRUE, thin = TRUE, prune = TRUE) %>%
  mutate(Site = ceiling(Point/10))
tabs <- frameSurvey(dat, default.species.params, pN ="Point", spName ="Species", base = "base", top = "top", he = "he", ht = "ht",
                    wid = "width", rec = "Site", sN = "SiteName", negEx = 2, max = 54.22, rate = 0.026, a = 3.35, b = 0.832, age = "Age", 
                    surf = 10, density = 300, nsH = 0.5, cover = 0.8, aQ = -0.0113, bQ = 0.9731, cQ = -1.5268, maxNS = NA, rateNS = NA)

# Model behaviour for each site
runs <- data.frame()
#IP <- data.frame()
#scorch <- data.frame()

for (r in unique(dat$Site, incomparables = FALSE)) {
  #cat("Modelling record", r)
  print(r)
  Structure <- filter(tabs[[2]], record == r)
  Flora <- filter(tabs[[1]], record == r)
  
  # Update suspNS
  T <- filter(default.species.params, name == "suspNS")
  suspNS <- as.numeric(filter(Flora, species == "suspNS")$weight[1])
  lengthS <- (0.6*((0.1 * suspNS) / (nsH * density))) / (pi * (T$leafThickness[1]/2)^2)
  sepS <- mean(sqrt(sqrt(nsH/lengthS)^2*2),sqrt(nsH/lengthS))
  default.species.params$leafSeparation[which(default.species.params$name == "suspNS")] <- sepS
  
  # Build params and model behaviour
  base.params <- suppressWarnings(buildParams(Structure = Structure, Flora = Flora, default.species.params, a = r,
                                              fLine = 100, slope = 0, temp = 30, dfmc = 0.1, wind = 10))
  
  weatherSet_Frame(base.params, weather, Structure = Structure, Flora = Flora, a = r, db.path = "out_mc.db", jitters = REPLICATES, l = 0.1,
                   Ms = 0.01, Pm = 1, Mr = 1.001, updateProgress = NULL)
  
  # MODEL CHENEY et al
  # Get fuel inputs from forest description
  S <- min(as.numeric(Flora$weight[Flora$species=="Litter"])/4,4)
  NS <- min(as.numeric(Flora$weight[Flora$species=="suspNS"]),4)
  NSh <- min(max(as.numeric(Flora$top[Flora$stratum==1]),na.rm = TRUE)*100,150)
  Eh <- min(max(as.numeric(Flora$top[Flora$stratum==2]),na.rm = TRUE)*100,500)
  
  #LOAD AND ORGANISE RESULTS, MODEL CHENEY
  res<-ffm_db_load("out_mc.db")
  runsA <- frameSummary(res$FlameSummaries, res$Sites, res$ROS, res$SurfaceResults) %>%
    mutate(S = S,
           NS = NS,
           NSh = NSh,
           Eh = Eh,
      Cheney_ROS = 30+1.5308*(max(0,wind_kph-5))^0.8576*S^0.9301*(NS*NSh)^0.6366*1.03*
             ((deadFuelMoistureProp*100)^-1.495/0.0545),
           Cheney_FH = 0.0193*Cheney_ROS^0.723*exp(0.0064*Eh)*1.07)
  
  # Collate results
  runs <- rbind(runs, runsA)
}
write.csv(runs, "Validation.csv")

# Summarise findings
F_res <- data.frame("Model" = "FRaME", "FH" = runs$fh) %>%
  mutate(diff = 10-FH)
C_res <- data.frame("Model" = "Cheney", "FH" = runs$Cheney_FH) %>%
  mutate(diff = 10-FH)
test <- t.test(F_res$diff,C_res$diff, paired = TRUE)
abs(mean(C_res$diff)) / abs(mean(F_res$diff))

F_perc <- round(length(which(F_res$FH < 10))/2,0)
F_mean <- round(mean(F_res$FH),1)
round(quantile(F_res$FH, c(0,0.25, 0.5, 0.75, 1)),1)
C_perc <- round(length(which(C_res$FH < 10))/2,0)
C_mean <- round(mean(C_res$FH),1)
round(quantile(C_res$FH, c(0,0.25, 0.5, 0.75, 1)),1)


cat(" Mean flame height (FRaME)", F_mean, "m", "\n",
    "Mean flame height (Cheney)", C_mean, "m", "\n",
    "Significance (p)", test$p.value,  "\n",
    "Mean S", round(mean(runs$S),1),"\n",
    "Mean NS", round(mean(runs$NS),1),"\n",
    "Mean NSh", round(mean(runs$NSh),0),"\n",
    "Mean Eh", round(mean(runs$Eh),0),"\n",
    "Min Eh", round(min(runs$Eh),0), "\n",
    "Max Eh", round(max(runs$Eh),0))
