# SCRIPT 3
# MODELS FLAME HEIGHT AND FIRE EFFECTS
gc()
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(frame)))
suppressWarnings(suppressMessages(library(berryFunctions)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(extraDistr)))
suppressWarnings(suppressMessages(library(splines)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(patchwork)))
source("tinglePar.r")


# Import data  
datAll <- read.csv("tingleSurvey.csv") # Survey data
default.species.params <- read.csv("TraitsTingleFin.csv") # Plant traits
weather <- read.csv("7_10March2001.csv") # Weather conditions
Dynamics <- read.csv("Dynamics.csv") # Forest dynamics
maxH <- max(datAll$top, na.rm = TRUE)

REPLICATES <- 5 # Replicates for each day of weather

print("Constructing pseudo-transects")
# Model species richness
pointRich <- rich(datAll, thres = 5, pnts = 10)

# Create pseudo-transects
dat <- data.frame()
for (Age in seq(from = 1, to = 100, by = 1)) {
  out <- pseudoTransect(Dynamics, pointRich, default.species.params, pnts = 100, Age = Age, 
                        growth = TRUE, thin = TRUE, prune = TRUE)
  dat<- rbind(dat,out)
}
# Format for modelling
tabs <- frameSurvey(dat, default.species.params, pN ="Point", spName ="Species", base = "base", top = "top", he = "he", ht = "ht",
                    wid = "width", rec = "Site", sN = "SiteName", negEx = 2, max = 54.22, rate = 0.026, a = 3.35, b = 0.832, age = "Age", 
                    surf = 10, density = 300, cover = 0.8, aQ = -0.0113, bQ = 0.9731, cQ = -1.5268, maxNS = NA, rateNS = NA, wNS = 1, thin = TRUE, sepSig = 0.1)

# Model behaviour for each age
Results <- behaviourDynamics(tabs, default.species.params, weather, REPLICATES = 5, l = 0.1,
                          Ms = 0.01, Pm = 1, Mr = 1.001, maxH = maxH, freeCores = 1)



frameSurvey <- function(dat, default.species.params, pN ="Point", spName ="Species", base = "base", top = "top", he = "he", ht = "ht",
                        wid = "width", rec = "Site", sN = "SiteName", negEx = 1, max = 54.22, rate = 0.026, a = 3.35, b = 0.832, age = NA, 
                        surf = 10, density = 300, cover = 0.8, aQ = NA, bQ = NA, cQ = NA, maxNS = NA, rateNS = NA, wNS = 1,
                        thin = TRUE, sLit  = TRUE, dec = TRUE, sepSig = 0.001) {
  
  # Find missing data
  entries <- which(is.na(dat[top]))
  if (length(entries)>0) {
    cat(" These rows were removed as they were missing top heights", "\n", entries, "\n", "\n")
    dat <- dat[-entries,] 
  }
  
  # Fill empty dimensions
  entries <- which(is.na(dat[ht])|is.na(dat[he]))
  if (length(entries)>0) {
    cat(" Empty values of ht & he were filled with top and base heights for these rows", "\n", entries, "\n", "\n")
    dat[ht][is.na(dat[ht])]<-dat[top][is.na(dat[ht])]
    dat[he][is.na(dat[he])]<-dat[base][is.na(dat[he])]
  }
  
  # Remove faulty data
  entries <- which(dat[ht]<dat[he]|dat[top]<dat[base])
  if (length(entries)>0) {
    cat(" These rows were removed as upper and lower heights conflicted", "\n", entries)
    dat <- dat[-entries,]
  }
  
  # Loop through records
  silist <- unique(dat$Site, incomparables = FALSE)
  Structure <- data.frame()
  Flora <- data.frame()
  
  for (rec in silist) {
    veg <- filter(dat, Site == rec)
    print(rec)
    # Find surface litter
    if (!is.na(age)) {
      AGE <- veg[1,age]
      
      # Control self-thinning
      (if (thin == TRUE && sLit == TRUE) {
        surf <- round(litter(negEx, max, rate, a, b, AGE),0)
      } else {
        preLit <- vector()
        for (x in 1:AGE) {
          preLit[x] <- litter(negEx, max, rate, a, b, x)
        }
        surf <- max(preLit)
      })
    }
    
    # Add suspended litter
    if (cover != 0) {
      if (thin == TRUE && dec == TRUE) {
        decline <- TRUE
      } else {
        decline <- FALSE
      }
      suspNS <- susp(default.species.params, density = density, cover = cover,
                     age = AGE, aQ = aQ, bQ = bQ, cQ = cQ, maxNS = maxNS, rate = rateNS, dec = decline)
      top <- suspNS[[1]]
      #Update tables
      if (top > 0) {
        row <- nrow(veg)
        rows <- round(cover*length(unique(veg$Point)),0)
        for (r in 1:rows) {
          veg[row+r,1] <- AGE
          veg[row+r,2] <- "suspNS"
          veg[row+r,3] <- 0
          veg[row+r,4] <- top
          veg[row+r,5] <- 0
          veg[row+r,6] <- top
          veg[row+r,7] <- wNS
          veg[row+r,8] <- AGE
          veg[row+r,9] <- rec
        }
      }
    }
    
    Struct <- buildStructure(veg, pN ="Point", spName ="Species", base = "base", top = "top", 
                             rec = "Site", sN = "SiteName", sepSig = sepSig)
    Flor <- buildFlora(veg, pN ="Point",  spName ="Species", base = "base", top = "top", he = "he", ht = "ht",
                       wid = "width", rec = "Site", sN = "SiteName", surf = surf, sepSig = sepSig)
    Structure <- rbind(Structure, Struct)
    Flora <- rbind(Flora, Flor)
  }
  
  return(list(Flora, Structure))
}











###Analysis stats
runs <- Results[[1]]
IP <- Results[[2]]
Runs <- apply(runs,2,as.character)
write.csv(Runs,"runsSept.csv")
write.csv(IP,"IPSept.csv")

x <- runs$Age
y <- runs$fh
#fit <- lm(y ~ bs(x, degree=4))
#pred <- as.vector(predict(fit))

fitControl=nls.control(maxiter=1000, tol=1e-7, minFactor = 1/999999999)
initBurr<-c(a=3,b=1)
  Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=runs,start=initBurr,trace=T, control = fitControl)
  BSum <- base::summary(Burr)
  Ba <- BSum$coefficients[1]
  Bb <- BSum$coefficients[2]
  pred <- as.vector(predict(Burr))

fhTr <- data.frame('Modelled' = pred) 
fhTr$rec <- seq(1:nrow(fhTr))
fhTr$Age <- ceiling(fhTr$rec/(REPLICATES*4))
runs <- left_join(runs, fhTr)
trend <- fhTr %>%
  group_by(Age) %>%
  summarise_all(mean)%>%
  select(Age, Modelled)

meanFh <- as.numeric(mean(trend$Modelled))
high <- as.numeric(which(trend$Modelled >= meanFh))
Young <- min(high)-1
Mature <- max(high)+1
runs$stage <- case_when(runs$Age < Mature ~ "Disturbed",
                        runs$Age >= Mature ~ "Mature")
runs$sub <- case_when(runs$Age <= Young ~ "Young",
                      runs$Age < Mature ~ "Regrowth",
                      runs$Age < (100-Young) ~ "Mature",
                      runs$Age >= (100-Young) ~ "Old")
distStats <- runs$fh[runs$Age < Mature]
matStats <- runs$fh[runs$Age >= Mature]
distMean <- mean(distStats)
matMean <- mean(matStats)
Test <- t.test(distStats,matStats)
FS <- round(distMean/matMean,1)
p <- Test$p.value

Y <- runs$fh[runs$sub == "Young"]
R <- runs$fh[runs$sub == "Regrowth"]
M <- runs$fh[runs$sub == "Mature"]
O <- runs$fh[runs$sub == "Old"]

# Indirect attack required
runs$indirect <- case_when(runs$fh >3 ~ 100,
                              TRUE ~ 0)

# Quokkas
quo <- runs %>%
  select(Age, quokkaOccupancy, quokkaExclusion) %>%
  group_by(Age) %>%
  summarise_all(mean)

windows(12,4)

ggplot(data = runs, aes(x = Age, y = fh, group = Age)) + 
  ggtitle("b.") +
  geom_boxplot(outlier.size = 0.1) +
  geom_smooth(method="lm", formula = y ~ bs(x, degree=12), color= "black",
              fill = "grey", alpha = 0.6, se = TRUE, aes(group=1)) +
  ylim(0,NA) +
  theme_bw() +
  labs(x = "Years since fire", y = "Flame height (m)" ) +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))+
  theme(strip.text = element_text(size = 12))

scPlot <- ggplot(data = runs, aes(x = Age, y = scHeight, group = Age)) + 
  ggtitle("c.") +
  geom_boxplot(outlier.size = 0.1) +
  geom_smooth(method="lm", formula = y ~ bs(x, degree=5), color= "black",
              fill = "grey", alpha = 0.6, aes(group=1)) +
  ylim(0,NA) +
  theme_bw() +
  labs(x = "Years since fire", y = "Scorch height (m)" ) +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))+
  theme(strip.text = element_text(size = 12))


quEx <- ggplot(data = runs, aes(x = Age, y = as.numeric(quokkaExclusion), group = Age)) + 
  ggtitle("d.") +
  geom_boxplot(outlier.size = 0.1) +
  ylim(0,100) +
  theme_bw() +
  labs(x = "Fire frequency", y = "Likelihood of quokka exclusion") +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))+
  theme(strip.text = element_text(size = 12))


cfLik <- ggplot(data = runs, aes(x = Age, y = b4, group = Age)) + 
  ggtitle("e.") +
  geom_boxplot(outlier.size = 0.1) +
  geom_smooth(method="lm", formula = y ~ bs(x, degree=5), color= "black",
              fill = "grey", alpha = 0.6, aes(group=1)) +
  ylim(0,100) +
  theme_bw() +
  labs(x = "Years since fire", y = "Likelihood of crown fire") +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))+
  theme(strip.text = element_text(size = 12))


indLik <- ggplot(data = runs, aes(x = Age, y = indirect, group = Age)) + 
  ggtitle("f.") +
  geom_boxplot(outlier.size = 0.1) +
  geom_smooth(method="lm", formula = y ~ bs(x, degree=5), color= "black",
              fill = "grey", alpha = 0.6, aes(group=1)) +
  ylim(0,100) +
  theme_bw() +
  labs(x = "Years since fire", y = "Likelihood that backburning required") +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))+
  theme(strip.text = element_text(size = 12))

windows(12,8)
scPlot + quEx + cfLik + indLik +
  plot_layout(ncol = 2)


# Quartiles in flame height for age sub-classes  
round(quantile(Y, c(0,0.25, 0.5, 0.75, 1)),1) # Young
round(quantile(R, c(0,0.25, 0.5, 0.75, 1)),1) # Regrowth
round(quantile(M, c(0,0.25, 0.5, 0.75, 1)),1) # Mature
round(quantile(O, c(0,0.25, 0.5, 0.75, 1)),1) # Old

cat("Young age", Young, "\n")
cat("Mature age", Mature, "\n")
cat("Old age", (100-Young), "\n")
cat("Young fh", round(mean(Y),1), "\n")
cat("Regrowth fh", round(mean(R),1), "\n")
cat("Disturbed fh", round(distMean,1), "\n")
cat("Mature fh", round(matMean,1), "\n")
cat("E-mature fh", round(mean(M),1), "\n")
cat("M-mature fh", round(mean(O),1), "\n")
cat("FS", round(FS,1), "\n")
cat("p", round(p,4), "\n")
