# SCRIPT 1
# Vegetation dynamics

rm(list=ls())
library(dplyr)
library(berryFunctions)
library(frame)
library(ggplot2)
library(data.table)
library(patchwork)

# Import and arrange data
datAll <- read.csv("tingleSurvey.csv")

# Model trends
Dynamics <- floraDynamics(dat = datAll, thres = 0, pnts = 10, p = 0.01, bTest  = 2, cTest  = 10, maxiter = 1000,
                      Sr = 0.026, Sk = 54.22, NSa = -0.0113, NSb = 0.9731, NSc = -1.5268)
write.csv(Dynamics, "Dynamics.csv")

# Models
Dyn <- data.frame('Age' = numeric(0), 'Species' = character(0), 'Cover'= numeric(0), 'Min' = numeric(0), 'Max' = numeric(0),
                  'Top'= numeric(0), 'MinT' = numeric(0), 'MaxT' = numeric(0),
                  'Base'= numeric(0), 'MinB' = numeric(0), 'MaxB' = numeric(0),
                  'He'= numeric(0), 'MinHe' = numeric(0), 'MaxHe' = numeric(0),
                  'Ht'= numeric(0), 'MinHt' = numeric(0), 'MaxHt' = numeric(0),
                  'Width'= numeric(0), 'MinW' = numeric(0), 'MaxW' = numeric(0))
for (Species in 1:(length(Dynamics$Species)-2)) {
  time <- c(1:85)
  sp <- rep(NA, 85)
  c <- rep(NA, 85)
  l <- rep(NA, 85)
  m <- rep(NA, 85)
  t <- rep(NA, 85)
  tl <- rep(NA, 85)
  tm <- rep(NA, 85)
  b <- rep(NA, 85)
  bl <- rep(NA, 85)
  bm <- rep(NA, 85)
  for (year in time) {
    sp[year] <- Dynamics$Species[Species]
    c[year] <- round(pCover(mods = Dynamics, sp=Dynamics$Species[Species], Age = year), 1)
    l[year] <- round(max((c[year] - (1.96* (as.numeric(Dynamics$C_RSE[Species])/100) * c[year])), 0), 1)
    m[year] <- min(round(max((c[year] + (1.96* (as.numeric(Dynamics$C_RSE[Species])/100) * c[year])), 0), 1),100)
    t[year] <- pTop(mods = Dynamics, sp=Dynamics$Species[Species], Age = year)
    tl[year] <- min(max((t[year] - (1.96* (as.numeric(Dynamics$T_RSE[Species])/100) * t[year])), 0),100)
    tm[year] <- min(t[year] + (1.96* (as.numeric(Dynamics$T_RSE[Species])/100) * t[year]),100)
    b[year] <- pBase(mods = Dynamics, sp=Dynamics$Species[Species], Age = year) * t[year]
    bl[year] <- max((b[year] - (1.96* as.numeric(Dynamics$B_RSE[Species]) * b[year])), 0)
    bm[year] <- b[year] + (1.96* as.numeric(Dynamics$B_RSE[Species]) * b[year])
    

    out <- data.frame('Age' = time, 'Species' = sp, 'Cover'= c, 'Min' = l, 'Max' = m, 'Top' = t, 'MinT' = tl, 'MaxT' = tm, 
                      'Base' = b, 'MinB' = bl, 'MaxB' = bm)
  }
  Dyn <- rbind(Dyn, out)
}

# Create plots
tall <- filter(Dyn, Species == "Corymbia calophylla" | Species == "Eucalyptus diversicolor" | Species == "Eucalyptus diversicolor shade" | Species == "Eucalyptus guilfoylei" 
                | Species == "Eucalyptus jacksonii"| Species == "Eucalyptus jacksonii e" | Species == "Eucalyptus marginata")
mid <- filter(Dyn, Species == "Acacia pentadenia" | Species == "Agonis flexuosa" | Species == "Allocasuarina decussata" | Species == "Chorilaena quercifolia" 
              | Species == "Corymbia calophylla shade"| Species == "Eucalyptus guilfoylei shade" | Species == "Eucalyptus jacksonii shade" |
                Species == "Eucalyptus marginata shade" | Species == "Taxandria parviceps" | Species == "Trymalium odoratissimum")
low <- filter(Dyn, Species != "Corymbia calophylla" & Species != "Eucalyptus diversicolor" & Species != "Eucalyptus diversicolor shade" & Species != "Eucalyptus guilfoylei" 
              & Species != "Eucalyptus jacksonii"& Species != "Eucalyptus jacksonii e" & Species != "Eucalyptus marginata" & 
                Species != "Acacia pentadenia" & Species != "Agonis flexuosa" & Species != "Allocasuarina decussata" & 
                Species != "Chorilaena quercifolia" & Species != "Corymbia calophylla shade"& Species != "Eucalyptus guilfoylei shade" & 
                Species != "Eucalyptus jacksonii shade" & Species != "Eucalyptus marginata shade" & Species != "Taxandria parviceps" & 
                Species != "Trymalium odoratissimum")
doms <- filter(Dyn, Species == "Acacia pentadenia" | Species == "Agonis flexuosa" | Species == "Allocasuarina decussata"
              | Species == "Chorilaena quercifolia" | Species == "Trymalium odoratissimum")

# Trees
vars <- c("Cover"="darkolivegreen", "Height"="darkorange4", "Base" = "darkorange3")
sc <- 0.6
windows(9,3.4)
ggplot(data = tall, 
       aes(x = Age, y = Cover)) + 
  geom_line(aes(x = Age, y = Cover, colour = "Cover"), size = 0.5) +
  geom_line(aes(x = Age, y = Min, colour = "Cover"), size = 0.5, alpha=0.3) +
  geom_line(aes(x = Age, y = Max, colour = "Cover"), size = 0.5, alpha=0.3) +
  geom_ribbon(aes(ymin=Min,ymax=Max), fill = "darkolivegreen", alpha=0.1) +
  geom_line(aes(x = Age, y = Top/sc, colour = "Height"), size = 0.5) +
  geom_line(aes(x = Age, y = MaxT/sc, colour = "Height"), size = 0.5, alpha=0.3) +
  geom_line(aes(x = Age, y = Base/sc, colour = "Base"), size = 0.5) +
  geom_line(aes(x = Age, y = MinB/sc, colour = "Base"), size = 0.5, alpha=0.3) +
  geom_ribbon(aes(ymin=MinB/sc,ymax=MaxT/sc), fill = "darkorange4", alpha=0.1) +
  scale_colour_manual(name="Legend", values=vars) +
  facet_wrap(~Species, ncol = 4) +
  theme_bw() +
  theme(strip.background = element_rect(colour="black",
                                        fill="white"))+
  labs(x = "Years since fire") + 
  scale_y_continuous(name = "Percent cover",
                     sec.axis = sec_axis( trans=~.*sc, name="Plant height (m)"),
                     limits=c(0,105), breaks = seq(0,100,by = 50)) +
  theme(axis.title = element_text(face="bold", size=12, vjust=5),
        axis.text  = element_text(angle=90, size=8, colour = "black"))+
  theme(strip.text = element_text(size = 8, face = "italic"))

# Midstorey
vars <- c("Cover"="darkolivegreen", "Height"="darkorange4", "Base" = "darkorange3")
sc <- 0.2
windows(9,4.8)
ggplot(data = mid, 
       aes(x = Age, y = Cover)) + 
  geom_line(aes(x = Age, y = Cover, colour = "Cover"), size = 0.5) +
  geom_line(aes(x = Age, y = Min, colour = "Cover"), size = 0.5, alpha=0.3) +
  geom_line(aes(x = Age, y = Max, colour = "Cover"), size = 0.5, alpha=0.3) +
  geom_ribbon(aes(ymin=Min,ymax=Max), fill = "darkolivegreen", alpha=0.1) +
  geom_line(aes(x = Age, y = Top/sc, colour = "Height"), size = 0.5) +
  geom_line(aes(x = Age, y = MaxT/sc, colour = "Height"), size = 0.5, alpha=0.3) +
  geom_line(aes(x = Age, y = Base/sc, colour = "Base"), size = 0.5) +
  geom_line(aes(x = Age, y = MinB/sc, colour = "Base"), size = 0.5, alpha=0.3) +
  geom_ribbon(aes(ymin=MinB/sc,ymax=MaxT/sc), fill = "darkorange4", alpha=0.1) +
  scale_colour_manual(name="Legend", values=vars) +
  facet_wrap(~Species, ncol = 4) +
  theme_bw() +
  theme(strip.background = element_rect(colour="black",
                                        fill="white"))+
  labs(x = "Years since fire") + 
  scale_y_continuous(name = "Percent cover",
                     sec.axis = sec_axis( trans=~.*sc, name="Plant height (m)"),
                     limits=c(0,105), breaks = seq(0,100,by = 50)) +
  theme(axis.title = element_text(face="bold", size=12, vjust=5),
        axis.text  = element_text(angle=90, size=8, colour = "black"))+
  theme(strip.text = element_text(size = 8, face = "italic"))

# Low plants
vars <- c("Cover"="darkolivegreen", "Height"="darkorange4", "Base" = "darkorange3")
sc <- 0.05
windows(9,9)
ggplot(data = low, 
       aes(x = Age, y = Cover)) + 
  geom_line(aes(x = Age, y = Cover, colour = "Cover"), size = 0.5) +
  geom_line(aes(x = Age, y = Min, colour = "Cover"), size = 0.5, alpha=0.3) +
  geom_line(aes(x = Age, y = Max, colour = "Cover"), size = 0.5, alpha=0.3) +
  geom_ribbon(aes(ymin=Min,ymax=Max), fill = "darkolivegreen", alpha=0.1) +
  geom_line(aes(x = Age, y = Top/sc, colour = "Height"), size = 0.5) +
  geom_line(aes(x = Age, y = MaxT/sc, colour = "Height"), size = 0.5, alpha=0.3) +
  geom_line(aes(x = Age, y = Base/sc, colour = "Base"), size = 0.5) +
  geom_line(aes(x = Age, y = MinB/sc, colour = "Base"), size = 0.5, alpha=0.3) +
  geom_ribbon(aes(ymin=MinB/sc,ymax=MaxT/sc), fill = "darkorange4", alpha=0.1) +
  scale_colour_manual(name="Legend", values=vars) +
  facet_wrap(~Species, ncol = 4) +
  theme_bw() +
  theme(strip.background = element_rect(colour="black",
                                        fill="white"))+
  labs(x = "Years since fire") + 
  scale_y_continuous(name = "Percent cover",
                     sec.axis = sec_axis( trans=~.*sc, name="Plant height (m)"),
                     limits=c(0,105), breaks = seq(0,100,by = 50)) +
  theme(axis.title = element_text(face="bold", size=12, vjust=5),
        axis.text  = element_text(angle=90, size=8, colour = "black"))+
  theme(strip.text = element_text(size = 8, face = "italic"))

#______________________________________________________
# Graph dominants for Fig. 2

vars <- c("Cover"="darkolivegreen", "Height"="darkorange4", "Base" = "darksalmon")
sc <- 0.2
windows(7,7.5)
ggplot(data = doms, 
       aes(x = Age, y = Cover)) + 
  ggtitle("b.") + 
  geom_line(aes(x = Age, y = Cover, colour = "Cover"), size = 0.5) +
  geom_line(aes(x = Age, y = Min, colour = "Cover"), size = 0.5, alpha=0.3) +
  geom_line(aes(x = Age, y = Max, colour = "Cover"), size = 0.5, alpha=0.3) +
  geom_ribbon(aes(ymin=Min,ymax=Max), fill = "darkolivegreen", alpha=0.1) +
  geom_line(aes(x = Age, y = Top/sc, colour = "Height"), size = 0.5) +
  geom_line(aes(x = Age, y = MaxT/sc, colour = "Height"), size = 0.5, alpha=0.3) +
  geom_line(aes(x = Age, y = Base/sc, colour = "Base"), size = 0.5) +
  geom_line(aes(x = Age, y = MinB/sc, colour = "Base"), size = 0.5, alpha=0.3) +
  geom_ribbon(aes(ymin=MinB/sc,ymax=MaxT/sc), fill = "darkorange4", alpha=0.1) +
  scale_colour_manual(name="Legend", values=vars) +
  facet_wrap(~Species, ncol = 2) +
  theme_bw() +
  theme(strip.background = element_rect(colour="black",
                                        fill="white")) +
  labs(x = "Years since fire") + 
  scale_y_continuous(name = "Percent cover",
                     sec.axis = sec_axis( trans=~.*sc, name=""),
                     limits=c(0,105), breaks = seq(0,100,by = 50)) +
  theme(axis.title = element_text(size=14, vjust=5),
        axis.text  = element_text(angle=90, size=13, colour = "black"),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))+
  theme(strip.text = element_text(size = 12, face = "italic"))