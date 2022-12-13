parDyn <- function(r) {
  sAge <- as.numeric(filter(tabs[[2]], record == r)$site[1])
  Structure <- filter(tabs[[2]], record == r)
  Flora <- filter(tabs[[1]], record == r)
  dbName <- paste("Age",r,".db",sep="")
  
  base.params <- suppressWarnings(buildParams(Structure = Structure, Flora = Flora, default.species.params, a = r,
                                              fLine = fLine, slope = slope, temp = 30, dfmc = 0.1, wind = 10))
  
  weatherSet_Frame(base.params, weather, Structure = Structure, Flora = Flora, a = r, db.path = dbName, jitters = REPLICATES, l = l,
                   Ms = Ms, Pm = Pm, Mr = Mr, updateProgress = NULL)
  #LOAD AND ORGANISE RESULTS
  res<-ffm_db_load(dbName)
  runsA <- frameSummaryBeta(res$FlameSummaries, res$Sites, res$ROS, res$SurfaceResults, res$IgnitionPaths)%>%
    mutate(time = ceiling(repId/REPLICATES),
           Site = r,
           Age = sAge)
  IP <- frame::repFlame(res$IgnitionPaths) %>%
    mutate(Age = sAge)
  scorch <- suppressMessages(frame::flora(runsA, IP, Param = base.params, Test = 80)) %>%
    select(!"wind_kph")
  runsA <- suppressMessages(left_join(runsA,scorch, by = "repId")) %>%
    mutate(scHeight = pmin(Height, maxH),
           quokkaOccupancy = NA,
           quokkaExclusion = NA)
  #QUOKKA IMPACTS
  for (n in 1:nrow(runsA)) {
    runsA$quokkaOccupancy[n] <- quokka(runsA$Age[n], runsA$scHeight[n])[1]
    runsA$quokkaExclusion[n] <- quokka(runsA$Age[n], runsA$scHeight[n])[2]
  }                
return(list(runsA, IP))
}


behaviourDynamics <- function(tabs, default.species.params, weather, REPLICATES = 5, l = 0.1,
                              Ms = 0.01, Pm = 1, Mr = 1.001, maxH = maxH, fLine = 100, slope = 0, freeCores = 1){
  
  cat("Modelling dynamics of fire behaviour and its impacts", "\n", "\n")
  
  # 1. Compile inputs
  r <- unique(as.numeric(tabs[[2]]$record), incomparables = FALSE)
  
  
  cat("Creating the quokka function", "\n", "\n")
  # Returns % likelihood of site occupancy, and likelihood of no quokkas
  quokka <- function(TSF, scorch) {
    # TSF effect on occupancy
    possibleOccupancy <- 0.1493*TSF^0.40659 
    La <- pnorm(0, possibleOccupancy, 0.1205)
    
    # Severity effect on exclusion
    scorchExclusion <- 0.15873 * scorch
    Lb <- 1 - pnorm(TSF, scorchExclusion, 0.5311)
    
    if (TSF < scorchExclusion) {
      occupancy <- 0
    } else {
      occupancy <- min(100,round(possibleOccupancy*100,1))}
    
    Likelihood <- round(max(La,Lb)*100,2)
    return(list(occupancy, Likelihood))
  }
  
  cat("Creating a cluster and loading the packages and inputs", "\n", "\n")
  
  # 2. Create a cluster of cores with replicated R on each
  nCores <- max(parallel::detectCores() - freeCores,1)
  cl <- parallel::makeCluster(nCores)
  parallel::clusterEvalQ(cl,
                         { library(dplyr)
                           library(berryFunctions)
                           library(frame)
                           library(data.table)
                           library(extraDistr)
                           library(splines)})
  # 4. Load the inputs
  parallel::clusterExport(cl,varlist=c('tabs', 'default.species.params', 'weather', 'REPLICATES',
                                       'l', 'Ms', 'Pm', 'Mr', 'maxH', 'fLine', 'slope', 'quokka'), environment())
  
  
  cat("Processing on", nCores, "threads", "\n", "\n")
  system.time(parT <- parallel::parLapply(cl, r, parDyn))
  parallel::stopCluster(cl)
  
  
  cat("Summarising outputs and finishing up.", "\n", "\n")
  runs <- data.frame()
  IP <- data.frame()
  for (n in 1:length(r)) {
    Na <- as.data.frame(parT[[n]][1])
    Nb <- as.data.frame(parT[[n]][2])
    runs <- rbind(runs,Na)
    IP <- rbind(IP,Nb)
  }
  return(list(runs,IP))
}
