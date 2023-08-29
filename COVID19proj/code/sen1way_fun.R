# This function, sen1way_analysis, is a function to summarize results of one-way
# sensitivity analysis on NPI effectiveness, threshold for NPI and adjustment time. 

sen1way_analysis <- function(name){
  
  
  calres_ls <- list()
  NPI_ls <- list()
  
  for(k in c(1:length(name))){
    
    calres_ls[[k]] <- cbind(Date= seq(set_start_date(), as.Date("2021-03-21"), by="day"), 
                            readRDS(paste0("data/sa/CalRes_", name[k] ,".rds")))
    NPI_ls[[k]] <- readRDS(paste0("data/sa/Trigger_", name[k] ,".rds"))
    
  }
  
  ## Days of triggering on & off: 
  get_days <- function(NPI){
    NPI <- as.vector(NPI)[-1]
    rl <- rle(NPI)
    on_lengths <- rl$lengths[rl$values=="80% effective contact reduction"]
    if (sum(rl$values=="80% effective contact reduction")==0){
      on <- 0
    } else{
      
      if (tail(rl$values,1) != "80% effective contact reduction") {
        on <- mean(on_lengths,na.rm = TRUE)
      } else{
        on <- mean(on_lengths[-length(on_lengths)], na.rm = TRUE)
      }
    }
    off <- mean(rl$lengths[rl$values=="50% effective contact reduction"],na.rm = TRUE)
    return(c(on, off))
  }
  
  # Triggering
  Trig <- data.frame()
  for (j in c(1:length(name))){
    DaysC<-as.data.frame(lapply(NPI_ls[[j]], get_days))
    Scn_C<- c( mean(unlist(DaysC[1,]), na.rm=TRUE), 
               quantile(unlist(DaysC[1,]), na.rm=TRUE,c(0.025, 0.975)),
               mean(unlist(DaysC[2,]), na.rm=TRUE),
               quantile(unlist(DaysC[2,]), na.rm=TRUE,c(0.025, 0.975)),
               mean(colSums(NPI_ls[[j]] == "80% effective contact reduction",na.rm = TRUE)/202),
               quantile(colSums(NPI_ls[[j]] == "80% effective contact reduction",na.rm = TRUE)/202, 
                        na.rm=TRUE,c(0.025, 0.975)))
    Trig <- rbind(Trig, Scn_C)
  }
  colnames(Trig) <- c("On_Avg", "On_lb", "On_ub", 
                      "Off_Avg", "Off_lb", "Off_ub", 
                      "pOn_Avg","pOn_lb", "pOn_ub")
  Trig <- data.frame(Name = name, Trig)
  
  # Triggering
  Hosp <- data.frame()
  for (i in c(1:length(name))){
    tmp<-c(mean(apply(calres_ls[[i]][164:365,-1], 2, max)),
           quantile(apply(calres_ls[[i]][164:365,-1], 2, max), na.rm=TRUE,c(0.025, 0.975)))
    Hosp <- rbind(Hosp, tmp)
  }
  colnames(Hosp) <- c("Prevanlent_Hospitalizations", "lb", "ub")
  
  res <- list(Trig = Trig, Hosp = Hosp)
  return(res)
}