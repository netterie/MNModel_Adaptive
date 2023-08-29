# This function, sen_analysis, is a function to summarize results of one-way
# sensitivity analysis on risk of transmissibility (beta), relative risk of 
# hospitalizations and different level of 
# fatigue contact reductions (shown in Figure 3)


sen_analysis <- function(name = "CR10"){
  
  ## Read in data
  if (name == "CR10") {
    
    for(k in c("A", "B", "C", "D", "E", "F")){
      
      assign(paste0("CalRes_", k), cbind(Date= seq(set_start_date(), as.Date("2021-03-21"), by="day"), 
                                         readRDS(paste0("data/CalRes_", k ,".rds"))))
      assign(paste0("NPI_", k), readRDS(paste0("data/Trigger_", k ,".rds")))
      
    }
    
  } else{
    for(k in c("A", "B", "C", "D", "E", "F")){
      
      assign(paste0("CalRes_", k), cbind(Date= seq(set_start_date(), as.Date("2021-03-21"), by="day"), 
                                         readRDS(paste0("data/sa/CalRes_", k , name,".rds"))))
      assign(paste0("NPI_", k), readRDS(paste0("data/sa/Trigger_", k , name,".rds")))
      
    }
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
  
  #C
  DaysC<-as.data.frame(lapply(NPI_C, get_days))
  Scn_C<- c( mean(unlist(DaysC[1,]), na.rm=TRUE), 
             quantile(unlist(DaysC[1,]), na.rm=TRUE,c(0.025, 0.975)),
             mean(unlist(DaysC[2,]), na.rm=TRUE),
             quantile(unlist(DaysC[2,]), na.rm=TRUE,c(0.025, 0.975)),
             mean(colSums(NPI_C == "80% effective contact reduction",na.rm = TRUE)/202),
             quantile(colSums(NPI_C == "80% effective contact reduction",na.rm = TRUE)/202, 
                      na.rm=TRUE,c(0.025, 0.975))
  )
  
  
  #D
  DaysD<-as.data.frame(lapply(NPI_D, get_days))
  Scn_D<- c( mean(unlist(DaysD[1,]), na.rm=TRUE), 
             quantile(unlist(DaysD[1,]), na.rm=TRUE,c(0.025, 0.975)),
             mean(unlist(DaysD[2,]), na.rm=TRUE),
             quantile(unlist(DaysD[2,]), na.rm=TRUE,c(0.025, 0.975)),
             mean(colSums(NPI_D == "80% effective contact reduction",na.rm = TRUE)/202),
             quantile(colSums(NPI_D == "80% effective contact reduction",na.rm = TRUE)/202, 
                      na.rm=TRUE,c(0.025, 0.975))
  )
  
  #E
  DaysE<-as.data.frame(lapply(NPI_E, get_days))
  Scn_E<- c( mean(unlist(DaysE[1,]), na.rm=TRUE), 
             quantile(unlist(DaysE[1,]), na.rm=TRUE,c(0.025, 0.975)),
             mean(unlist(DaysE[2,]), na.rm=TRUE),
             quantile(unlist(DaysE[2,]), na.rm=TRUE,c(0.025, 0.975)),
             mean(colSums(NPI_E == "80% effective contact reduction",na.rm = TRUE)/202),
             quantile(colSums(NPI_E == "80% effective contact reduction",na.rm = TRUE)/202, 
                      na.rm=TRUE,c(0.025, 0.975))
  )
  
  #F
  DaysF<-as.data.frame(lapply(NPI_F, get_days))
  Scn_F<- c( mean(unlist(DaysF[1,]), na.rm=TRUE), 
             quantile(unlist(DaysF[1,]), na.rm=TRUE,c(0.025, 0.975)),
             mean(unlist(DaysF[2,]), na.rm=TRUE),
             quantile(unlist(DaysF[2,]), na.rm=TRUE,c(0.025, 0.975)),
             mean(colSums(NPI_F == "80% effective contact reduction",na.rm = TRUE)/202),
             quantile(colSums(NPI_F == "80% effective contact reduction",na.rm = TRUE)/202, 
                      na.rm=TRUE,c(0.025, 0.975))
  )
  
  Table_2 <- rbind(Scn_C, Scn_D, Scn_E, Scn_F)
  colnames(Table_2) <- c("On_Avg", "On_lb", "On_ub", 
                         "Off_Avg", "Off_lb", "Off_ub", 
                         "pOn_Avg","pOn_lb", "pOn_ub")
  
  # Peak Hospitalizations
  A <- c(mean(apply(CalRes_A[164:365,-1], 2, max)),
         quantile(apply(CalRes_A[164:365,-1], 2, max), na.rm=TRUE,c(0.025, 0.975)))
  
  B <- c(mean(apply(CalRes_B[164:365,-1], 2, max)),
         quantile(apply(CalRes_B[164:365,-1], 2, max), na.rm=TRUE,c(0.025, 0.975)))
  
  C <- c(mean(apply(CalRes_C[164:365,-1], 2, max)),
         quantile(apply(CalRes_C[164:365,-1], 2, max), na.rm=TRUE,c(0.025, 0.975)))
  
  D <- c(mean(apply(CalRes_D[164:365,-1], 2, max)),
         quantile(apply(CalRes_D[164:365,-1], 2, max), na.rm=TRUE,c(0.025, 0.975)))
  
  E <- c(mean(apply(CalRes_E[164:365,-1], 2, max)),
         quantile(apply(CalRes_E[164:365,-1], 2, max), na.rm=TRUE,c(0.025, 0.975)))
  
  F <- c(mean(apply(CalRes_F[164:365,-1], 2, max)),
         quantile(apply(CalRes_F[164:365,-1], 2, max), na.rm=TRUE,c(0.025, 0.975)))
  
  Hosp <- rbind(A,B,C,D,E,F)
  res <- list(Trig = Table_2, Hosp = Hosp)
  return(res)
}