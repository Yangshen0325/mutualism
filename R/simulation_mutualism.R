# Internal function of the DAISIE simulation
 library(reshape2)
 time <-0.3
 M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
 p <- 1
 p_status<-c(0,1,1,0)
 a_status<-c(1,0,0,0,1)
 gam0_par <- c(0.6,0.6)
 lac0_par <- c(0.3,0.3)
 K_par <- c(100,0.6,100,0.6)
 mu_par <- c(0.02,0.01,0.02,0.01)
 laa_par <- c(0.1,0.2,0.1,0.2)
 lambda1 <- 0.1
 qloss <- 0.1
 qgain <- 0.1

#  DAISIE_sim_core_mutualism(time=time,
# M0=M0,
 # p=p,
 # p_status,
 # a_status,gam0_,par,lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain)
DAISIE_sim_core_mutualism <- function(
  time,
  M0,
  p,
  p_status,
  a_status,
  gam0_lac0_par,
  K_par,
  mu_par,
  laa_par,
  lambda1,
  qloss,
  qgain){
  
  #### Initialization ####
  timeval <- 0 
  totaltime <- time
  Mt <- M0
  maxplantID <- NROW(M0)
  maxanimalID <- NCOL(M0)
  
  island_spec_plant <- matrix(ncol = 7)
  island_spec_plant[1,] <- c(2,2,0,"I",NA,NA,NA)
  island_spec_plant <- rbind(island_spec_plant, c(3,3,0,"A",NA,NA,NA))
  
  island_spec_animal <- matrix(ncol = 7)
  island_spec_animal[1,] <- c(1,1,0,"I",NA,NA,NA)
  island_spec_animal <- rbind(island_spec_animal,c(5,5,0,"A",NA,NA,NA))
  
  stt_table <- matrix(ncol = 7)
  colnames(stt_table) <- c("Time","nIp","nAp","nCp","nIa","nAa","nCa")
  stt_table[1,] <- c(totaltime,1,1,0,1,1,0)
  
  # rates <-update_rates_mutual(gam0_lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain,
  #                M0,Mt,p_status,a_status,island_spec_plant,island_spec_animal)
  #load("rates.Rdata")
  
  
  #### Start Monte Carlo iterations ####
  while (timeval < totaltime) {
    
    rates <- update_rates_mutual(
      gam0_par = gam0_par,
      lac0_par = lac0_par,
      K_par = K_par,
      mu_par = mu_par,
      laa_par = laa_par,
      lambda1 = lambda1,
      qloss = qloss,
      qgain = qgain,
      M0 = M0,
      Mt = Mt,
      p_status = p_status,
      a_status = a_status,
      island_spec_plant = island_spec_plant,
      island_spec_animal = island_spec_animal)
    
    timeval_and_dt <- calc_next_timeval(rates = rates,timeval = timeval)
    timeval <- timeval_and_dt$timeval
    # next time
    
    possible_event <- event(rates = rates)
    # next event
    
    updated_state <- update_state(timeval = timeval,
                                 totaltime = totaltime,
                                 possible_event = possible_event,
                                 p=p,
                                 Mt=Mt,
                                 p_status=p_status,
                                 a_status=a_status,
                                 maxplantID = maxplantID,
                                 maxanimalID = maxanimalID,
                                 island_spec_plant = island_spec_plant,
                                 island_spec_animal = island_spec_animal,
                                 stt_table = stt_table)
    
    Mt <- updated_state$Mt
    p_status <- updated_state$p_status
    a_status <- updated_state$a_status
    island_spec_plant <- updated_state$island_spec_plant
    island_spec_animal <- updated_state$island_spec_animal
    maxplantID <- updated_state$maxplantID
    maxanimalID <- updated_state$maxanimalID
    stt_table <- updated_state$stt_table
    #num_spec <- length(island_spec[, 1])
  }
  
  #### Finalize STT ####
  stt_table <- rbind(
    stt_table,
    c(
      0,
      stt_table[nrow(stt_table), 2],
      stt_table[nrow(stt_table), 3],
      stt_table[nrow(stt_table), 4],
      stt_table[nrow(stt_table), 5],
      stt_table[nrow(stt_table), 6],
      stt_table[nrow(stt_table), 7]
    )
  )
}

# save(rates, file="rates.Rdata")
# save(updated_state, file="updated_state.Rdata")

















