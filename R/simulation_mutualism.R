DAISIE_sim_core_mutualism < function(
  time,
  M0,
  
){
  #### Initialization ####
  timeval <- 0
  totaltime <- time
  mainland_nplant <- NROW(M0)
  mainland_nanimal <- NCOL(M0)
  mainland_ntotal <- mainland_nplant + mainland_nanimal
  if(mainland_n != 0){
    mainland_spec <- seq(1, mainland_ntotal, 1)
  }else{
    mainland_spec = c()
  }
  maxspecID <- mainland_ntotal
  
  island_spec <- c()
  stt_table <- matrix(ncol = 7)
  colnames(stt_table) <- c("Time","nIp","nAp","nCp","nIa","nAa","nCa")
  stt_table[1,] <- c(totaltime,0,0,0,0,0,0)
  
  # update_rates_mutual(gam0_lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain,
  #                  M0,Mt,p_status,a_status)
  #load("rates.Rdata")
  
  num_spec <- length(island_spec[, 1])
  num_spec_plant <- length(which(island_spec[,8] == "p"))
  num_spec_animal <- length(which(island_spec[,8] == "a"))
  
  #### Start Monte Carlo iterations ####
  while (timeval < totaltime) {
    
    rates <- update_rates_mutual(
      gam0_lac0_par = gam0_lac0_par,
      K_par = K_par,
      mu_par = mu_par,
      laa_par = laa_par,
      lambda1 = lambda1,
      qloss = qloss,
      qgain = qgain,
      M0 = M0,
      Mt = Mt,
      p_status = p_status,
      a_status = a_status
    )
    
    timeval_and_dt <- calc_next_timeval(rates = rates,timeval = timeval)
    timeval <- timeval_and_dt$timeval
    
    possible_event <- possible_event(rates = rates)
    
    update_state <- update_state(timeval = timeval,
                                 totaltime = totaltime,
                                 possible_event = possible_event,
                                 
                                 maxspecID = maxspecID,
                                 mainland_spec = mainland_spec,
                                 island_spec = island_spec,
                                 stt_table = stt_table)
    
    island_spec <- updated_state$island_spec
    maxspecID <- updated_state$maxspecID
    stt_table <- updated_state$stt_table
    num_spec <- length(island_spec[, 1])
  }
}
#### Finalize STT ####





















