
DAISIE_sim_core_mutualism <- function(
  time=1,
  mainland_n=1,
  pars=c(0.5,0.2,Inf,0.05,1),
  M=100,
  mutualism_pars = NULL
){
  
  
  # mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
  
  #### Initialization ####
  timeval <- 0
  totaltime <- time
  M0 <- mutualism_pars$M0
  maxspecID <- 1
  
  Mt <- M0
  
  island_spec <- c()
  stt_table <- matrix(ncol = 7)
  colnames(stt_table) <- c("Time", "nI", "nA", "nC", "nIa","nAa","nCa")
  
  stt_table[1, ] <- c(totaltime, 0, 0, 0, 0, 0, 0)

  
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

















