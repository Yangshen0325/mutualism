
DAISIE_sim_core_mutualism <- function(time=1,
                                      mutualism_pars){
  
  #### Initialization ####
  timeval <- 0
  totaltime <- time
  M0 <- mutualism_pars$M0
  Mt <- M0
  p_status <- rep(0, NROW(M0))
  a_status <- rep(0, NCOL(M0))
  
  island_spec_plant <- c()
  island_spec_animal <- c()
  
  stt_table <- matrix(ncol = 7)
  colnames(stt_table) <- c("Time", "nIp", "nAp", "nCp", "nIa","nAa","nCa")
  stt_table[1, ] <- c(totaltime, 0, 0, 0, 0, 0, 0)
  
  
  #### Start Monte Carlo iterations ####
  while (timeval < totaltime) {
    rates <- update_rates_mutualism(Mt = Mt,
                                    p_status = p_status,
                                    a_status = a_status,
                                    mutualism_pars = mutualism_pars,
                                    island_spec_plant = island_spec_plant,
                                    island_spec_animal = island_spec_animal)
    
    timeval_and_dt <- calc_next_timeval(rates = rates,timeval = timeval)
    timeval <- timeval_and_dt$timeval
    # next time
    if(timeval <= totaltime){
      rates <- update_rates_mutualism(Mt = Mt,
                                      p_status = p_status,
                                      a_status = a_status,
                                      mutualism_pars = mutualism_pars)
      
      possible_event <- DAISIE_sample_event_mutualism(rates = rates)
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


















