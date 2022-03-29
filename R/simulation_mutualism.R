
DAISIE_sim_core_mutualism <- function(totaltime,
                                      mutualism_pars){
  
  #### Initialization ####
  timeval <- 0
  totaltime <- 1 
  M0 <- mutualism_pars$M0
  maxplantID <- NROW(M0)
  maxanimalID <- NCOL(M0)
  Mt <- M0
  p_status <- rep(0, NROW(M0))
  a_status <- rep(0, NCOL(M0))
  
  island_spec <- c()
 
  stt_table <- matrix(ncol = 7)
  colnames(stt_table) <- c("Time", "nIp", "nAp", "nCp", "nIa","nAa","nCa")
  stt_table[1, ] <- c(totaltime, 0, 0, 0, 0, 0, 0)
  
  #### Start Monte Carlo iterations ####
  while (timeval < totaltime) {
    rates <- update_rates_mutualism(Mt = Mt,
                                    p_status = p_status,
                                    a_status = a_status,
                                    mutualism_pars = mutualism_pars,
                                    island_spec = island_spec)
    # next time
    timeval_and_dt <- calc_next_timeval_mutualism(rates = rates,timeval = timeval)
    timeval <- timeval_and_dt$timeval
    
    # next event
    if(timeval <= totaltime){
      rates <- update_rates_mutualism(Mt = Mt,
                                      p_status = p_status,
                                      a_status = a_status,
                                      mutualism_pars = mutualism_pars,
                                      island_spec = island_spec)
      
      possible_event <- DAISIE_sample_event_mutualism(rates = rates)
      print(possible_event$L1)
      updated_state <- DAISIE_sim_update_state_mutualism(timeval = timeval,
                                                         totaltime = totaltime,
                                                         possible_event = possible_event,
                                                         Mt = Mt,
                                                         p_status = p_status,
                                                         a_status = a_status,
                                                         maxplantID = maxplantID,
                                                         maxanimalID = maxanimalID,
                                                         island_spec = island_spec,
                                                         stt_table = stt_table,
                                                         mutualism_pars = mutualism_pars)
  
      Mt <- updated_state$Mt
      p_status <- updated_state$p_status
      a_status <- updated_state$a_status
      island_spec <- updated_state$island_spec
      maxplantID <- updated_state$maxplantID
      maxanimalID <- updated_state$maxanimalID
      stt_table <- updated_state$stt_table
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
  return(list(island_spec = island_spec,
              stt_table = stt_table))
  #island <- DAISIE_create_island_mutualism(
  #  stt_table = stt_table,
  #  totoaltime = totaltime,
  #  island_spec = island_spec,
  #  M0 = M0)
 }
#set.seed(1)
#result <- DAISIE_sim_core_mutualism(totaltime =1,mutualism_pars = mutualism_pars)

















