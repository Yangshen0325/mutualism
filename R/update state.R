# update state of island given sampled event

# [1]: immigration event with plant species
# [2]: extinction event with plant species
# [3]: cladogenesis event with plant species
# [4]: anagenesis event with plant species
# [5]: immigration event with animal species
# [6]: extinction event with animal species
# [7]: cladogenesis event with animal species
# [8]: anagenesis event with animal species
# [9]: cospeciation event between pairs
# [10]: gain links event between pairs
# [11]: loss links event between pairs

#Update state of island given sampled event
DAISIE_sim_update_state_mutualism <- function(timeval,
                                              totaltime,
                                              possible_event,
                                              Mt,
                                              p_status,
                                              a_status,
                                              maxplantID,
                                              maxanimalID,
                                              island_spec,
                                              stt_table,
                                              mutualism_pars ){

  
  ##########################################
  # [1]: immigration event with plant species 
  if(possible_event$L1 == 1){
    
    colonist = possible_event$Var1
    p_status[colonist] <- 1
    
    if (length(island_spec[,1]) != 0)
    {
      isitthere <- which(island_spec[,1] == colonist)
    }else{
      isitthere <- c()
    }
    
    if(length(isitthere) == 0)
    {
      island_spec <- rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA,"plant"))
    }
    if (length(isitthere) != 0) {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA,"plant")
    }
  } 
  
  ##########################################
  # [2]: extinction event with plant species
  if(possible_event$L1 == 2){
    
    extinct =  possible_event$Var1
    p_status[extinct] <- 0
    
    ind <- which(island_spec[,1]==extinct)
    typeofspecies = island_spec[ind,4] 
    
    if(typeofspecies == "I")
    {
      island_spec = island_spec[-ind, ] #remove immigrant
    }
    
    if(typeofspecies == "A")
    {
      island_spec = island_spec[-ind, ] #remove anagenetic
    }
    
    if(typeofspecies == "C")##############################
    {
      #remove cladogenetic
      #first find species with same ancestor AND arrival totaltime
      sisters = intersect(which(island_spec[,2] == island_spec[ind,2]),
                          which(island_spec[,3] == island_spec[ind,3]))
      survivors = sisters[which(sisters != ind)]
      
      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic
        island_spec[survivors,4] = "A"
        island_spec[survivors,c(5,6)] = c(NA,NA)
        island_spec[survivors,7] = "Clado_extinct"
        island_spec = island_spec[-ind, ]
      }
      
      if(length(sisters) >= 3)
      {
        numberofsplits = nchar(island_spec[ind,5])
        mostrecentspl = substring(island_spec[ind,5],numberofsplits)
        
        if(mostrecentspl=="B")
        {
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }
        
        motiftofind = paste(substring(island_spec[ind,5],1,numberofsplits-1),
                            sistermostrecentspl,sep = "")
        possiblesister = survivors[which(substring(island_spec[survivors,5],1,
                                                   numberofsplits) == motiftofind)]
        
        if(mostrecentspl == "A")
        {
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec[possiblesister,6] == 
                                            min(as.numeric(island_spec[possiblesister,6])))]
          island_spec[tochange,6] = island_spec[ind, 6]
        }
        #remove the offending A/B from these species
        island_spec[possiblesister,5] = paste(substring(island_spec[possiblesister,5],1,numberofsplits - 1),
                                                    substring(island_spec[possiblesister,5],numberofsplits + 1,
                                                              nchar(island_spec[possiblesister,5])),sep = "")
        island_spec = island_spec[-ind, ]
      }
    }
    island_spec = rbind(island_spec)
  }
  
  ##########################################
  # [3]: cladogenesis event with plant species
  if(possible_event$L1 == 3){
    
    tosplit <- possible_event$Var1
    ind <- which(island_spec[,1]==tosplit)
    
    Mt <- new_Mt_clado_p(Mt=Mt, possible_event=possible_event, mutualism_pars = mutualism_pars)
    p_status[tosplit] <- 0
    p_status <- c(p_status,1,1)
    
    #if the species that speciates is cladogenetic
    if(island_spec[ind,4] == "C")
    {
      #for daughter A
      island_spec[ind, 1] <- maxplantID + 1
      oldstatus <- island_spec[ind, 5]
      island_spec[ind, 5] = paste(oldstatus,"A",sep = "")
      island_spec[ind, 6] = timeval
      island_spec[ind, 7] = NA
      island_spec[ind,8] = "plant"
      
      #for daughter B
      island_spec <- rbind(island_spec,c(maxplantID + 2,island_spec[ind, 2],island_spec[ind, 3],
                                                     "C",paste(oldstatus,"B",sep = ""),timeval,NA,"plant"))
      maxplantID = maxplantID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      island_spec[ind, 4] = "C"
      island_spec[ind, 1] = maxplantID + 1
      island_spec[ind, 5] = "A"
      island_spec[ind, 6] = island_spec[ind, 3]
      island_spec[ind, 7] = NA
      island_spec[ind,8] = "plant"
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxplantID + 2,island_spec[ind, 2],
                                                    island_spec[ind,3],"C","B",timeval,NA,"plant"))
      maxplantID = maxplantID + 2
    }
  } 
  
  ##########################################
  # [4]: anagenesis event with plant species 
  if(possible_event$L1 == 4){
    
    anagenesis = possible_event$Var1
    
    Mt <- new_Mt_ana_p(Mt=Mt, possible_event=possible_event,mutualism_pars=mutualism_pars)
    p_status[anagenesis] <- 0
    p_status <- c(p_status,1)
    
    ind <- which(island_spec[,1]==anagenesis)
    island_spec[ind, 4] = "A"
    island_spec[ind, 1] = maxplantID + 1
    island_spec[ind, 7] = "Immig_parent"
    island_spec[ind,8] = "plant"
    
    maxplantID <- maxplantID + 1
  } 
  
  ##########################################
  # [5]: immigration event with animal species
  if(possible_event$L1 == 5){
    
    colonistID <- possible_event$Var1
    a_status[colonistID] <- 1
    colonist <- colonistID + length(p_status)
    
    if (length(island_spec[,1]) != 0)
    {
      isitthere <- which(island_spec[,1] == colonist)
    }else{
      isitthere <- c()
    }
    if(length(isitthere) == 0) {
      island_spec <- rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA,"animal"))
    }
    if(length(isitthere) != 0) {
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA,"animal")
    }
  }
  
  ##########################################
  # [6]: extinction event with animal species
  if(possible_event$L1 == 6){
    
    extinctID <- possible_event$Var1
    a_status[extinctID] <- 0
    extinct <-  extinctID + length(p_status)
    
    ind <- which(island_spec[,1]==extinct)
    typeofspecies = island_spec[ind,4]
    
    if(typeofspecies == "I")
    {
      island_spec = island_spec[-ind, ] #remove immigrant
    }
    
    if(typeofspecies == "A")
    {
      island_spec = island_spec[-ind, ] #remove anagenetic
    }
    
    if(typeofspecies == "C")
    {
      #remove cladogenetic
      #first find species with same ancestor AND arrival totaltime
      sisters = intersect(which(island_spec[,2] == island_spec[ind,2]),
                          which(island_spec[,3] == island_spec[ind,3]))
      survivors = sisters[which(sisters != ind)]
      
      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic
        island_spec[survivors,4] = "A"
        island_spec[survivors,c(5,6)] = c(NA,NA)
        island_spec[survivors,7] = "Clado_extinct"
        island_spec = island_spec[-ind, ]
      }
      
      if(length(sisters) >= 3)
      {
        numberofsplits = nchar(island_spec[ind,5])
        mostrecentspl = substring(island_spec[ind,5],numberofsplits)
        
        if(mostrecentspl=="B")
        {
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }
        
        motiftofind = paste(substring(island_spec[ind,5],1,numberofsplits-1),
                            sistermostrecentspl,sep = "")
        possiblesister = survivors[which(substring(island_spec[survivors,5],1,
                                                   numberofsplits) == motiftofind)]
        
        if(mostrecentspl == "A")
        {
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec[possiblesister,6] == 
                                            min(as.numeric(island_spec[possiblesister,6])))]
          island_spec[tochange,6] = island_spec[ind, 6]
        }
        #remove the offending A/B from these species
        island_spec[possiblesister,5] = paste(substring(island_spec[possiblesister,5],1,numberofsplits - 1),
                                                     substring(island_spec[possiblesister,5],numberofsplits + 1,
                                                               nchar(island_spec[possiblesister,5])),sep = "")
        island_spec = island_spec[-int, ]
      }
    }
    island_spec = rbind(island_spec)
  }
  
  ##########################################
  # [7]: cladogenesis event with animal species
  if(possible_event$L1 == 7){
    
    tosplitID <- possible_event$Var1

    Mt <- new_Mt_clado_a(Mt=Mt,possible_event=possible_event,mutualism_pars=mutualism_pars)
    a_status[tosplitID] <- 0
    a_status <- c(a_status,1,1)
    tosplit <- tosplitID + length(p_status)
    ind <- which(island_spec[,1]==tosplit)
    #if the species that speciates is cladogenetic
    if(island_spec[ind,4] == "C")
    {
      #for daughter A
      island_spec[ind, 1] <- maxanimalID + 1
      oldstatus <- island_spec[ind, 5]
      island_spec[ind, 5] = paste(oldstatus,"A",sep = "")
      island_spec[ind, 6] = timeval
      island_spec[ind, 7] = NA
      island_spec[ind,8] = "animal"
      #for daughter B
      island_spec <- rbind(island_spec,c(maxanimalID + 2,island_spec[ind, 2],island_spec[ind, 3],
                                                       "C",paste(oldstatus,"B",sep = ""),timeval,NA,"animal"))
      maxanimalID = maxanimalID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      island_spec[ind, 4] = "C"
      island_spec[ind, 1] = maxanimalID + 1
      island_spec[ind, 5] = "A"
      island_spec[ind, 6] = island_spec[ind, 3]
      island_spec[ind, 7] = NA
      island_spec[ind,8] = "animal"
      #for daughter B
        island_spec = rbind(island_spec,c(maxanimalID + 2,island_spec[ind, 2],
                                                      island_spec[ind,3],"C","B",timeval,NA,"animal"))
        maxanimalID = maxanimalID + 2
    }
  }
  
  ##########################################
  # [8]: anagenesis event with animal species
  if(possible_event$L1 == 8){
    
    anagenesisID <- possible_event$Var1
    
    Mt <- new_Mt_ana_a(Mt=Mt,possible_event=possible_event,mutualism_pars=mutualism_pars)
    a_status[anagenesisID] <- 0
    a_status <- c(a_statusID,1)
    anagenesis <- anagenesisID + length(p_status)
    
    ind <- which(island_spec[,1]==anagenesis)
    maxanimalID <- maxanimalID + 1
    island_spec[ind, 4] = "A"
    island_spec[ind, 1] = maxanimalID
    island_spec[ind, 7] = "Immig_parent"
    island_spec[ind,8] = "animal"
    
  }
  
  ##########################################
  # [9]: cospeciation event between pairs
  if(possible_event$L1 == 9){
    
    cospec_plant <- possible_event$Var1
    cospec_animalID <- possible_event$Var2
    cospec_animal <- cospec_animalID + length(p_status)
    ind1 <- which(island_spec[, 1] == cospec_plant)
    ind2 <- which(island_spec[, 1] == cospec_animal)
    
    Mt <- new_Mt_cospec(Mt=Mt,possible_event=possible_event,mutualism_pars=mutualism_pars)
    p_status[cospec_plant] <- 0
    a_status[cospec_animalID] <- 0
    p_status <- c(p_status,1,1)
    a_status <- c(a_status,1,1)
    
    #for plant species, if the species that speciates is cladogenetic
    if(island_spec[ind1,4] == "C")
    {
      #for daughter A
      island_spec[ind1, 1] <- maxplantID + 1
      oldstatus <- island_spec[ind1, 5]
      island_spec[ind1, 5] = paste(oldstatus,"A",sep = "")
      island_spec[ind1, 6] = timeval
      island_spec[ind1, 7] = NA
      island_spec[ind1, 8] = "plant"
      #for daughter B
      island_spec <- rbind(island_spec,c(maxplantID + 2,island_spec[ind1, 2],island_spec[ind1, 3],
                                                     "C",paste(oldstatus,"B",sep = ""),timeval,NA,"plant"))
      maxplantID = maxplantID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      island_spec[ind1, 4] = "C"
      island_spec[ind1, 1] = maxplantID + 1
      island_spec[ind1, 5] = "A"
      island_spec[ind1, 6] = island_spec[ind1, 3]
      island_spec[ind1, 7] = NA
      island_spec[ind1, 8] = "plant"
      
      #for daughter B
      island_spec = rbind(island_spec,c(maxplantID + 2,island_spec[ind1, 2],
                                                    island_spec[ind1,3],"C","B",timeval,NA,"plant"))
      maxplantID = maxplantID + 2
    }
    #for animal species, if the species that speciates is cladogenetic
    
    if(island_spec[ind2,4] == "C")
    {
      #for daughter A
      island_spec[ind2, 1] <- maxanimalID + 1
      oldstatus <- island_spec[ind2, 5]
      island_spec[ind2, 5] = paste(oldstatus,"A",sep = "")
      island_spec[ind2, 6] = timeval
      island_spec[ind2, 7] = NA
      island_spec[ind1, 8] = "animal"
      #for daughter B
      island_spec <- rbind(island_spec,c(maxanimalID + 2,island_spec[ind2, 2],island_spec[ind2, 3],
                                                       "C",paste(oldstatus,"B",sep = ""),timeval,NA,"animal"))
      maxanimalID = maxanimalID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      island_spec[ind2, 4] = "C"
      island_spec[ind2, 1] = maxanimalID + 1
      island_spec[ind2, 5] = "A"
      island_spec[ind2, 6] = island_spec[ind2, 3]
      island_spec[ind2, 7] = NA
      island_spec[ind1, 8] = "animal"
      #for daughter B
      island_spec = rbind(island_spec,c(maxanimalID + 2,island_spec[ind2, 2],
                                                      island_spec[ind2,3],"C","B",timeval,NA,"animal"))
      maxanimalID = maxanimalID + 2
    }
  }
  
  ##########################################
  # [10]: gain links event between pairs
  if(possible_event$L1 == 10){
    
    togain_plant <- possible_event$Var1
    togain_animal <- possible_event$Var2
    Mt[togain_plant,togain_animal] == 1
  } 
  
  ##########################################
  # [11]: lose links event between pairs
  if(possible_event$L1 == 11){# [11]: loss links event between pairs
    
    tolose_plant <- possible_event$Var1
    tolose_animal <- possible_event$Var2
    Mt[tolose_plant,tolose_animal] == 0
  }
  
  ##########################################
  stt_table <- rbind(stt_table,
                     c(totaltime - timeval,
                       length(intersect(which(island_spec[,4] == "I"),which(island_spec[,8] == "plant"))),   
                       length(intersect(which(island_spec[,4] == "A"),which(island_spec[,8] == "plant"))),    
                       length(intersect(which(island_spec[,4] == "C"),which(island_spec[,8] == "plant"))),    
                       length(intersect(which(island_spec[,4] == "I"),which(island_spec[,8] == "animal"))),   
                       length(intersect(which(island_spec[,4] == "A"),which(island_spec[,8] == "animal"))),   
                       length(intersect(which(island_spec[,4] == "C"),which(island_spec[,8] == "animal")))))   
  
  updated_state <- list(Mt = Mt,
                        p_status = p_status,
                        a_status = a_status,
                        island_spec = island_spec,
                        maxplantID = maxplantID,
                        maxanimalID = maxanimalID,
                        stt_table = stt_table)
  return(updated_state)
}

