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
DAISIE_sim_update_state_mutualism <- function(
  timeval,
  totaltime,
  possible_event,
  Mt,
  p_status,
  a_status,
  maxplantID,
  maxanimalID,
  island_spec_plant,
  island_spec_animal,
  stt_table,
  mutualism_pars){
  
  ##########################################
  # [1]: immigration event with plant species 
  if(possible_event$L1 == 1){
    
    colonist = possible_event$Var1
    p_status[colonist] <- 1
    
    if (length(island_spec_plant[,1]) != 0)
    {
      isitthere <- which(island_spec_plant[,1] == colonist)
    }else{
      isitthere <- c()
    }
    if(length(isitthere) == 0)
    {
      island_spec_plant <- rbind(island_spec_plant,c(colonist,colonist,timeval,"I",NA,NA,NA))
    }
    if (length(isitthere) != 0) {
      island_spec_plant[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)
    }
  } 
  
  ##########################################
  # [2]: extinction event with plant species
  if(possible_event$L1 == 2){
    
    extinct =  possible_event$Var1
    p_status[extinct] <- 0
    
    ind <- which(island_spec_plant[,1]==extinct)
    typeofspecies = island_spec_plant[ind,4] 
    
    if(typeofspecies == "I")
    {
      island_spec_plant = island_spec_plant[-ind, ] #remove immigrant
    }
    
    if(typeofspecies == "A")
    {
      island_spec_plant = island_spec_plant[-ind, ] #remove anagenetic
    }
    
    if(typeofspecies == "C")##############################
    {
      #remove cladogenetic
      #first find species with same ancestor AND arrival totaltime
      sisters = intersect(which(island_spec_plant[,2] == island_spec_plant[ind,2]),
                          which(island_spec_plant[,3] == island_spec_plant[ind,3]))
      survivors = sisters[which(sisters != ind)]
      
      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic
        island_spec_plant[survivors,4] = "A"
        island_spec_plant[survivors,c(5,6)] = c(NA,NA)
        island_spec_plant[survivors,7] = "Clado_extinct"
        island_spec_plant = island_spec_plant[-ind, ]
      }
      
      if(length(sisters) >= 3)
      {
        numberofsplits = nchar(island_spec_plant[ind,5])
        mostrecentspl = substring(island_spec_plant[ind,5],numberofsplits)
        
        if(mostrecentspl=="B")
        {
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }
        
        motiftofind = paste(substring(island_spec_plant[ind,5],1,numberofsplits-1),
                            sistermostrecentspl,sep = "")
        possiblesister = survivors[which(substring(island_spec_plant[survivors,5],1,
                                                   numberofsplits) == motiftofind)]
        
        if(mostrecentspl == "A")
        {
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec_plant[possiblesister,6] == 
                                            min(as.numeric(island_spec_plant[possiblesister,6])))]
          island_spec_plant[tochange,6] = island_spec_plant[ind, 6]
        }
        #remove the offending A/B from these species
        island_spec_plant[possiblesister,5] = paste(substring(island_spec_plant[possiblesister,5],1,numberofsplits - 1),
                                                    substring(island_spec_plant[possiblesister,5],numberofsplits + 1,
                                                              nchar(island_spec_plant[possiblesister,5])),sep = "")
        island_spec_plant = island_spec_plant[-int, ]
      }
    }
    island_spec_plant = rbind(island_spec_plant)
  }
  
  ##########################################
  # [3]: cladogenesis event with plant species
  if(possible_event$L1 == 3){
    
    tosplit <- possible_event$Var1
    ind <- which(island_spec_plant[,1]==tosplit)
    
    Mt <- new_Mt_clado_p(Mt=Mt, possible_event=possible_event, p=p)
    p_status[tosplit] <- 0
    p_status <- c(p_status,1,1)
    
    #if the species that speciates is cladogenetic
    if(island_spec_plant[ind,4] == "C")
    {
      #for daughter A
      island_spec_plant[ind, 1] <- maxplantID + 1
      oldstatus <- island_spec_plant[ind, 5]
      island_spec_plant[ind, 5] = paste(oldstatus,"A",sep = "")
      island_spec_plant[ind, 6] = timeval
      island_spec_plant[ind, 7] = NA
      
      #for daughter B
      island_spec_plant <- rbind(island_spec_plant,c(maxplantID + 2,island_spec_plant[ind, 2],island_spec_plant[ind, 3],
                                                     "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      maxplantID = maxplantID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      island_spec_plant[ind, 4] = "C"
      island_spec_plant[ind, 1] = maxplantID + 1
      island_spec_plant[ind, 5] = "A"
      island_spec_plant[ind, 6] = island_spec_plant[ind, 3]
      island_spec_plant[ind, 7] = NA
      
      #for daughter B
      island_spec_plant = rbind(island_spec_plant,c(maxplantID + 2,island_spec_plant[ind, 2],
                                                    island_spec_plant[ind,3],"C","B",timeval,NA))
      maxplantID = maxplantID + 2
    }
  } 
  
  ##########################################
  # [4]: anagenesis event with plant species 
  if(possible_event$L1 == 4){
    
    anagenesis = possible_event$Var1
    
    Mt <- new_Mt_ana_p(Mt=Mt, possible_event=possible_event, p=p)
    p_status[anagenesis] <- 0
    p_status <- c(p_status,1)
    
    ind <- which(island_spec_plant[,1]==anagenesis)
    island_spec_plant[ind, 4] = "A"
    island_spec_plant[ind, 1] = maxplantID + 1
    island_spec_plant[ind, 7] = "Immig_parent"
    
    maxplantID <- maxplantID + 1
  } 
  
  ##########################################
  # [5]: immigration event with animal species
  if(possible_event$L1 == 5){
    
    colonist <- possible_event$Var1
    a_status[colonist] <- 1
    
    if (length(island_spec_animal[,1]) != 0)
    {
      isitthere <- which(island_spec_animal[,1] == colonist)
    }else{
      isitthere <- c()
    }
    if(length(isitthere) == 0) {
      island_spec_animal <- rbind(island_spec_animal,c(colonist,colonist,timeval,"I",NA,NA,NA))
    }
    if(length(isitthere) != 0) {
      island_spec_animal[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA)
    }
  }
  
  ##########################################
  # [6]: extinction event with animal species
  if(possible_event$L1 == 6){
    
    extinct <- possible_event$Var1
    a_status[extinct] <- 0
    
    ind <- which(island_spec_animal[,1]==extinct)
    typeofspecies = island_spec_animal[ind,4]
    
    if(typeofspecies == "I")
    {
      island_spec_animal = island_spec_animal[-ind, ] #remove immigrant
    }
    
    if(typeofspecies == "A")
    {
      island_spec_animal = island_spec_animal[-ind, ] #remove anagenetic
    }
    
    if(typeofspecies == "C")
    {
      #remove cladogenetic
      #first find species with same ancestor AND arrival totaltime
      sisters = intersect(which(island_spec_animal[,2] == island_spec_animal[ind,2]),
                          which(island_spec_animal[,3] == island_spec_animal[ind,3]))
      survivors = sisters[which(sisters != ind)]
      
      if(length(sisters) == 2)
      {
        #survivors status becomes anagenetic
        island_spec_animal[survivors,4] = "A"
        island_spec_animal[survivors,c(5,6)] = c(NA,NA)
        island_spec_animal[survivors,7] = "Clado_extinct"
        island_spec_animal = island_spec_animal[-ind, ]
      }
      
      if(length(sisters) >= 3)
      {
        numberofsplits = nchar(island_spec_animal[ind,5])
        mostrecentspl = substring(island_spec_animal[ind,5],numberofsplits)
        
        if(mostrecentspl=="B")
        {
          sistermostrecentspl = "A"
        }
        if(mostrecentspl=="A")
        {
          sistermostrecentspl = "B"
        }
        
        motiftofind = paste(substring(island_spec_animal[ind,5],1,numberofsplits-1),
                            sistermostrecentspl,sep = "")
        possiblesister = survivors[which(substring(island_spec_animal[survivors,5],1,
                                                   numberofsplits) == motiftofind)]
        
        if(mostrecentspl == "A")
        {
          #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
          tochange = possiblesister[which(island_spec_animal[possiblesister,6] == 
                                            min(as.numeric(island_spec_animal[possiblesister,6])))]
          island_spec_animal[tochange,6] = island_spec_animal[ind, 6]
        }
        #remove the offending A/B from these species
        island_spec_animal[possiblesister,5] = paste(substring(island_spec_animal[possiblesister,5],1,numberofsplits - 1),
                                                     substring(island_spec_animal[possiblesister,5],numberofsplits + 1,
                                                               nchar(island_spec_animal[possiblesister,5])),sep = "")
        island_spec_animal = island_spec_animal[-int, ]
      }
    }
    island_spec_animal = rbind(island_spec_animal)
  }
  
  ##########################################
  # [7]: cladogenesis event with animal species
  if(possible_event$L1 == 7){
    
    tosplit <- possible_event$Var1
    ind <- which(island_spec_animal[,1]==tosplit)
    
    Mt <- new_Mt_clado_a(Mt=Mt,possible_event=possible_event,p=p)
    a_status[tosplit] <- 0
    a_status <- c(a_status,1,1)
    #if the species that speciates is cladogenetic
    if(island_spec_animal[ind,4] == "C")
    {
      #for daughter A
      island_spec_animal[ind, 1] <- maxanimalID + 1
      oldstatus <- island_spec_animal[ind, 5]
      island_spec_animal[ind, 5] = paste(oldstatus,"A",sep = "")
      island_spec_animal[ind, 6] = timeval
      island_spec_animal[ind, 7] = NA
      
      #for daughter B
      island_spec_animal <- rbind(island_spec_animal,c(maxanimalID + 2,island_spec_animal[ind, 2],island_spec_animal[ind, 3],
                                                       "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      maxanimalID = maxanimalID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      island_spec_animal[ind, 4] = "C"
      island_spec_animal[ind, 1] = maxanimalID + 1
      island_spec_animal[ind, 5] = "A"
      island_spec_animal[ind, 6] = island_spec_animal[ind, 3]
      island_spec_animal[ind, 7] = NA
      
      #for daughter B
        island_spec_animal = rbind(island_spec_animal,c(maxanimalID + 2,island_spec_animal[ind, 2],
                                                      island_spec_animal[ind,3],"C","B",timeval,NA))
      maxanimalID = maxanimalID + 2
    }
  }
  
  ##########################################
  # [8]: anagenesis event with animal species
  if(possible_event$L1 == 8){
    
    anagenesis <- possible_event$Var1
    
    Mt <- new_Mt_ana_a(Mt=Mt,possible_event=possible_event,p=p)
    a_status[anagenesis] <- 0
    a_status <- c(a_status,1)
    
    ind <- which(island_spec_animal[,1]==anagenesis)
    island_spec_animal[ind, 4] = "A"
    island_spec_animal[ind, 1] = maxanimalID + 1
    island_spec_animal[ind, 7] = "Immig_parent"
    
    maxanimalID <- maxanimalID + 1
  }
  
  ##########################################
  # [9]: cospeciation event between pairs
  if(possible_event$L1 == 9){
    
    cospec_plant <- possible_event$Var1
    cospec_animal <- possible_event$Var2
    ind1 <- which(island_spec_plant[, 1] == cospec_plant)
    ind2 <- which(island_spec_animal[, 1] == cospec_animal)
    
    Mt <- new_Mt_cospec(Mt=Mt,possible_event=possible_event,p=p)
    p_status[cospec_plant] <- 0
    a_status[cospec_animal] <- 0
    p_status <- c(p_status,1,1)
    a_status <- c(a_status,1,1)
    
    #for plant species, is the species that speciates is cladogenetic
    if(island_spec_plant[ind1,4] == "C")
    {
      #for daughter A
      island_spec_plant[ind1, 1] <- maxplantID + 1
      oldstatus <- island_spec_plant[ind1, 5]
      island_spec_plant[ind1, 5] = paste(oldstatus,"A",sep = "")
      island_spec_plant[ind1, 6] = timeval
      island_spec_plant[ind1, 7] = NA
      
      #for daughter B
      island_spec_plant <- rbind(island_spec_plant,c(maxplantID + 2,island_spec_plant[ind1, 2],island_spec_plant[ind1, 3],
                                                     "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      maxplantID = maxplantID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      island_spec_plant[ind1, 4] = "C"
      island_spec_plant[ind1, 1] = maxplantID + 1
      island_spec_plant[ind1, 5] = "A"
      island_spec_plant[ind1, 6] = island_spec_plant[ind1, 3]
      island_spec_plant[ind1, 7] = NA
      
      #for daughter B
      island_spec_plant = rbind(island_spec_plant,c(maxplantID + 2,island_spec_plant[ind1, 2],
                                                    island_spec_plant[ind1,3],"C","B",timeval,NA))
      maxplantID = maxplantID + 2
    }
    #for animal species, if the species that speciates is cladogenetic
    
    if(island_spec_animal[ind2,4] == "C")
    {
      #for daughter A
      island_spec_animal[ind2, 1] <- maxanimalID + 1
      oldstatus <- island_spec_animal[ind2, 5]
      island_spec_animal[ind2, 5] = paste(oldstatus,"A",sep = "")
      island_spec_animal[ind2, 6] = timeval
      island_spec_animal[ind2, 7] = NA
      
      #for daughter B
      island_spec_animal <- rbind(island_spec_animal,c(maxanimalID + 2,island_spec_animal[ind2, 2],island_spec_animal[ind2, 3],
                                                       "C",paste(oldstatus,"B",sep = ""),timeval,NA))
      maxanimalID = maxanimalID + 2
    } else {
      #if the species that speciates is not cladogenetic
      
      #for daughter A
      island_spec_animal[ind2, 4] = "C"
      island_spec_animal[ind2, 1] = maxanimalID + 1
      island_spec_animal[ind2, 5] = "A"
      island_spec_animal[ind2, 6] = island_spec_animal[ind2, 3]
      island_spec_animal[ind2, 7] = NA
      
      #for daughter B
      island_spec_animal = rbind(island_spec_animal,c(maxanimalID + 2,island_spec_animal[ind2, 2],
                                                      island_spec_animal[ind2,3],"C","B",timeval,NA))
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
                       length(which(island_spec_plant[, 4] == "I")),
                       length(which(island_spec_plant[, 4] == "A")),
                       length(which(island_spec_plant[, 4] == "C")),
                       length(which(island_spec_animal[, 4] == "I")),
                       length(which(island_spec_animal[, 4] == "A")),
                       length(which(island_spec_animal[, 4] == "C"))))
  
  updated_state <- list(Mt = Mt,
                        p_status = p_status,
                        a_status = a_status,
                        island_spec_plant = island_spec_plant,
                        island_spec_animal = island_spec_animal,
                        maxplantID = maxplantID,
                        maxanimalID = maxanimalID,
                        stt_table = stt_table)
  return(updated_state)
}

