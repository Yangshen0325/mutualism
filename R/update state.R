# update Mt p_status a_status

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

# Mt = {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# M0 = {set.seed(2);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(1,0,1,1)
# a_status<-c(1,1,1,1,1)
# event <- possible_event(rates) #   Var1 Var2 value L1
#    1    1   0.7  8
# new Mt if cladogenesis happends with plant species
new_Mt_clado_p <- function(Mt,event,p){
  
  possible_output <- list(c(1,1), c(1,0), c(0,1))
  newrows <- list()
  h <- event$Var1
  newrows[c(which(Mt[h,]==0))] <- list(c(0,0))
  newrows[c(which(Mt[h,]==1))] <- sample(possible_output, 
                                         size=length(which(Mt[h,]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)/2,(1-p)/2))
  newrows<-matrix(unlist(newrows),nrow=2,ncol=NCOL(Mt))
  
  Mt <- rbind(Mt,newrows)
  return(Mt)
}

# new Mt if anagenesis happends with plant species
new_Mt_ana_p <- function(Mt,event,p){
  
  newrows <- c()
  h <- event$Var1
  newrows[c(which(Mt[h,]==0))] <- 0
  newrows[c(which(Mt[h,]==1))] <- sample(c(1,0), 
                                         size=length(which(Mt[h,]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)))
  Mt <- rbind(Mt,newrows)
  return(Mt)
}
# if cladogenesis happends with animal species
new_Mt_clado_a <- function(Mt,event,p){
  
  possible_output <- list(c(1,1), c(1,0), c(0,1))
  newcols <- list()
  h <- event$Var1
  newcols[c(which(Mt[,h]==0))] <- list(c(0,0))
  newcols[c(which(Mt[,h]==1))] <- sample(possible_output, 
                                         size=length(which(Mt[,h]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)/2,(1-p)/2))
  
  newcols<-t(matrix(unlist(newcols),nrow=2,ncol=NROW(Mt)))
  Mt <- cbind(Mt,newcols)
  return(Mt)
}

# if anagenesis happends with animal species
# Mt = {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
new_Mt_ana_a <- function(Mt,event,p){  
  
  newcols <- c()
  h <- event$Var1
  newcols[c(which(Mt[,h]==0))] <- 0
  newcols[c(which(Mt[,h]==1))] <- sample(c(1,0), 
                                         size=length(which(Mt[,h]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)))
  
  Mt <- cbind(Mt,newcols)
  return(Mt)
}

# if cospeciation happends with plant species i and animal species j
# Mt = {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
new_Mt_cospec <- function(Mt,event,p){
  
  possible_output <- list(c(1,1), c(1,0), c(0,1))
  newrows <- list()
  newcols <- list()
  h <- event$Var1
  k <- event$Var2
  
  newcols[c(which(Mt[,k]==0))] <- list(c(0,0))
  newcols[c(which(Mt[,k]==1))] <- sample(possible_output, 
                                         size=length(which(Mt[,k]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)/2,(1-p)/2))
  
  newcols<-t(matrix(unlist(newcols),nrow=2,ncol=NROW(Mt)))
  
  newrows[c(which(Mt[h,]==0))] <- list(c(0,0))
  newrows[c(which(Mt[h,]==1))] <- sample(possible_output, 
                                         size=length(which(Mt[h,]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)/2,(1-p)/2))
  
  newrows <- matrix(unlist(newrows),nrow=2,ncol=NCOL(Mt))
  newrows <- cbind(newrows, diag(1,2,2))
  Mt <- cbind(Mt,newcols)
  Mt <- rbind(Mt,newrows)
  
  return(Mt)
}

#Updates state of island given sampled event
update_state <- function(
  timeval,
  totaltime,
  possible_event,
  p,
  Mt,
  p_status,
  a_status){
  
  if(possible_event$L1 == 1){# [1]: immigration event with plant species
    
    colonist = possible_event$Var1
    p_status[colonist] <- 1
    
    if (length(island_spec[,1]) != 0)
    {
      isitthere = which(island_spec[,1] == colonist)
      island_spec[isitthere,] = c(colonist,colonist,timeval,"I",NA,NA,NA,"p")
    } else
    {
      island_spec = rbind(island_spec,c(colonist,colonist,timeval,"I",NA,NA,NA,"p"))
    }
    
  } else
    if(possible_event$L1 == 2){# [2]: extinction event with plant species
      
      island_spec_plant = which(island_spec[,8] == "p")
      extinctp = DDD::sample2(island_spec_plant,1)
      p_status[extinctp] <- 0
      
      typeofspecies = island_spec[extinctp,4]
      
      if(typeofspecies == "I")
      {
        island_spec = island_spec[-extinctp,] #remove immigrant
      }
      
      if(typeofspecies == "A")
      {
        island_spec = island_spec[-extinctp,] #remove anagenetic
      }
      
      if(typeofspecies == "C")
      {
        #remove cladogenetic
        #first find species with same ancestor AND arrival totaltime
        sisters = intersect(which(island_spec[,2] == island_spec[extinctp,2]),which(island_spec[,3] == island_spec[extinctp,3]))
        survivors = sisters[which(sisters != extinctp)]
        
        if(length(sisters) == 2)
        {
          #survivors status becomes anagenetic
          island_spec[survivors,4] = "A"
          island_spec[survivors,c(5,6)] = c(NA,NA)
          island_spec[survivors,7] = "Clado_extinct"
          island_spec = island_spec[-extinctp,]
        }
        
        if(length(sisters) >= 3)
        {
          numberofsplits = nchar(island_spec[extinctp,5])
          mostrecentspl = substring(island_spec[extinct,5],numberofsplits)
          
          if(mostrecentspl=="B")
          {
            sistermostrecentspl = "A"
          }
          if(mostrecentspl=="A")
          {
            sistermostrecentspl = "B"
          }
          
          motiftofind = paste(substring(island_spec[extinct,5],1,numberofsplits-1),sistermostrecentspl,sep = "")
          possiblesister = survivors[which(substring(island_spec[survivors,5],1,numberofsplits) == motiftofind)]
          
          if(mostrecentspl == "A")
          {
            #change the splitting date of the sister species so that it inherits the early splitting that used to belong to A.
            tochange = possiblesister[which(island_spec[possiblesister,6] == min(as.numeric(island_spec[possiblesister,6])))]
            island_spec[tochange,6] = island_spec[extinctp,6]
          }
          
          #remove the offending A/B from these species
          island_spec[possiblesister,5] = paste(substring(island_spec[possiblesister,5],1,numberofsplits - 1),
                                                substring(island_spec[possiblesister,5],numberofsplits + 1,
                                                          nchar(island_spec[possiblesister,5])),sep = "")
          island_spec = island_spec[-extinctp,]
        }
      }
      island_spec = rbind(island_spec)
    }else
      
      if(possible_event$L1 == 3){# [3]: cladogenesis event with plant species
        
        island_spec_plant = which(island_spec[,8] == "p")
        tosplitp = DDD::sample2(island_spec_plant,1)
        
        Mt <- new_Mt_clado_p(Mt=Mt,possible_event=possible_event,p=p)
        p_status[tosplitp] <- 0
        p_status <- c(p_status,1,1)
        
        if(island_spec[tosplitp,4] == "C")
        {
          #for daughter A
          
          island_spec[tosplitp,4] = "C"
          island_spec[tosplitp,1] = maxspecID + 1
          oldstatus = island_spec[tosplitp,5]
          island_spec[tosplitp,5] = paste(oldstatus,"A",sep = "")
          #island_spec[tosplit,6] = timeval
          island_spec[tosplitp,7] = NA
          island_spec[tosplitp,8] = "1"
          
          #for daughter B
          island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplitp,2],island_spec[tosplitp,3],
                                            "C",paste(oldstatus,"B",sep = ""),timeval,NA,1))
          maxspecID = maxspecID + 2
        } else {
          #if the species that speciates is not cladogenetic
          
          #for daughter A
          
          island_spec[tosplitp,4] = "C"
          island_spec[tosplitp,1] = maxspecID + 1
          island_spec[tosplitp,5] = "A"
          island_spec[tosplitp,6] = island_spec[tosplitp,3]
          island_spec[tosplitp,7] = NA
          island_spec[tosplitp,8] = "1"
          
          #for daughter B
          island_spec = rbind(island_spec,c(maxspecID + 2,island_spec[tosplitp,2],island_spec[tosplitp,3],"C","B",timeval,NA,1))
          maxspecID = maxspecID + 2
        }
        
      } else
        if(possible_event$L1 == 4){# [4]: anagenesis event with plant species
          
          immi_specs = intersect(which(island_spec[,4] == "I"), which(island_spec[,8] == "p"))
          anagenesisp = DDD::sample2(immi_specs,1)
          
          Mt <- new_Mt_ana_p(Mt=Mt,possible_event,p=p)
          p_status[anagenesisp] <- 0
          p_status <- c(p_status,1)
          
          maxspecID = maxspecID + 1
          island_spec[anagenesisp,4] = "A"
          island_spec[anagenesisp,1] = maxspecID
          island_spec[anagenesisp,7] = "Immig_parent"
          island_spec[anagenesisp,8] = "p"
          
        } else
          if(possible_event$L1 == 5){# [5]: immigration event with animal species
            
            Mt <- Mt # can be deleted
            h <- event$Var1
            p_status <- p_status # can be deleted
            a_status[h] <- 1
            
          } else
            if(event$L1 == 6){# [6]: extinction event with animal species
              
              Mt <- Mt # can be deleted
              h <- event$Var1
              p_status <- p_status # can be deleted
              a_status[h] <- 0
              
            } else 
              if(event$L1 == 7){# [7]: cladogenesis event with animal species
                
                Mt <- new_Mt_clado_a(Mt=Mt,event=event,p=p)
                p_status <- p_status # can be deleted
                h <- event$Var1
                a_status[h] <- 0
                a_status <- c(a_status,1,1)
                
              } else
                if(event$L1 == 8){# [8]: anagenesis event with animal species
                  
                  Mt <- new_Mt_ana_a(Mt=Mt,event=event,p=p)
                  p_status <- p_status # can be deleted
                  h <- event$Var1
                  a_status[h] <- 0
                  a_status <- c(a_status,1)
                  
                } else 
                  if(event$L1 == 9){# [9]: cospeciation event between pairs
                    
                    Mt <- new_Mt_cospec(Mt=Mt,event=event,p=p)
                    h <- event$Var1
                    k <- event$Var2
                    p_status[h] <- 0
                    a_status[k] <- 0
                    p_status <- c(p_status,1,1)
                    a_status <- c(a_status,1,1)
                    
                  } else 
                    if(event$L1 == 10){# [10]: gain links event between pairs
                      
                      h <- event$Var1
                      k <- event$Var2
                      Mt[h,k] == 1
                      p_status <- p_status # can be deleted
                      a_status <- a_status # can be deleted
                      
                    } else 
                      if(event$L1 == 11){# [11]: loss links event between pairs
                        
                        h <- event$Var1
                        k <- event$Var2
                        Mt[h,k] == 1
                        p_status <- p_status # can be deleted
                        a_status <- a_status # can be deleted
                        
                      }
  new_state <- list(Mt = Mt,
                    p_status = p_status,
                    a_status = a_status)
  return(new_state)
}

