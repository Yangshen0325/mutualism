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
new_Mt_ana_a <- function(Mt,event,p){  #??? what if c(0,0,0,0)
  
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

#
update_state <- function(
  p,
  Mt,
  p_status,
  a_status,
  rates){
  
  event <- possible_event(rates=rates)
  
  if(event$L1 == 1){# [1]: immigration event with plant species
    
    Mt <- Mt
    h <- event$Var1
    p_status[h] <- 1
    a_status <- a_status
    
  } else
    if(event$L1 == 2){# [2]: extinction event with plant species
      
      Mt <- Mt
      h <- event$Var1
      p_status[h] <- 0
      a_status <- a_status
      
    } else
      if(event$L1 == 3){# [3]: cladogenesis event with plant species
        
        Mt <- new_Mt_clado_p(Mt=Mt,event=event,p=p)
        h <- event$Var1
        p_status[h] <- 0
        p_status <- c(p_status,1,1)
        a_status <- a_status
        
      } else
        if(event$L1 == 4){# [4]: anagenesis event with plant species
          
          Mt <- new_Mt_ana_p(Mt=Mt,event=event,p=p)
          h <- event$Var1
          p_status[h] <- 0
          p_status <- c(p_status,1)
          a_status <- a_status
          
        } else
          if(event$L1 == 5){# [5]: immigration event with animal species
            
            Mt <- Mt
            h <- event$Var1
            p_status <- p_status
            a_status[h] <- 1
            
          } else
            if(event$L1 == 6){# [6]: extinction event with animal species
              
              Mt <- Mt
              h <- event$Var1
              p_status <- p_status
              a_status[h] <- 0
              
            } else 
              if(event$L1 == 7){# [7]: cladogenesis event with animal species
                
                Mt <- new_Mt_clado_a(Mt=Mt,event=event,p=p)
                p_status <- p_status
                h <- event$Var1
                a_status[h] <- 0
                a_status <- c(a_status,1,1)
                
              } else
                if(event$L1 == 8){# [8]: anagenesis event with animal species
                  
                  Mt <- new_Mt_ana_a(Mt=Mt,event=event,p=p)
                  p_status <- p_status
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
                      p_status <- p_status
                      a_status <- a_status
                      
                    } else 
                      if(event$L1 == 11){# [11]: loss links event between pairs
                        
                        h <- event$Var1
                        k <- event$Var2
                        Mt[h,k] == 1
                        p_status <- p_status
                        a_status <- a_status
                        
                      }
  new_state <- list(Mt = Mt,
                    p_status = p_status,
                    a_status = a_status)
  return(new_state)
}

