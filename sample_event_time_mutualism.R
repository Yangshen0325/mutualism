# Samples what event to happen next

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

#example:
# library(reshape2)
# load("rates.RData")
# update_rates_mutual(gam0_lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain,
#                   M0,Mt,p_status,a_status)
# event <- possible_event(rates) 
possible_event <- function(rates){
  
  output <- melt(setNames(rates,seq_along(rates)))
  
  x <- sample(1:dim(output)[1],
         size = 1,
         replace = FALSE,
         prob=unlist(rates))
  
  event <- output[x,]

  return(event)
}

# Calculates when the next timestep will be.

# rates <- update_rates_mutual(gam0_par,mu_par,laa_par,lac0_par,lambda1,K_par,w,
#                   M0,Mt,p_status,a_status,qloss,qgain,num_immigrants)
# timeval <- 1
# calc_next_timeval(rates,timeval)

# update_rates_mutual(gam0_lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain,
#                   M0,Mt,p_status,a_status)

calc_next_timeval <- function(rates, timeval) {
  
  output <- melt(setNames(rates,seq_along(rates)))
  
  totalrate <- sum(output$value)
  
  dt <- stats::rexp(1, totalrate)
  timeval <- timeval + dt
  
  return(list(timeval = timeval, dt = dt))
}


