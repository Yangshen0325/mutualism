# Calculate the gain rate.

# mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
# qgain: a probability that a pair of plant and animal species could get a link 
#

# qgain <- 0.1
# gain_rate <- get_gain_rate(Mt,p_status,a_status,mutualism_pars)

get_gain_rate <- function(Mt,
                          p_status,
                          a_status,
                          mutualism_pars){
  
  expand_matrix_list <- get_expand_matrix(Mt = Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  qgain <- mutualism_pars$qgain
  gain_rate <-  qgain * (1-Mt) * expand_matrix_list[[1]] * expand_matrix_list[[2]]
  
  return(gain_rate)
}

####seprate function