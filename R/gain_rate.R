# Calculate the gain rate.

# mutualism_pars <- list(lac_par,mu_par,K_par,gam_par,laa_par,M0,qgain,qloss,lambda1) 
# qgain: a probability that a pair of plant and animal species could get a link 
#
# example
#
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
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