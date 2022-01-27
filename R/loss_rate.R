# Calculate the gain rate.

# mutualism_pars <- list(lac_par,mu_par,K_par,gam_par,laa_par,M0,qgain,qloss,lambda1) 
# qloss: a probability that a pair of plant and animal species could loss a link 
#
# example
#
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# qloss <- 0.1
# loss_rate <- get_loss_rate(Mt,p_status,a_status, mutualism_pars)

get_loss_rate <- function(Mt,
                        p_status,
                        a_status,
                        mutualism_pars){
  
  expand_matrix_list <- get_expand_matrix(Mt = Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  qloss <- mutualism_pars$qloss
  
  loss_rate <- qloss * Mt * 
    (expand_matrix_list[[1]]*expand_matrix_list[[2]] + 
       expand_matrix_list[[2]]*(1-expand_matrix_list[[1]]) + 
       expand_matrix_list[[1]]*(1-expand_matrix_list[[2]]))
  
  return(loss_rate)
}

####seprate function