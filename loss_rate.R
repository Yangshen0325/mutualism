# Calculate the gain rate.
#
# Mt: interaction matrix on island at time t, t=0, Mt is identical to the initial matrix in mainland, M0
# p_status: to show whether the plant species present on island, p_status at time 0 is like rep(0,NROW(M0))
# a_status: to show whether the animal species present on island,a_status at time 0 is like rep(0,NCOL(M0))
# qgain: a coefficient to mediate the gain process
#
# example
#
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status <- c(1,0,0,1)
# a_status <- c(1,1,1,0,0)
# qloss <- 0.1
# get_loss_rate(Mt,p_status,a_status,qloss)

get_loss_rate <- function(Mt,
                        p_status,
                        a_status,
                        qloss){
  
  expand_matrix_list <- get_expand_matrix(Mt = Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  
  loss_rate <- qloss * Mt * 
    (expand_matrix_list[[1]]*expand_matrix_list[[2]] + 
       expand_matrix_list[[2]]*(1-expand_matrix_list[[1]]) + 
       expand_matrix_list[[1]]*(1-expand_matrix_list[[2]]))
  
  return(loss_rate)
}

####seprate function