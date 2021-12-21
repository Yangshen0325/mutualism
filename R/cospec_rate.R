
# Calculate cospeciation rate
#

# lambda1: a coefficient that mediate the effects from mutualism
# K_par includes: K_par <- c(KP0,KP1,KA0,KA1)
# KP0: the initial carrying capacity of plant species without any mutualists
# KA0: the initial carrying capacity of animal species without any mutualists
# KP1: a coefficient showing the influence from mutualism to plant species
# KA1: a coefficient showing the influence from mutualism to animal species
# Mt: the possible interaction matrix on island at time t
# p_status: to show whether the plant species present on island, p_status at time 0 is like rep(0,NROW(M0))
# a_status: to show whether the animal species present on island,a_status at time 0 is like rep(0,NCOL(M0))
#
# example:
#
# lambda1 <- 0.1
# K_par <- c(20,0.6,20,0.6)
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(1,0,0,1)
# a_status<-c(1,0,0,0,1)
#
# get_cospec_rate(lambda1,K_par,Mt,p_status,a_status)

get_cospec_rate <- function(lambda1,
                            K_par,
                            M0,
                            Mt,
                            p_status,
                            a_status){
  
  part_compe_list <- get_part_compe(M0=M0,
                                    Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
  
  NK_list <- get_NK(K_par=K_par,
                    M0=M0,
                    Mt=Mt,
                    p_status= p_status,
                    a_status= a_status)
  
  expand_matrix_list <- get_expand_matrix(Mt = Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  
  cospec_rate <- lambda1 * Mt * expand_matrix_list[[1]] * expand_matrix_list[[2]] *
    matrix(rep(NK_list[[3]],NCOL(Mt)), ncol = NCOL(Mt)) *
    t(matrix(rep(NK_list[[4]],NROW(Mt)),ncol = NROW(Mt)))
  
  return(cospec_rate)
}
