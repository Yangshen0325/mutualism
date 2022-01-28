
# Calculate cospeciation rate
# mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
# lambda1: a coefficient that mediate the effects from mutualism
# Mt: the possible interaction matrix on island at time t
# p_status: to show whether the plant species present on island, p_status at time 0 is like rep(0,NROW(M0))
# a_status: to show whether the animal species present on island,a_status at time 0 is like rep(0,NCOL(M0))
#
#
# cospec_rate <- get_cospec_rate(lambda1,K_par,Mt,p_status,a_status)

get_cospec_rate <- function(K,
                            Mt,
                            p_status,
                            a_status,
                            mutualism_pars){
  
  NK_list <- get_NK(K = K,
                    Mt = Mt,
                    p_status = p_status,
                    a_status = a_status,
                    mutualism_pars = mutualism_pars)
  
  expand_matrix_list <- get_expand_matrix(Mt = Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  lambda1 <- mutualism_pars$lambda1
  cospec_rate <- lambda1 * Mt * expand_matrix_list[[1]] * expand_matrix_list[[2]] *
    matrix(rep(NK_list[[3]],NCOL(Mt)), ncol = NCOL(Mt)) *
    t(matrix(rep(NK_list[[4]],NROW(Mt)),ncol = NROW(Mt)))
  
  return(cospec_rate)
}
