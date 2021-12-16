# some functions will be used
get_part_compe <- function(Mt,
                           p_status,
                           a_status){
  tMt <- t(Mt) 
  plant_part <- Mt %*% a_status
  animal_part <- tMt %*% p_status
  plant_compe <- matrix()
  animal_compe <- matrix()
  for (x in seq(ncol(tMt))){
    plant_compe[x] <- sum((colSums(tMt[,x]*tMt[,-x]* a_status)>=1) * as.matrix(p_status)[-x,])
  }
  for (x in seq(ncol(Mt))){
    animal_compe[x] <- sum((colSums(Mt[,x] * Mt[,-x]* p_status)>=1) * as.matrix(a_status)[-x,])
  }
  
  part_compe_list <- list(plant_part = plant_part,
                          animal_part = animal_part,
                          plant_compe = plant_compe,
                          animal_compe = animal_compe)
  
  return(part_compe_list)
}

get_NK <- function(K_par,
                   Mt,
                   p_status,
                   a_status){
  
  part_compe_list <- get_part_compe(Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
  
  indp <- which(part_compe_list[[1]]==0)
  inda <- which(part_compe_list[[2]]==0)
  
  plant_NK <- exp(-(sum(p_status)/K_par[1])+
                    part_compe_list[[3]]/(K_par[2]*part_compe_list[[1]]))
  plant_NK[indp] <- exp(-(sum(p_status)/K_par[1]))
  
  
  animal_NK <- exp(-(sum(a_status)/K_par[3])+
                     part_compe_list[[4]]/(K_par[4]*part_compe_list[[2]])) 
  animal_NK[inda] <- exp(-(sum(a_status)/K_par[3]))
  
  
  NK_list <- list(plant_NK = plant_NK,
                  animal_NK = animal_NK)
  
  return(NK_list)
}

get_expand_matrix <- function(Mt,
                              p_status,
                              a_status){
  expd_p_sta <- matrix(rep(p_status,NCOL(Mt)), ncol = NCOL(Mt))
  expd_a_sta <- t(matrix(rep(a_status,NROW(Mt)),ncol = NROW(Mt)))
  expand_matrix_list <- list(expd_p_sta = expd_p_sta,
                             expd_a_sta = expd_a_sta)
  return(expand_matrix_list)
}

#immigration rate and cladogenesis rate
get_immig_clado_rate <- function(gam0_lac0_par,
                                 K_par,
                                 Mt,
                                 p_status,
                                 a_status){
  
  NK_list <- get_NK(K_par=K_par,
                    Mt=Mt,
                    p_status= p_status,
                    a_status= a_status)
  
  plant_immig_rate <- gam0_lac0_par[1] * NK_list[[1]]
  animal_immig_rate <- gam0_lac0_par[2] * NK_list[[2]]
  
  plant_clado_rate <- gam0_lac0_par[3] * NK_list[[1]]
  animal_clado_rate <- gam0_lac0_par[4] * NK_list[[2]]
  
  immig_clado_list <- list(plant_immig_rate = plant_immig_rate,
                           animal_immig_rate = animal_immig_rate,
                           plant_clado_rate = plant_clado_rate,
                           animal_clado_rate = animal_clado_rate)
  return(immig_clado_list)
}

#extinction rate
get_ext_rate <- function(mu_par,
                         Mt,
                         p_status,
                         a_status) {
  
  part_compe_list <- get_part_compe(Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
  
  plant_ext_rate <- as.matrix(pmax(0, mu_par[1] - mu_par[2] * part_compe_list[[1]])) 
  animal_ext_rate <- as.matrix(pmax(0, mu_par[3] - mu_par[4] * part_compe_list[[2]])) 
  
  ext_list <- list(plant_ext_rate = plant_ext_rate,
                   animal_ext_rate = animal_ext_rate)
  return(ext_list)
}
#anagenetic rate
get_ana_rate <- function(laa_par,
                         M0,
                         Mt,
                         p_status,
                         a_status) {  
  
  plant_ana_rate =  laa_par[1] + 
    laa_par[2] *  abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0) %*% a_status[1:NCOL(M0)]
  animal_ana_rate =  laa_par[3] + 
    laa_par[4] * t(abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0)) %*% p_status[1:NROW(M0)]
  
  ana_list <- list(plant_ana_rate = plant_ana_rate,
                   animal_ana_rate = animal_ana_rate)
  return(ana_list)
}

#cospeciation rate
get_cospec_rate <- function(lambda1,
                            K_par,
                            Mt,
                            p_status,
                            a_status){
  
  part_compe_list <- get_part_compe(Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
  
  NK_list <- get_NK(K_par=K_par,
                    Mt=Mt,
                    p_status= p_status,
                    a_status= a_status)
  
  expand_matrix_list <- get_expand_matrix(Mt = Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  
  cospec_rate <- lambda1 * Mt * expand_matrix_list[[1]] * expand_matrix_list[[2]] *
    matrix(rep(NK_list[[1]],NCOL(Mt)), ncol = NCOL(Mt)) *
    t(matrix(rep(NK_list[[2]],NROW(Mt)),ncol = NROW(Mt)))
  
  return(cospec_rate)
}

#gain rate
get_gain_rate <- function(Mt,
                          p_status,
                          a_status,
                          qgain){
  
  expand_matrix_list <- get_expand_matrix(Mt = Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  
  gain_rate <-  qgain * (1-Mt) * expand_matrix_list[[1]] * expand_matrix_list[[2]]
  
  return(gain_rate)
}

#loss rate
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
#
#update all rates
#
#example

# gam0_lac0_par <- c(0.6,0.6,0.3,0.3)
# mu_par <- c(0.02,0.01,0.02,0.01)
# laa_par <- c(0.1,0.2,0.1,0.2)
# lac0_par <- c(0.3,0.3)
# lambda1 <- 0.1
# K_par <- c(20,0.6,20,0.6)
# M0 <- {set.seed(2);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(1,1,1,1)
# a_status<-c(1,1,1,1,1)
# qloss <- 0.1
# qgain <- 0.1
# #
# update_rates_mutual(gam0_lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain,
#                  M0,Mt,p_status,a_status)

update_rates_mutual <- function(gam0_lac0_par,
                                K_par,
                                mu_par,
                                laa_par,
                                lambda1,
                                qgain,
                                qloss,
                                M0,
                                Mt,
                                p_status,
                                a_status) {
  
  part_compe_list <- get_part_compe(Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
  
  NK_list <- get_NK(K_par=K_par,
                    Mt=Mt,
                    p_status= p_status,
                    a_status= a_status)
  
  expand_matrix_list <- get_expand_matrix(Mt = Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  
  #immigration rate and cladogenesis rate
  immig_clado_list <- get_immig_clado_rate(
    gam0_lac0_par = gam0_lac0_par,
    K_par = K_par,
    Mt = Mt,
    p_status = p_status,
    a_status = a_status)
  #testit::assert(is.list(immig_clado_list)) 
  
  #extinction rate
  ext_rate <- get_ext_rate(
    mu_par = mu_par,
    Mt = Mt,
    p_status = p_status,
    a_status = a_status)
  #testit::assert(is.list(ext_rate))
  
  #anagenetic rate
  ana_rate <- get_ana_rate(
    laa_par = laa_par,
    M0 = M0,
    Mt = Mt,
    p_status = p_status,
    a_status = a_status)
  #testit::assert(is.list(ana_rate))
  
  #cospeciation rate
  cospec_rate <- get_cospec_rate(
    lambda1 = lambda1,
    K_par = K_par,
    Mt = Mt,
    p_status = p_status,
    a_status = a_status)
  #testit::assert(is.matrix(cospec_rate))
  
  #gain rate
  gain_rate <- get_gain_rate(
    Mt = Mt,
    p_status = p_status,
    a_status = a_status,
    qgain = qgain)
  #testit::assert(is.matrix(gain_rate))
  
  #loss rate
  loss_rate <- get_loss_rate(
    Mt = Mt,
    p_status = p_status,
    a_status = a_status,
    qloss = qloss)
  #testit::assert(is.matrix(loss_rate))
  
  rates <- list(
    p_immig_rate = immig_clado_list$plant_immig_rate,
    p_ext_rate = ext_rate$plant_ext_rate,
    p_clado_rate = immig_clado_list$plant_clado_rate,
    p_ana_rate = ana_rate$plant_ana_rate,
    a_immig_rate = immig_clado_list$animal_immig_rate,
    a_ext_rate = ext_rate$animal_ext_rate,
    a_clado_rate = immig_clado_list$animal_clado_rate,
    a_ana_rate = ana_rate$animal_ana_rate,
    cospec_rate = cospec_rate,
    gain_rate = gain_rate,
    loss_rate = loss_rate)
  
  return(rates)
}

# save(rates, file="rates.Rdata")
