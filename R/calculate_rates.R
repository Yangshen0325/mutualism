
#immigration rate
get_immig_rate <- function(gam0_par,
                           K_par,
                           Mt,
                           M0,
                           p_status,
                           a_status){
  
  NK_list <- get_NK(K_par=K_par,
                    M0=M0,
                    Mt=Mt,
                    p_status= p_status,
                    a_status= a_status)
  
  plant_immig_rate <- gam0_par[1] * NK_list[[1]] 
  animal_immig_rate <- gam0_par[2] * NK_list[[2]]
  
  
  immig_list <- list(plant_immig_rate = plant_immig_rate,
                     animal_immig_rate = animal_immig_rate)
  
  return(immig_list)
}

#extinction rate
get_ext_rate <- function(mu_par,
                         Mt,
                         p_status,
                         a_status) {
  
  part_compe_list <- get_part_compe(M0=M0,
                                    Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
  
  plant_ext_rate <- as.matrix(pmax(0, mu_par[1] - mu_par[2] * part_compe_list[[3]])) * p_status
  animal_ext_rate <- as.matrix(pmax(0, mu_par[3] - mu_par[4] * part_compe_list[[4]])) *a_status
  
  ext_list <- list(plant_ext_rate = plant_ext_rate,
                   animal_ext_rate = animal_ext_rate)
  return(ext_list)
}
#anagenetic rate
get_ana_rate <- function(laa_par,
                         M0,
                         Mt,
                         p_status,
                         a_status,
                         island_spec_plant,
                         island_spec_animal) { 
  
  possible_ana_p <- matrix(0,nrow=NROW(M0))
  possible_ana_a <- matrix(0,nrow=NCOL(M0))
  
  plant_ind <-as.numeric(island_spec_plant[which(island_spec_plant[,4]=="I"),1])
  animal_ind <-as.numeric(island_spec_animal[which(island_spec_animal[,4]=="I"),1])
  possible_ana_p[plant_ind] = 1
  possible_ana_a[animal_ind] = 1
  
  plant_ana_rate =  (laa_par[1] + 
                       laa_par[2] *  abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0) %*% a_status[1:NCOL(M0)]) *
    p_status[1:NROW(M0)] * possible_ana_p
  animal_ana_rate =  (laa_par[3] + 
                        laa_par[4] * t(abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0)) %*% p_status[1:NROW(M0)]) * 
    a_status[1:NCOL(M0)] * possible_ana_a
  
  ana_list <- list(plant_ana_rate = plant_ana_rate,
                   animal_ana_rate = animal_ana_rate)
  return(ana_list)
}

#cladogenetic  rate
get_clado_rate <- function(lac0_par,
                           K_par,
                           Mt,
                           M0,
                           p_status,
                           a_status){
  
  NK_list <- get_NK(K_par=K_par,
                    M0=M0,
                    Mt=Mt,
                    p_status= p_status,
                    a_status= a_status)
  
  plant_clado_rate <- lac0_par[1] * NK_list[[3]] * p_status
  animal_clado_rate <- lac0_par[2] * NK_list[[4]] *a_status
  
  
  clado_list <- list(plant_clado_rate = plant_clado_rate,
                     animal_clado_rate = animal_clado_rate)
  
  return(clado_list)
}


#cospeciation rate
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

# gam0_par <- c(0.6,0.6)
# lac0_par <- c(00.3,0.3)
# mu_par <- c(0.02,0.01,0.02,0.01)
# laa_par <- c(0.1,0.2,0.1,0.2)
# lambda1 <- 0.1
# K_par <- c(20,0.6,20,0.6)
# M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# qloss <- 0.1
# qgain <- 0.1
# 

rates <-update_rates_mutual(gam0_par,lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain,M0,
               Mt,p_status,a_status,island_spec_plant,island_spec_animal)

update_rates_mutual <- function(gam0_par,
                                lac0_par,
                                K_par,
                                mu_par,
                                laa_par,
                                lambda1,
                                qgain,
                                qloss,
                                M0,
                                Mt,
                                p_status,
                                a_status,
                                island_spec_plant,
                                island_spec_animal) {
  
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
  
  #immigration rate
  immig_rate <- get_immig_rate(
    gam0_par = gam0_par,
    K_par = K_par,
    Mt = Mt,
    M0 = M0,
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
    a_status = a_status,
    island_spec_plant = island_spec_plant,
    island_spec_animal = island_spec_animal)
  #testit::assert(is.list(ana_rate))
  
  #cladogenetic rate
  clado_rate <- get_clado_rate(
    lac0_par = lac0_par,
    K_par = K_par,
    Mt = Mt,
    M0 =M0,
    p_status = p_status,
    a_status = a_status)
  
  #cospeciation rate
  cospec_rate <- get_cospec_rate(
    lambda1 = lambda1,
    K_par = K_par,
    M0 = M0,
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
    p_immig_rate = immig_rate$plant_immig_rate,
    p_ext_rate = ext_rate$plant_ext_rate,
    p_clado_rate = clado_rate$plant_clado_rate,
    p_ana_rate = ana_rate$plant_ana_rate,
    a_immig_rate = immig_rate$animal_immig_rate,
    a_ext_rate = ext_rate$animal_ext_rate,
    a_clado_rate = clado_rate$animal_clado_rate,
    a_ana_rate = ana_rate$animal_ana_rate,
    cospec_rate = cospec_rate,
    gain_rate = gain_rate,
    loss_rate = loss_rate)
  
  return(rates)
}

# save(rates, file="rates.Rdata")
