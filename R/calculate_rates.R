# Calculates algorithm rates
update_rates_mutualism <- function(Mt,
                                   p_status,
                                   a_status,
                                   mutualism_pars,
                                   island_spec){
  
  immig_rate <- get_immig_rate(
    Mt = Mt,
    p_status = p_status,
    a_status = a_status,
    mutualism_pars = mutualism_pars
  )
  
  ext_rate <- get_ext_rate(
    Mt = Mt,
    p_status = p_status,
    a_status = a_status
  )
  
  ana_rate <- get_ana_rate(
    Mt = Mt,
    p_status = p_status,
    a_status = a_status,
    mutualism_pars = mutualism_pars,
    island_spec = island_spec
  )
  
  clado_rate <- get_clado_rate(
    Mt = Mt,
    p_status = p_status,
    a_status = a_status,
    mutualism_pars = mutualism_pars
  )
  
  gain_rate <- get_gain_rate(Mt = Mt,
                             p_status = p_status,
                             a_status = a_status,
                             mutualism_pars = mutualism_pars)
  
  loss_rate <- get_loss_rate(Mt = Mt,
                             p_status = p_status,
                             a_status = a_status,
                             mutualism_pars = mutualism_pars)
  
  cospec_rate <- get_cospec_rate(Mt = Mt,
                                 p_status = p_status,
                                 a_status = a_status,
                                 mutualism_pars = mutualism_pars)
  
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
#immigration rate
get_immig_rate <- function(Mt,
                           p_status,
                           a_status,
                           mutualism_pars) {
  
  NK_list <- get_NK(Mt = Mt,
                    p_status = p_status,
                    a_status = a_status,
                    mutualism_pars = mutualism_pars)
  
  M0 <- mutualism_pars$M0
  gam_par <- mutualism_pars$gam_par
  
  plant_immig_rate <- gam_par[1] * NK_list[[1]]
  plant_immig_rate <- as.matrix(plant_immig_rate[1:NROW(M0)])
  animal_immig_rate <- gam_par[2] * NK_list[[2]]
  animal_immig_rate <- as.matrix(animal_immig_rate[1:NCOL(M0)])
  
  immig_list <- list(plant_immig_rate = plant_immig_rate,
                     animal_immig_rate = animal_immig_rate)
  return(immig_list)
}

#extinction rate
get_ext_rate <- function(Mt,
                         p_status,
                         a_status) {
  
  part_compe_list <- get_part_compe(Mt = Mt,
                                    p_status = p_status,
                                    a_status = a_status)
  mu_par <- mutualism_pars$mu_par
  
  plant_ext_rate <- as.matrix(pmax(0, mu_par[1] - mu_par[3] * part_compe_list[[1]])) * p_status
  animal_ext_rate <- as.matrix(pmax(0, mu_par[2] - mu_par[4] * part_compe_list[[2]])) *a_status
  
  ext_list <- list(plant_ext_rate = plant_ext_rate,
                   animal_ext_rate = animal_ext_rate)
  return(ext_list)
}

#anagenetic rate
get_ana_rate <- function(Mt,
                         p_status,
                         a_status,
                         mutualism_pars,
                         island_spec) {
  
  laa_par <- mutualism_pars$laa_par
  M0 <- mutualism_pars$M0
  
  possible_ana_p <- matrix(0,nrow=NROW(M0))
  possible_ana_a <- matrix(0,nrow=NCOL(M0))
  
  plant_ind <-as.numeric(island_spec[intersect(which(island_spec[,4] == "I"),
                                              which(island_spec[,8] == "plant")),1])
  animal_ind <-as.numeric(island_spec[intersect(which(island_spec[,4] == "I"),
                                    which(island_spec[,8] == "animal")),1])
  possible_ana_p[plant_ind] = 1
  possible_ana_a[animal_ind] = 1
  
  plant_ana_rate =  (laa_par[1] +
                       laa_par[3] * abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0) %*% a_status[1:NCOL(M0)]) *
    p_status[1:NROW(M0)] * possible_ana_p
  animal_ana_rate =  (laa_par[2] +
                        laa_par[4] * t(abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0)) %*% p_status[1:NROW(M0)]) *
    a_status[1:NCOL(M0)] * possible_ana_a
  
  ana_list <- list(plant_ana_rate = plant_ana_rate,
                   animal_ana_rate = animal_ana_rate)
  return(ana_list)
}

#cladogenetic  rate
get_clado_rate <- function(Mt,
                           p_status,
                           a_status,
                           mutualism_pars) {
  
  NK_list <- get_NK(Mt = Mt,
                    p_status = p_status,
                    a_status = a_status,
                    mutualism_pars = mutualism_pars)
  
  lac_par <- mutualism_pars$lac_par
  
  plant_clado_rate <- lac_par[1] * NK_list[[1]] * p_status
  animal_clado_rate <- lac_par[2] * NK_list[[2]] *a_status
  
  clado_list <- list(plant_clado_rate = plant_clado_rate,
                     animal_clado_rate = animal_clado_rate)
  return(clado_list)
}

#cospeciation rate
get_cospec_rate <- function(Mt,
                            p_status,
                            a_status,
                            mutualism_pars){
  
  NK_list <- get_NK(Mt = Mt,
                    p_status = p_status,
                    a_status = a_status,
                    mutualism_pars = mutualism_pars)
  
  expand_matrix_list <- get_expand_matrix(Mt=Mt,
                                          p_status = p_status,
                                          a_status = a_status)
  
  lambda1 <- mutualism_pars$lambda1
  
  cospec_rate <- lambda1 * Mt * expand_matrix_list[[1]] * expand_matrix_list[[2]] *
    matrix(rep(NK_list[[1]],NCOL(Mt)), ncol = NCOL(Mt)) *
    t(matrix(rep(NK_list[[2]],NROW(Mt)),ncol = NROW(Mt)))
  
  return(cospec_rate)
}

#gain rate
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

#loss rate
get_loss_rate <- function(Mt,
                          p_status,
                          a_status,
                          mutualism_pars){
  
  qloss <- mutualism_pars$qloss
  loss_rate <- qloss * Mt 
  
  return(loss_rate)
}


