
# Calculate extinction rates
# mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
# mu_par <- c(mu_A0, mu_P1, mu_A1)
# mu_P1: a coefficient mediate mutualism
# mu_A0: the initial extinction rate of animal species without mutualistic partners
# mu_A1: a coefficient mediate mutualism

# example:
# mu_par <- c(0.2,0.01,0.01)
# ext_list <- get_ext_rate(mu=0.2,
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)},
#                      p_status<-c(0,1,1,0),
#                      a_status<-c(1,0,0,0,1),
#                        mutualism_pars = mutualism_pars)

get_ext_rate <- function(mu,
                         Mt,
                         p_status,
                         a_status,
                         mutualism_pars) {
  
  if (is.null(mutualism_pars)){
    num_spec <- sum(p_status)+sum(a_status)
    ext_rate <- max(0, mu * num_spec, na.rm = TRUE)
    return(ext_rate)
  } else {
  part_compe_list <- get_part_compe(Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
  mu_par <- mutualism_pars$mu_par
  
  plant_ext_rate <- as.matrix(pmax(0, mu - mu_par[2] * part_compe_list[[1]])) * p_status
  animal_ext_rate <- as.matrix(pmax(0, mu_par[1] - mu_par[3] * part_compe_list[[2]])) *a_status
  
  ext_list <- list(plant_ext_rate = plant_ext_rate,
                   animal_ext_rate = animal_ext_rate)
  return(ext_list)
  }
}
# plant_ext_rate is a list, so is animal_ext_rate