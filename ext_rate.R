
# Calculate extinction rate
# 
# mu_par includes: mu_par <- c(muP0,muP1,muA0,muA1)
# muP0: the initial extinction rate of plant species without mutualistic partners
# muP1: a coefficient mediate mutualism
# muA0: the initial extinction rate of animal species without mutualistic partners
# muA1: a coefficient mediate mutualism
# Mt: the possible interaction matrix on island at time t
# p_status: to show whether the plant species present on island
# a_status: to show whether the animal species present on island
#
# example:
#
# mu_par <- c(0.02,0.01,0.02,0.01)
# 
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(1,0,1,0)
# a_status<-c(1,1,0,0,1)
#
# get_ext_rate(mu_par,Mt,p_status,a_status)

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
# plant_ext_rate is a list, so is animal_ext_rate