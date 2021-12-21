
# Calculate immigration rate #
# gam0_lac0_par includes: gam0_lac0_par <- c(gamP0,gamA0,lacP0,laca0)
# gamP0: the initial immigration rate of plant species without mutualistic partners 
# gamA0: the initial immigration rate of animal species without mutualistic partners 
# lacP0: the initial cladogenesis rate of plant species without mutualistic partners
# lacA0: the initial cladogenesis rate of animal species without mutualistic partners
# K_par includes: K_par <- c(KP0,KP1,KA0,KA1)
# KP0: the initial carrying capacity of plant species without any mutualists
# KA0: the initial carrying capacity of animal species without any mutualists
# KP1: a coefficient showing the influence from mutualism to plant species
# KA1: a coefficient showing the influence from mutualism to animal species
# Mt: interaction matrix on island at time t
# p_status: to show whether the plant species present on island
# a_status: to show whether the animal species present on island
#
# example: gam0_lac0_par <- c(gamP0,gamA0,lacP0,laca0)

# gam0_par <- c(0.6,0.6)
# K_par <- c(20,0.6,20,0.6)
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
#
#get_immig_rate(gam0_lac0_par,K_par,Mt,M0,p_status,a_status)
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


