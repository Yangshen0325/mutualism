get_part_compe2 <- function(Mt,
                           p_status,
                           a_status){
  tMt <- t(Mt) 
  
  plant_part <- Mt %*% a_status #the number of partners that each plant species has
  animal_part <- tMt %*% p_status #the number of partners that each animal species has
  
  plant_compe <- matrix()
  animal_compe <- matrix()
  for (x in seq(ncol(tMt))){
    plant_compe[x] <- sum((colSums(tMt[,x]*tMt[,-x]* a_status)>=1) * as.matrix(p_status)[-x,])
  }  
  # the number of competitors of each plant species
  
  for (x in seq(ncol(Mt))){
    animal_compe[x] <- sum((colSums(Mt[,x] * Mt[,-x]* p_status)>=1) * as.matrix(a_status)[-x,])
  } 
  # the number of competitors of each animal species
  
  part_compe_list <- list(plant_part = plant_part,
                          animal_part = animal_part,
                          plant_compe = plant_compe,
                          animal_compe = animal_compe)
  
  return(part_compe_list)
}
# N/K
# K_par includes: K_par <- c(KP0,KP1,KA0,KA1)
# KP0: the initial carrying capacity of plant species without any mutualists
# KA0: the initial carrying capacity of animal species without any mutualists
# KP1: a coefficient showing the influence from mutualism to plant species
# KA1: a coefficient showing the influence from mutualism to animal species
# K_par <- c(20,0.6,20,0.6)
# get_NK(K_par,M0, p_status,a_status)
get_NK2 <- function(K_par,
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