# get mutualistic partners and competitors for animal and plant species
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(1,1,0,1)
# a_status<-c(1,1,1,0,1)
# get_part_compe(Mt, p_status,a_status)

get_part_compe <- function(Mt,
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
  plant_compe[which(p_status==0)] <- 0
  # the number of competitors of each plant species, elements equal to 0 means no competitors or
  # the species is not shown up in island
  
  for (x in seq(ncol(Mt))){
    animal_compe[x] <- sum((colSums(Mt[,x] * Mt[,-x]* p_status)>=1) * as.matrix(a_status)[-x,])
  } 
  animal_compe[which(a_status==0)] <- 0
  # the number of competitors of each animal species,elements equal to 0 means no competitors or
  #the species is not shown up in island
  
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
get_NK <- function(K_par,
                   Mt,
                   p_status,
                   a_status){
  
  part_compe_list <- get_part_compe(Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
 
  part_compe_list[[1]][-which(part_compe_list[[1]]==0)]
  part_compe_list[[2]][-which(part_compe_list[[2]]==0)]
  
  plant_NK <-  exp(-(sum(p_status)/K_par[1])+
                     part_compe_list[[3]]/(K_par[2]*part_compe_list[[1]]))
  animal_NK <- exp(-(sum(a_status)/K_par[3])+
                     part_compe_list[[4]]/(K_par[4]*part_compe_list[[2]]))
  
  NK_list <- list(plant_NK = plant_NK,
                  animal_NK = animal_NK)
  
  return(NK_list)
}
# get p_status and a_status expanded
get_expand_matrix <- function(Mt,
                              p_status,
                              a_status){
  
  expd_p_sta <- matrix(rep(p_status,NCOL(Mt)), ncol = NCOL(Mt))
  expd_a_sta <- t(matrix(rep(a_status,NROW(Mt)),ncol = NROW(Mt)))
  
  expand_matrix_list <- list(expd_p_sta = expd_p_sta,
                             expd_a_sta = expd_a_sta)
  return(expand_matrix_list)
}

