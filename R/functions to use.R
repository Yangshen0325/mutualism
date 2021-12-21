# get mutualistic partners and competitors for animal and plant species
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# get_part_compe(M0,Mt, p_status,a_status)

get_part_compe <- function(M0,
                           Mt,
                           p_status,
                           a_status){
  
  tMtforimmig <- t(Mt[1:NROW(M0),1:NCOL(M0)]) 
  tMtforspecia <-t(Mt)
  
  plant_part_forimmig <- Mt[1:NROW(M0),1:NCOL(M0)] %*% a_status[1:NCOL(M0)] #the number of partners that each mainland plant species has
  animal_part_forimmig <- tMtforimmig %*% p_status[1:NROW(M0)] #the number of partners that each mainland animal species has
  
  plant_part_forspecia <- Mt %*% a_status
  animal_part_forspecia <- tMtforspecia %*% p_status
  
  plant_compe_forimmig <- matrix()
  animal_compe_forimmig <- matrix()
  
  plant_compe_forspecia <- matrix()
  animal_compe_forspecia <- matrix()
  
  for (x in seq(ncol(tMtforimmig))){
    plant_compe_forimmig[x] <- sum((colSums(tMtforimmig[,x]*tMtforimmig[,-x]* a_status[1:NCOL(M0)])>=1) * as.matrix(p_status[1:NROW(M0)])[-x,])
  }  
  # the number of competitors of each plant species that comes from mainland
  
  for (x in seq(ncol(tMtforspecia))){
    plant_compe_forspecia[x] <- sum((colSums(tMtforspecia[,x]*tMtforspecia[,-x]* a_status)>=1) * as.matrix(p_status)[-x,])
  } 
  # the number of competitors of each plant species that are shown on island
  
  for (x in seq(ncol(Mt[1:NROW(M0),1:NCOL(M0)]))){
    animal_compe_forimmig[x] <- sum((colSums(Mt[1:NROW(M0),1:NCOL(M0)][,x] * Mt[1:NROW(M0),1:NCOL(M0)][,-x]* p_status[1:NROW(M0)])>=1) * as.matrix(a_status[1:NCOL(M0)])[-x,])
  } 
  # the number of competitors of each animal species that comes from mainland
  for (x in seq(ncol(Mt))){
    animal_compe_forspecia[x] <- sum((colSums(Mt[,x] * Mt[,-x]* p_status)>=1) * as.matrix(a_status)[-x,])
  } 
  
  part_compe_list <- list(plant_part_forimmig = plant_part_forimmig,     #[[1]] plant partners used for immigration
                          animal_part_forimmig = animal_part_forimmig,   #[[2]] aniaml partners used for immigration
                          plant_part_forspecia = plant_part_forspecia,   #[[3]] plant partners used for speciation
                          animal_part_forspecia = animal_part_forspecia, #[[4]] aniaml partners used for speciation
                         
                          plant_compe_forimmig = plant_compe_forimmig,   #[[5]] plant competitors used for immigration
                          animal_compe_forimmig = animal_compe_forimmig, #[[6]] animal competitors used for immigration
                          plant_compe_forspecia = plant_compe_forspecia, #[[7]] plant competitors used for speciation
                          animal_compe_forspecia = animal_compe_forspecia)#[[8]] animal competitors used for speciation
  
  return(part_compe_list)
}
# N/K
# K_par includes: K_par <- c(KP0,KP1,KA0,KA1)
# KP0: the initial carrying capacity of plant species without any mutualists
# KA0: the initial carrying capacity of animal species without any mutualists
# KP1: a coefficient showing the influence from mutualism to plant species
# KA1: a coefficient showing the influence from mutualism to animal species
# K_par <- c(20,0.6,20,0.6)
# get_NK(K_par,M0,Mt, p_status,a_status)
get_NK <- function(K_par,
                   M0,
                   Mt,
                   p_status,
                   a_status){
  part_compe_list <- get_part_compe(M0=M0,
                                    Mt=Mt,
                                    p_status= p_status,
                                    a_status= a_status)
  
  indp_immig <- which(part_compe_list[[1]]==0)
  inda_immig <- which(part_compe_list[[2]]==0)
  
      plant_NK_forimmig <- exp(-(sum(p_status)/K_par[1])+
                           part_compe_list[[5]]/(K_par[2]*part_compe_list[[1]]))
      plant_NK_forimmig[indp_immig] <- exp(-(sum(p_status)/K_par[1]))
 
 
      animal_NK_forimmig <- exp(-(sum(a_status)/K_par[3])+
                            part_compe_list[[6]]/(K_par[4]*part_compe_list[[2]])) 
      animal_NK_forimmig[inda_immig] <- exp(-(sum(a_status)/K_par[3]))
      
      ##############
      
      indp_specia <- which(part_compe_list[[3]]==0)
      inda_specia <- which(part_compe_list[[4]]==0)
      
      plant_NK_forspecia <- exp(-(sum(p_status)/K_par[1])+
                                 part_compe_list[[7]]/(K_par[2]*part_compe_list[[3]]))
      plant_NK_forspecia[indp_specia] <- exp(-(sum(p_status)/K_par[1]))
      
      
      animal_NK_forspecia <- exp(-(sum(a_status)/K_par[3])+
                                  part_compe_list[[8]]/(K_par[4]*part_compe_list[[4]])) 
      animal_NK_forspecia[inda_specia] <- exp(-(sum(a_status)/K_par[3]))
  
  NK_list <- list(plant_NK_forimmig = plant_NK_forimmig,
                  animal_NK_forimmig = animal_NK_forimmig,
                  plant_NK_forspecia = plant_NK_forspecia,
                  animal_NK_forspecia = animal_NK_forspecia)
  
  return(NK_list)
}
# get p_status and a_status expanded
# get_expand_matrix <- (Mt, p_status, a_status)
get_expand_matrix <- function(Mt,
                              p_status,
                              a_status){
  expd_p_sta <- matrix(rep(p_status,NCOL(Mt)), ncol = NCOL(Mt))
  expd_a_sta <- t(matrix(rep(a_status,NROW(Mt)),ncol = NROW(Mt)))
  
  expand_matrix_list <- list(expd_p_sta = expd_p_sta,
                             expd_a_sta = expd_a_sta)
  return(expand_matrix_list)
}

# new Mt if cladogenesis happends with plant species
new_Mt_clado_p <- function(Mt, possible_event, p){
  
  possible_output <- list(c(1,1), c(1,0), c(0,1))
  newrows <- list()
  h <- possible_event$Var1
  newrows[c(which(Mt[h,]==0))] <- list(c(0,0))
  newrows[c(which(Mt[h,]==1))] <- sample(possible_output, 
                                         size=length(which(Mt[h,]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)/2,(1-p)/2))
  newrows<-matrix(unlist(newrows),nrow=2,ncol=NCOL(Mt))
  
  Mt <- rbind(Mt,newrows)
  return(Mt)
}

# new Mt if anagenesis happends with plant species
new_Mt_ana_p <- function(Mt, possible_event, p){
  
  newrows <- c()
  h <- possible_event$Var1
  newrows[c(which(Mt[h,]==0))] <- 0
  newrows[c(which(Mt[h,]==1))] <- sample(c(1,0), 
                                         size=length(which(Mt[h,]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)))
  Mt <- rbind(Mt,newrows)
  return(Mt)
}
# if cladogenesis happends with animal species
new_Mt_clado_a <- function(Mt, possible_event, p){
  
  possible_output <- list(c(1,1), c(1,0), c(0,1))
  newcols <- list()
  h <- possible_event$Var1
  newcols[c(which(Mt[,h]==0))] <- list(c(0,0))
  newcols[c(which(Mt[,h]==1))] <- sample(possible_output, 
                                         size=length(which(Mt[,h]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)/2,(1-p)/2))
  
  newcols<-t(matrix(unlist(newcols),nrow=2,ncol=NROW(Mt)))
  Mt <- cbind(Mt,newcols)
  return(Mt)
}

# if anagenesis happends with animal species
# Mt = {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
new_Mt_ana_a <- function(Mt, possible_event, p){  
  
  newcols <- c()
  h <- possible_event$Var1
  newcols[c(which(Mt[,h]==0))] <- 0
  newcols[c(which(Mt[,h]==1))] <- sample(c(1,0), 
                                         size=length(which(Mt[,h]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)))
  
  Mt <- cbind(Mt,newcols)
  return(Mt)
}

# if cospeciation happends with plant species i and animal species j
# Mt = {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
new_Mt_cospec <- function(Mt, possible_event, p){
  
  possible_output <- list(c(1,1), c(1,0), c(0,1))
  newrows <- list()
  newcols <- list()
  h <- possible_event$Var1
  k <- possible_event$Var2
  
  newcols[c(which(Mt[,k]==0))] <- list(c(0,0))
  newcols[c(which(Mt[,k]==1))] <- sample(possible_output, 
                                         size=length(which(Mt[,k]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)/2,(1-p)/2))
  
  newcols<-t(matrix(unlist(newcols),nrow=2,ncol=NROW(Mt)))
  
  newrows[c(which(Mt[h,]==0))] <- list(c(0,0))
  newrows[c(which(Mt[h,]==1))] <- sample(possible_output, 
                                         size=length(which(Mt[h,]==1)),
                                         replace = TRUE,
                                         prob=c(p,(1-p)/2,(1-p)/2))
  
  newrows <- matrix(unlist(newrows),nrow=2,ncol=NCOL(Mt))
  newrows <- cbind(newrows, diag(1,2,2))
  Mt <- cbind(Mt,newcols)
  Mt <- rbind(Mt,newrows)
  
  return(Mt)
}

