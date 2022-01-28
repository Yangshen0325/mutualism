# get mutualistic partners and competitors for animal and plant species
# mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
# example:
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# part_compe_list <-get_part_compe(Mt, p_status,a_status)

get_part_compe <- function(Mt,
                           p_status,
                           a_status){
  tMt <- t(Mt)
  
  plant_partners <- Mt %*% a_status #the number of partners that each plant species has
  animal_partners <- tMt %*% p_status #the number of partners that each animal species has
  
  plant_compe <- matrix()
  animal_compe <- matrix()
  
  for (x in seq(nrow(Mt))){
    plant_compe[x] <- sum((colSums(tMt[,x]*tMt[,-x]* a_status)>=1) * as.matrix(p_status)[-x,])
  }  # the number of competitors of each plant species
  
  
  for (x in seq(ncol(Mt))){ 
    animal_compe[x] <- sum((colSums(Mt[,x] * Mt[,-x]* p_status)>=1) * as.matrix(a_status)[-x,])
  } # the number of competitors of each animal species
  
  part_compe_list <- list(plant_partners = plant_partners,     
                          animal_partners = animal_partners,   
                          plant_compe = plant_compe,   
                          animal_compe = animal_compe) 

  return(part_compe_list)
}
# N/K
# mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
# K_par <- c(K_A0, K_P1, K_A1)
# K_A0: the initial carrying capacity of animal species without any mutualists
# K_P1: a coefficient showing the influence from mutualism to plant species
# K_A1: a coefficient showing the influence from mutualism to animal species
# example:
# K <- Inf;
# K_par <- c(Inf,0.6,0.6)
# NK_list <- get_NK(K,Mt, p_status,a_status,mutualism_pars)
get_NK <- function(K,
                   Mt,
                   p_status,
                   a_status,
                   mutualism_pars){
  part_compe_list <- get_part_compe(Mt = Mt,
                                    p_status = p_status,
                                    a_status = a_status)
  
  indp <- which(part_compe_list[[1]]==0)
  inda <- which(part_compe_list[[2]]==0)
  
  K_par <- mutualism_pars$K_par
  
      plant_NK <- exp(-(sum(p_status)/K)+
                           part_compe_list[[3]]/(K_par[2]*part_compe_list[[1]]))
      plant_NK[indp] <- exp(-(sum(p_status)/K)) # N/K for plant species
 
      animal_NK <- exp(-(sum(a_status)/K_par[1])+
                            part_compe_list[[4]]/(K_par[3]*part_compe_list[[2]])) 
      animal_NK[inda] <- exp(-(sum(a_status)/K_par[1]))
  
  NK_list <- list(plant_NK = plant_NK,
                  animal_NK = animal_NK)
  
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

