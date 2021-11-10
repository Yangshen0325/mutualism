# get mutualistic partners and competitors for animal and plant species
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(1,1,1,1)
# a_status<-c(1,1,1,1,1)
# get_part_compe(Mt, p_status,a_status)

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
  
  plant_NK <-  exp(-(sum(p_status)/K_par[1])+
                     part_compe_list[[3]]/(K_par[2]*part_compe_list[[1]]))
  animal_NK <- exp(-(sum(a_status)/K_par[3])+
                     part_compe_list[[4]]/(K_par[4]*part_compe_list[[2]]))
  
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

