##############
# mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
# lac_animal: the initial cladogenesis rates for animal species without mutualism

# example:
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# lac_par <- c(0.3,0.3)
# clado_list <- get_clado_rate(lac=0.3,
#                              K = Inf,
# M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)},  
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)},
#                    p_status<-c(0,1,1,0),
#                    a_status<-c(1,0,0,0,1),
#                          mutualism_pars = mutualism_pars)

get_clado_rate <- function(lac,
                           K,
                           Mt,
                           p_status,
                           a_status,
                           mutualism_pars){
  
  if (is.null(mutualism_pars)){
    num_spec <- sum(p_status)+sum(a_status)
    clado_rate <- max(
      0, lac * num_spec  * (1 - num_spec /K), na.rm = TRUE)
    return(clado_rate)
  }else{
    
    NK_list <- get_NK(K = K,
                      Mt = Mt,
                      p_status = p_status,
                      a_status = a_status,
                      mutualism_pars = mutualism_pars)
    
    lac_animal <- mutualism_pars$lac_animal
    
    plant_clado_rate <- lac * NK_list[[1]] * p_status # "* p_status" assures the                                                           # species is on the island 
    animal_clado_rate <- lac_animal * NK_list[[2]] *a_status
    #only species presented on island have cladogenestic rates
    
    clado_list <- list(plant_clado_rate = plant_clado_rate,
                       animal_clado_rate = animal_clado_rate)
    
    return(clado_list)
  }
}