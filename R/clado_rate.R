##############
# mutualism_pars <- list(lac_par,mu_par,K_par,gam_par,laa_par,M0,qgain,qloss,lambda1)
# lac_par <- c(lac_plant, lac_animal)
# lac_plant: the initial cladogenesis rates for plant species without mutualism
# lac_animal: the initial cladogenesis rates for animal species without mutualism
# example:
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# lac_par <- c(0.3,0.3)
# clado_list <- get_clado_rate(Mt,p_status,a_status,mutualism_pars)

get_clado_rate <- function(Mt,
                           p_status,
                           a_status,
                           mutualism_pars){
 
  NK_list <- get_NK(Mt = Mt,
                    p_status = p_status,
                    a_status = a_status)
  
  lac_par <- mutualism_pars$lac_par
  
  plant_clado_rate <- lac_par[1] * NK_list[[1]] * p_status 
  animal_clado_rate <- lac_par[2] * NK_list[[2]] *a_status
  #only species presented on island have cladogenestic rates
  
  clado_list <- list(plant_clado_rate = plant_clado_rate,
                     animal_clado_rate = animal_clado_rate)
  
  return(clado_list)
}