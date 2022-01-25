
# Calculate immigration rate #
# mutualism_pars <- list(lac_par,mu_par,K_par,gam_par,laa_par,M0,qgain,qloss,lambda1)
# gam_par <- c(gam_plant, gam_animal)
# gam_plant: the initial immigration rate of plant species without mutualistic partners 
# gam_animal: the initial immigration rate of animal species without mutualistic partners 
# example: 
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# gam_par <- c(0.6,0.6)
# immig_list <- get_immig_rate(Mt,p_status,a_status,mutualism_pars)
get_immig_rate <- function(Mt,
                           p_status,
                           a_status,
                           mutualism_pars){
  
  NK_list <- get_NK(Mt = Mt,
                    p_status = p_status,
                    a_status = a_status)
  
  M0 <- mutualism_pars$M0
  gam_par <- mutualism_pars$gam_par
  
  plant_immig_rate <- gam_par[1] * NK_list[[1]][1:nrow(M0)] 
  animal_immig_rate <- gam_par[2] * NK_list[[2]][1:ncol(M0)] 
  
  
  immig_list <- list(plant_immig_rate = plant_immig_rate,
                           animal_immig_rate = animal_immig_rate)
                          
  return(immig_list)
}


