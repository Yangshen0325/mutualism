
# Calculate immigration rate #
# mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
# mu_par <- c(mu_A0,mu_P1,mu_A1)
# K_par <- c(K_A0, K_P1, K_A1)
# laa_par <- c(laa_A0,laa_P1, laa_A1)

# gam_animal: the initial immigration rate of animal species without mutualistic partners 
# example: 
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# gam_par <- c(0.6,0.6)
# immig_list <- get_immig_rate(Mt,p_status,a_status,mutualism_pars)
get_immig_rate <- function(gam,
                           K,
                           M0,
                           Mt,
                           p_status,
                           a_status,
                           mutualism_pars){
  
  if (is.null(mutualism_pars)){
    mainland_n <- sum(NROW(M0)+NCOL(M0))
    num_spec <- sum(p_status)+sum(a_status)
    immig_rate <- max(c( mainland_n * gam * (1 - (num_spec / K)),
                         0), na.rm = TRUE)
    return(immig_rate)
  } else {
    NK_list <- get_NK(Mt = Mt,
                      p_status = p_status,
                      a_status = a_status)
    
    
    gam_par <- mutualism_pars$gam_par
    
    plant_immig_rate <- gam_par[1] * NK_list[[1]][1:nrow(M0)] 
    animal_immig_rate <- gam_par[2] * NK_list[[2]][1:ncol(M0)] 
    
    
    immig_list <- list(plant_immig_rate = plant_immig_rate,
                       animal_immig_rate = animal_immig_rate)
    
    return(immig_list)
  }
}


