# Calculate anagenesis rates
# mutualism_pars <- list(lac_animal,mu_par,K_par,gam_animal,laa_par,qgain,qloss,lambda1)
# laa_par <- c(laa_A0, laa_P1, laa_A1)
# laa_A0: the initial anagenesis rate of animal species without mutualistic partners
# laa_P1: a coefficient mediate the effects from mutualism
# laa_A1: a coefficient mediate the effects from mutualism

# example:
# island_spec_plant <- matrix(ncol = 7)
# island_spec_plant[1,] <- c(2,2,0,"I",NA,NA,NA)
# island_spec_plant <- rbind(island_spec_plant, c(3,3,0,"A",NA,NA,NA))
# 
# island_spec_animal <- matrix(ncol = 7)
# island_spec_animal[1,] <- c(1,1,0,"I",NA,NA,NA)
# island_spec_animal <- rbind(island_spec_animal,c(5,5,0,"A",NA,NA,NA))
#
# laa_par <- c(0.1,0.01,0.01)
# ana_list <- get_ana_rate(laa=0.1,
#                        num_immigrant=2,
# M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)},
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)},
#                      p_status<-c(0,1,1,0),
#                      a_status<-c(1,0,0,0,1),
#                      mutualism_pars = mutualism_pars,
#                     island_spec_plant = island_spec_plant,
#                     island_spec_animal = island_spec_animal)
get_ana_rate <- function(laa,
                         num_immigrant,
                         M0,
                         Mt,
                         p_status,
                         a_status,
                         mutualism_pars) {
  
  if(is.null(mutualism_pars)){
    ana_rate <- laa * num_immigrants
    return(ana_rate)
  }else{
    
    laa_par <- mutualism_pars$laa_par
    
    plant_ana_rate = laa + laa_par[2] * (abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0) %*% a_status[1:NCOL(M0)]) 
    # "[1:NROW(M0),1:NCOL(M0)]" keeps only species ID from mainland can happend to anagenesis.                                                                                     
    animal_ana_rate =  laa_par[1] + laa_par[3] * (t(abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0)) %*%p_status[1:NROW(M0)]) 
    
    ana_list <- list(plant_ana_rate = plant_ana_rate,
                     animal_ana_rate = animal_ana_rate)
    return(ana_list)
  }
}

