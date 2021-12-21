# Calculate anagenesis rates
# 
# laa_par includes: laa_par <- c(p_laa0,p_laa1,a_laa0,a_laa1)
# p_laa0: the initial anagenesis rate of plant species without mutualistic partners
# a_laa0: the initial anagenesis rate of animal species without mutualistic partners
# p_laa1: a coefficient mediate the effects from mutualism
# a_laa1: a coefficient mediate the effects from mutualism
# M0: the initial matrix in the mainland
# Mt: the possible interaction matrix on island at time t
# p_status: to show whether the plant species present on island
# a_status: to show whether the animal species present on island
#
# example:
# island_spec_plant <- matrix(ncol = 7)
# island_spec_plant[1,] <- c(2,2,0,"I",NA,NA,NA)
# island_spec_plant <- rbind(island_spec_plant, c(3,3,0,"A",NA,NA,NA))
# 
# island_spec_animal <- matrix(ncol = 7)
# island_spec_animal[1,] <- c(1,1,0,"I",NA,NA,NA)
# island_spec_animal <- rbind(island_spec_animal,c(5,5,0,"A",NA,NA,NA))
#
# laa_par <- c(0.1,0.2,0.1,0.2)
# M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# 
# get_ana_rate(laa_par,M0,Mt,p_status,a_status,island_spec_plant,island_spec_animal)

get_ana_rate <- function(laa_par,
                         M0,
                         Mt,
                         p_status,
                         a_status,
                         island_spec_plant,
                         island_spec_animal) { 
  
  possible_ana_p <- matrix(0,nrow=NROW(M0))
  possible_ana_a <- matrix(0,nrow=NCOL(M0))
  
  plant_ind <-as.numeric(island_spec_plant[which(island_spec_plant[,4]=="I"),1])
  animal_ind <-as.numeric(island_spec_animal[which(island_spec_animal[,4]=="I"),1])
  possible_ana_p[plant_ind] = 1
  possible_ana_a[animal_ind] = 1
  
  plant_ana_rate =  (laa_par[1] + 
                       laa_par[2] *  abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0) %*% a_status[1:NCOL(M0)]) *
    p_status[1:NROW(M0)] * possible_ana_p
  animal_ana_rate =  (laa_par[3] + 
                        laa_par[4] * t(abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0)) %*% p_status[1:NROW(M0)]) * 
    a_status[1:NCOL(M0)] * possible_ana_a
  
  ana_list <- list(plant_ana_rate = plant_ana_rate,
                   animal_ana_rate = animal_ana_rate)
  return(ana_list)
}

