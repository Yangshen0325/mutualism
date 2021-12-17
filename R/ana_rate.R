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
#
# laa_par <- c(0.1,0.2,0.1,0.2)
# M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,1,1)
# 
# rep.row<-function(x,n){
#   matrix(rep(x,each=n),nrow=n)}
# island_spec_plant <-c("NA","NA","NA","I","NA","NA","NA")
# island_spec_plant <- rep.row(c(island_spec_plant),4)
# island_spec_plant[c(1,4),4] <- "A"
# 
# island_spec_animal <-c("NA","NA","NA","I","NA","NA","NA")
# island_spec_animal <- rep.row(c(island_spec_animal),5)
# island_spec_animal[c(2,3,4),4] <- "A"

# get_ana_rate(laa_par,M0,Mt,p_status,a_status)

get_ana_rate <- function(laa_par,
                         M0,
                         Mt,
                         p_status,
                         a_status,
                         island_spec_plant,
                         island_spec_animal) {  
  
  plant_ana_rate =  (laa_par[1] + 
                       laa_par[2] *  abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0) %*% a_status[1:NCOL(M0)]) *
                        p_status * (island_spec_plant[, 4] == "I")
  animal_ana_rate =  (laa_par[3] + 
                        laa_par[4] * t(abs(Mt[1:NROW(M0),1:NCOL(M0)]-M0)) %*% p_status[1:NROW(M0)]) * 
                          a_status * (island_spec_animal[, 4] == "I")
  
  ana_list <- list(plant_ana_rate = plant_ana_rate,
                   animal_ana_rate = animal_ana_rate)
  return(ana_list)
}

