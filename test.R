# example:
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# part_compe_list <-get_part_compe(Mt, p_status,a_status)


# K_par <- c(K_P0, K_A0, K_P1, K_A1)
# K_A0: the initial carrying capacity of animal species without any mutualists
# K_P1: a coefficient showing the influence from mutualism to plant species
# K_A1: a coefficient showing the influence from mutualism to animal species
# example:
# K <- Inf;
# K_par <- c(inf,Inf,Inf,Inf)
# NK_list <- get_NK(K,Mt, p_status,a_status,mutualism_pars)

# get_expand_matrix <- (Mt, p_status, a_status)

#' mutualism_pars <- list(lac_par=lac_par,mu_par=mu_par,K_par=K_par,gam_par=gam_par,laa_par=laa_par,qgain=qgain,qloss=qloss,lambda1=lambda1,M0=M0,pro=pro)
# lac_par <- c(0.5, 0)           
# mu_par <- c(0.02, 0, 0, 0)          
# K_par <- c(Inf, Inf, Inf, Inf)           
# gam_par <- c(0.05, 0)          
# laa_par <- c(1, 0, 0, 0)
# qgain <- 0
# qloss <- 0
# lambda1 <- 0
# M0 <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# Mt <- {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# pro=1
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# update_rates_mutualism(Mt=Mt,p_status=p_status,a_status=a_status,mutualism_pars=mutualism_pars)

#example:
# library(reshape2)
# load("rates.Rdata")
# update_rates_mutual(gam0_lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain,
#                   M0,Mt,p_status,a_status)
# event <- possible_event(rates) 



# rates <- update_rates_mutual(gam0_par,mu_par,laa_par,lac0_par,lambda1,K_par,w,
#                   M0,Mt,p_status,a_status,qloss,qgain,num_immigrants)
# timeval <- 0
# calc_next_timeval(rates,timeval)
# update_rates_mutual(gam0_lac0_par,K_par,mu_par,laa_par,lambda1,qloss,qgain,
#                   M0,Mt,p_status,a_status)



# M0= {set.seed(1);matrix(sample(c(0,1),20,replace = TRUE),ncol=5,nrow=4)}
# p_status<-c(0,1,1,0)
# a_status<-c(1,0,0,0,1)
# possible_event <- possible_event(rates)  #could be like   Var1 Var2 value L1
###########################################################  1    1   0.7  8

######################################################
# island_spec_plant <- matrix(ncol = 7)
# island_spec_plant[1,] <- c(2,2,0,"I",NA,NA,NA)
# island_spec_plant <- rbind(island_spec_plant, c(3,3,0,"A",NA,NA,NA))
# 
# island_spec_animal <- matrix(ncol = 7)
# island_spec_animal[1,] <- c(1,1,0,"I",NA,NA,NA)
# island_spec_animal <- rbind(island_spec_animal,c(5,5,0,"A",NA,NA,NA))

# stt_table <- matrix(ncol = 7)
# colnames(stt_table) <- c("Time","nIp","nAp","nCp","nIa","nAa","nCa")
# stt_table[1,] <- c(totaltime,1,1,0,1,1,0)






















