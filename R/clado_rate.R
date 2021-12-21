##############
# lac0_par <- c(0.3,0.3)
#get_clado_rate(gam0_lac0_par,K_par,Mt,p_status,a_status)

get_clado_rate <- function(lac0_par,
                           K_par,
                           Mt,
                           M0,
                           p_status,
                           a_status){
  
  NK_list <- get_NK(K_par=K_par,
                    M0=M0,
                     Mt=Mt,
                     p_status= p_status,
                     a_status= a_status)
  
  plant_clado_rate <- lac0_par[1] * NK_list[[3]] * p_status
  animal_clado_rate <- lac0_par[2] * NK_list[[4]] *a_status
  
  
  clado_list <- list(plant_clado_rate = plant_clado_rate,
                     animal_clado_rate = animal_clado_rate)
  
  return(clado_list)
}