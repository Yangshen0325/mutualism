DAISIE_create_island_mutualism <- function(stt_table,
                                           totoaltime,
                                           island_spec,
                                           M0){
  ### if there are no species on the island branching_times = island_age,
  ### stac = 0, missing_species = 0
  mainland_n <- NROW(M0) + NCOL(M0)
  
  if(length(island_spec[, 1]) == 0){
    island <- list(stt_table = stt_table,
                   branching_times = totoaltime,
                   stac = 0,
                   missing_species = 0)
  } else {
    cnames <- c("Species",
                "Mainland Ancestor",
                "Colonisation time (BP)",
                "Species type",
                "branch_code",
                "branching time (BP)",
                "Anagenetic_origin",
                "Species state" )
    colnames(island_spec) <- cnames
    ### set ages as counting backwards from present
    island_spec[, "branching time (BP)"] <- totaltime -
      as.numeric(island_spec[, "branc hing time (BP)"])
    island_spec[, "Colonisation time (BP)"] <- totaltime -
      as.numeric(island_spec[, "Colonisation time (BP)"])
     
    if (mainland_n == 1){
      island <- DAISIE_ONEcolonist(totaltime,
                                   island_spec,
                                   stt_table)
    } else if (mainland_n > 1){
     
    ### number of colonists present
    island_spec_plant <- island_spec[which(island_spec[,8] == "plant"), ]
    island_spec_animal <- island_spec[which(island_spec[,8] == "animal"), ]
    present_plant <- sort(as.numeric(unique(island_spec_plant[, "Mainland Ancestor"])))
    present_animal <- sort(as.numeric(unique(island_spec_animal[, "Mainland Ancestor"])))
    #num_present_plant <- length(present_plant) 
    #num_present_animal <- length(present_animal)
    
      island_clades_info_plant <- list()
      for (i in 1:num_present_plant) {
        subset_island_plant <- island_spec_plant[which(island_spec_plant[, "Mainland Ancestor"] ==
                                             present_plant[i]), ]
        if (!is.matrix(subset_island_plant)) {
          subset_island_plant <- rbind(subset_island_plant[1:8])
          colnames(subset_island_plant) <- cnames
        }
        island_clades_info_plant[[i]] <- DAISIE_ONEcolonist(
          totaltime,
          island_spec = subset_island_plant,
          stt_table = NULL)
        island_clades_info_plant[[i]]$stt_table <- NULL
      }
      
      island_clades_info_animal <- list()
      for (i in 1:num_present_animal) {
        subset_island_animal <- island_spec_animal[which(island_spec_animal[, "Mainland Ancestor"] ==
                                                         present_animal[i]), ]
        if (!is.matrix(subset_island_animal)) {
          subset_island_animal <- rbind(subset_island_animal[1:8])
          colnames(subset_island_animal) <- cnames
        }
        island_clades_info_animal[[i]] <- DAISIE_ONEcolonist(
          totaltime,
          island_spec = subset_island_animal,
          stt_table = NULL)
        island_clades_info_animal[[i]]$stt_table <- NULL
      }
      island_clades_info <- c(island_clades_info_plant,island_clades_info_animal)  
      island <- list(stt_table = stt_table,
                     taxon_list = island_clades_info)
    }
  }
  return(island)
}

