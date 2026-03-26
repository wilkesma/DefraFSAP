source(here::here("paths.R"))

##Species names
get_species <- function(group){
  load(file.path(PATH_PROCESSED, paste0(group, "_data.RData")))
  if(group=="invs"){
    invs_abun.sc2[invs_abun.sc2>0] <- 1 #Convert to binary just to calculate prevalence
    data.frame(group=group, species=colnames(invs_abun.sc2), prev=colSums(invs_abun.sc2)/nrow(invs_abun.sc2))
  } else{
    fish_abun[fish_abun>0] <- 1 #Convert to binary just to calculate prevalence
    data.frame(group=group, species=colnames(fish_abun), prev=colSums(fish_abun)/nrow(fish_abun))
  }
}
species_names <- rbind(get_species("invs"), get_species("fish"))
species_names$minor <- FALSE
species_names$minor[species_names$species %in% c("Barbatula barbatula", "Cobitis taenia", "Cottus gobio", "Gasterosteus aculeatus", "Phoxinus phoxinus", "Pungitius pungitius")] <- TRUE

mod_args <- data.frame(task_id=1:nrow(species_names), species_names)

write.csv(mod_args, file.path(PATH_PROCESSED, "inla_mod_args.csv"))

#Prepare subdirectories for modelling
dir.create(file.path(PATH_PROCESSED, "inla_mods", "mods"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PATH_PROCESSED, "inla_mods", "posts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PATH_PROCESSED, "inla_mods", "logs"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PATH_PROCESSED, "inla_future_pred", "preds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PATH_PROCESSED, "inla_future_pred", "logs"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PATH_PROCESSED, "inla_spatial_pred", "preds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PATH_PROCESSED, "inla_spatial_pred", "logs"),  recursive = TRUE, showWarnings = FALSE)
