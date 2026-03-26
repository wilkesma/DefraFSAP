source(here::here("paths.R"))

####Packages
library(INLA)
library(parallel)

mod_args <- read.csv(file.path(PATH_PROCESSED, "inla_mod_args.csv", row.names = 1))

get_post <- function(mod) {
  inla.posterior.sample(n = 200, result = mod, num.threads=1)
}

run_posts <- function(task_id){
  message(task_id)
  mod_file <- file.path(PATH_PROCESSED, "inla_mods", "mods", paste0("mod_", task_id, ".rds"))
  
  no_metals <- NULL
  metals <- NULL
  
  if(file.exists(mod_file)) {
    mods <- readRDS(mod_file)
    
    if(length(mods[[1]])>1){
      if(class(mods[[1]]$fit)=="inla"){
        no_metals <- get_post(mods[[1]]$fit)
      }
    }
    
    if(length(mods[[2]])>1){
      if(class(mods[[2]]$fit)=="inla"){
        metals <- get_post(mods[[2]]$fit)
      }
    }
  }
  
  saveRDS(list(no_metals=no_metals, metals=metals), file.path(PATH_PROCESSED, "inla_mods", "posts", paste0(task_id, ".rds")))
}

mclapply(mod_args$task_id, run_posts, mc.cores=133)
