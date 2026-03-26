#!/usr/bin/Rscript
task_id <- as.integer(Sys.getenv("SGE_TASK_ID")) #If moving to grid job

source(here::here("paths.R"))

####Packages
library(INLA)
#inla.upgrade(testing=TRUE)
#inla.binary.install()
#1 #Interactively choose CentOS Linux-7
library(inlatools)
library(dplyr)

mod_args <- read.csv(file.path(PATH_PROCESSED, "inla_mod_args.csv"), row.names = 1)

##Predictor variables
vars <- c("PC1", "PC2", "sewage", "rainfall", "ASR", "CRI", "mining", "CAMS", "barrier_density",
          paste0(rep(c("pH", "TIN", "PO4", "wT", "Cu_d", "Zn_d"), each=4), c("_12M", "_9M", "_6M", "_3M")),
          "SITE_ID", "basin", "year", "julian", "season")

invs_env <- readRDS(file.path(PATH_PROCESSED, "invs_env_data_complete_2003_2023.rds"))
fish_env <- readRDS(file.path(PATH_PROCESSED, "fish_env_data_complete_2003_2023.rds"))

##Scaling
invs_env <- invs_env %>%
  mutate(
    PC1 = scale(PC1)[,1],
    PC2 = scale(PC2)[,1],
    CRI = scale(CRI)[,1],
    ASR = scale(ASR)[,1],
    wT_12M = scale(wT_12M)[,1],
    TIN_12M = scale(TIN_12M)[,1],
    PO4_12M = scale(PO4_12M)[,1],
    pH_12M = scale(pH_12M)[,1],
    Cu_d_12M = scale(Cu_d_12M)[,1],
    Zn_d_12M = scale(Zn_d_12M)[,1],
    wT_9M = scale(wT_9M)[,1],
    TIN_9M = scale(TIN_9M)[,1],
    PO4_9M = scale(PO4_9M)[,1],
    pH_9M = scale(pH_9M)[,1],
    Cu_d_9M = scale(Cu_d_9M)[,1],
    Zn_d_9M = scale(Zn_d_9M)[,1],
    wT_6M = scale(wT_6M)[,1],
    TIN_6M = scale(TIN_6M)[,1],
    PO4_6M = scale(PO4_6M)[,1],
    pH_6M = scale(pH_6M)[,1],
    Cu_d_6M = scale(Cu_d_6M)[,1],
    Zn_d_6M = scale(Zn_d_6M)[,1],
    wT_3M = scale(wT_3M)[,1],
    TIN_3M = scale(TIN_3M)[,1],
    PO4_3M = scale(PO4_3M)[,1],
    pH_3M = scale(pH_3M)[,1],
    Cu_d_3M = scale(Cu_d_3M)[,1],
    Zn_d_3M = scale(Zn_d_3M)[,1],
    rainfall = scale(rainfall)[,1]
  )

fish_env <- fish_env %>%
  mutate(
    PC1 = scale(PC1)[,1],
    PC2 = scale(PC2)[,1],
    CRI = scale(CRI)[,1],
    ASR = scale(ASR)[,1],
    wT_12M = scale(wT_12M)[,1],
    TIN_12M = scale(TIN_12M)[,1],
    PO4_12M = scale(PO4_12M)[,1],
    pH_12M = scale(pH_12M)[,1],
    Cu_d_12M = scale(Cu_d_12M)[,1],
    Zn_d_12M = scale(Zn_d_12M)[,1],
    wT_9M = scale(wT_9M)[,1],
    TIN_9M = scale(TIN_9M)[,1],
    PO4_9M = scale(PO4_9M)[,1],
    pH_9M = scale(pH_9M)[,1],
    Cu_d_9M = scale(Cu_d_9M)[,1],
    Zn_d_9M = scale(Zn_d_9M)[,1],
    wT_6M = scale(wT_6M)[,1],
    TIN_6M = scale(TIN_6M)[,1],
    PO4_6M = scale(PO4_6M)[,1],
    pH_6M = scale(pH_6M)[,1],
    Cu_d_6M = scale(Cu_d_6M)[,1],
    Zn_d_6M = scale(Zn_d_6M)[,1],
    wT_3M = scale(wT_3M)[,1],
    TIN_3M = scale(TIN_3M)[,1],
    PO4_3M = scale(PO4_3M)[,1],
    pH_3M = scale(pH_3M)[,1],
    Cu_d_3M = scale(Cu_d_3M)[,1],
    Zn_d_3M = scale(Zn_d_3M)[,1],
    rainfall = scale(rainfall)[,1]
  )

##Prep model data frame
load(file.path(PATH_PROCESSED, "invs_data.RData")) #invs_abun.sc2
rm(invs_info.sc2, invs_meta)
load(file.path(PATH_PROCESSED, "fish_data.RData")) #fish_abun
rm(fish_meta)

#e.g. dropvar=c("Cu_d", "Zn_d") if you want to model more data without metals
get_data <- function(species, group, dropvar=NULL){
  if(!is.null(dropvar)){
    vars_drop <- vars[-which(vars %in% c(dropvar, paste0(rep(dropvar,each=4), c("_12M", "_9M", "_6M", "_3M"))))]
  } else{
    vars_drop <- vars
  }
  if(group=="invs"){
    env <- invs_env[complete.cases(invs_env[, vars_drop]), c("ANALYSIS_ID", vars, "fold")]
    env$ANALYSIS_ID <- as.character(format(env$ANALYSIS_ID, scientific=FALSE)) #Was converting some long ANALYSIS IDs to scientific format
    env <- env[which(env$ANALYSIS_ID %in% row.names(invs_abun.sc2)),] #1 ANALYSIS ID wasn't in the abundance data
    env <- env[complete.cases(env[, vars_drop]), ]
    output <- data.frame(y=invs_abun.sc2[env$ANALYSIS_ID, species], env)
  } else{
    env <- fish_env[complete.cases(fish_env[, vars_drop]), c("SURVEY_ID", "FISHED_AREA", "REGION", vars, "fold")]
    env$SURVEY_ID <- as.character(format(env$SURVEY_ID, scientific=FALSE)) #Was converting some long ANALYSIS IDs to scientific format
    env <- env[which(env$SURVEY_ID %in% row.names(fish_abun)),]
    env <- env[complete.cases(env[, vars_drop]), ]
    output <- data.frame(y=fish_abun[env$SURVEY_ID, species], env)
  }
  output
}

do_cv <- function(data, species, group, minor, fold, form, fam="zeroinflatednbinomial1"){
  message(paste0("Fold ", fold))
  #Prepare prediction data
  pred_data <- data[data$fold == fold, ]
  pred_data$y <- NA
  
  combined_data <- bind_rows(
    data[data$fold != fold, ], #training data
    pred_data #prediction data with y = NA
  )
  
  #Fit model
  site_prior <- list(prec = list(prior = "pc.prec", param = c(3, 0.05))) #P(SD > 3) = 0.05
  fold_fit <- inla(as.formula(form), family = fam, data = combined_data, control.fixed=list(mean = 0, prec = 0.33), control.inla=list(cmin=0))
  
  #Extract fixed effects only predictions and convert to count scale
  n_train <- sum(data$fold != fold)
  predictions <- as.data.frame(apply(fold_fit$summary.linear.predictor[(n_train + 1):nrow(combined_data), ], 2, exp))
  
  #Get observed values
  observed <- data[data$fold == fold, "y"]
  if(group=="fish"){
    observed <- observed/data$FISHED_AREA[data$fold==fold]
  }
  
  #Calculate evaluation metrics
  mae <- mean(abs(predictions$mean - observed))
  rmse <- sqrt(mean((predictions$mean - observed)^2))
  me <- mean(predictions$mean - observed)  #Mean Error (bias)
  
  data.frame(species=species, fold=fold, me=me, mae=mae, rmse=rmse, status=fold_fit$mode$mode.status)
}

fit_mod <- function(task_id, metals){
  start_time <- Sys.time()
  species <- mod_args$species[mod_args$task_id==task_id]
  message(paste0(species, ": Task ID ", task_id, " of ", max(mod_args$task_id)))
  group <- mod_args$group[mod_args$task_id==task_id]
  minor <- mod_args$minor[mod_args$task_id==task_id]
  
  if(metals==TRUE){
    data <- get_data(species, group) #Don't drop any variables
    f_base <- "y ~ PC1 + PC2 + CRI + CAMS + ASR + sewage + barrier_density +
    wT_12M + TIN_12M + PO4_12M + wT_12M:TIN_12M + wT_12M:PO4_12M +
    pH_12M + Cu_d_12M + Zn_d_12M + pH_12M:Cu_d_12M + pH_12M:Zn_d_12M +
    rainfall + f(SITE_ID, model = 'iid', hyper = site_prior) +
    f(julian_week, model = 'rw2') + f(SITE_ID_slope, scale(rainfall), model = 'iid')"
  } else{
    data <- get_data(species, group, dropvar=c("Cu_d", "Zn_d"))
    f_base <- "y ~ PC1 + PC2 + CRI + CAMS + ASR + sewage + barrier_density +
    wT_12M + TIN_12M + PO4_12M + wT_12M:TIN_12M + wT_12M:PO4_12M +
    pH_12M + rainfall + f(SITE_ID, model = 'iid', hyper = site_prior) +
    f(julian_week, model = 'rw2') + f(SITE_ID_slope, scale(rainfall), model = 'iid')"
  }
  data$y[data$y<0] <- 0 #1 observation for a single fish species
  
  data$prob <- 1 #Setting up data for zero inflation
  data$julian_week <- floor((data$julian - 59) / 7) + 1 #Simplify julian to a weekly scale
  data$SITE_ID_slope <- data$SITE_ID #For random slope
  site_prior <- list(prec = list(prior = "pc.prec", param = c(3, 0.05))) #P(SD > 3) = 0.05
  
  if(group=="fish" & minor==FALSE){
    f_final <- paste(f_base, "offset(log(FISHED_AREA))", sep=" + ")
    fit <- inla(as.formula(f_final),
                family = "poisson",
                data = data,
                control.fixed=list(mean = 0, prec = 0.33),
                control.inla=list(cmin=0),
                control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    fit_dis <- dispersion_check(fit)
    if(mean(fit_dis$model <= fit_dis$data)==1){ #Overdispersed
      fit <- inla(as.formula(f_final),
                  family = "nbinomial2",
                  data = data,
                  control.fixed=list(mean = 0, prec = 0.33),
                  control.inla=list(cmin=0),
                  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
      cv_result <- do.call(rbind, lapply(1:5, function(x) do_cv(data, species, group, minor, x, form=f_final, fam="nbinomial2")))
    } else{
      cv_result <- do.call(rbind, lapply(1:5, function(x) do_cv(data, species, group, minor, x, form=f_final, fam="poisson")))
    }
  }
  if(minor==TRUE){
    data$REGION_year <- interaction(data$REGION, data$year) #For nested effect of year within region on zi prob
    f_final <- paste(f_base, "offset(log(FISHED_AREA))",
                     "f(REGION, model='iid', replicate = prob)",
                     "f(REGION_year, model = 'iid', replicate = prob)", sep=" + ")
    fit <- inla(as.formula(f_final),
                family = "zeroinflatednbinomial1",
                data = data,
                control.fixed=list(mean = 0, prec = 0.33),
                control.inla=list(cmin=0),
                control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    cv_result <- do.call(rbind, lapply(1:5, function(x) do_cv(data, species, group, minor, x, form=f_final, fam="zeroinflatednbinomial1")))
  }
  if(group=="invs"){
    f_final <- paste(f_base, "f(year, model = 'iid', replicate = prob)", sep=" + ")
    fit <- inla(as.formula(f_final),
                family = "zeroinflatednbinomial1",
                data = data,
                control.fixed=list(mean = 0, prec = 0.33),
                control.inla=list(cmin=0),
                control.compute = list(dic = TRUE, waic = TRUE, config = TRUE))
    cv_result <- do.call(rbind, lapply(1:5, function(x) do_cv(data, species, group, minor, x, form=f_final, fam="zeroinflatednbinomial1")))
  }
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  list(fit=fit, cv_result=cv_result, runtime=runtime)
}

fit_all <- function(task_id){
  mod_path <- file.path(PATH_PROCESSED, "inla_mods", "mods", paste0("mod_", task_id, ".rds"))
  if(!file.exists(mod_path)){
    mod_nometals <- try(fit_mod(task_id, metals=FALSE))
    mod_metals <- try(fit_mod(task_id, metals=TRUE))
    mods <- list(mod_nometals, mod_metals)
    saveRDS(mods, mod_path, compress = "xz")
  } else{
    NULL
  }
}

##For non-grid job
fit_all(task_id)





