#!/usr/bin/Rscript
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
#task_id=266

source(here::here("paths.R"))

####Packages
library(INLA)
library(dplyr)

##Load args
mod_args <- read.csv(file.path(PATH_PROCESSED, "inla_mod_args_filter.csv"), row.names = 1)
mod_args <- mod_args[mod_args$task_id==task_id,] #Continue job only if $metals OR $no_metals == TRUE

if(mod_args$no_metals==TRUE | mod_args$metals==TRUE){
  pred_args <- read.csv(file.path(PATH_PROCESSED, "future_pred_args.csv"))
  pred_args <- list(
    no_metals=pred_args[pred_args$model_type=="no_metals",],
    metals=pred_args[pred_args$model_type=="metals",]
    )
  
  ##Load data for right group
  group <- mod_args$group
  
  if(group=="invs"){
    fut <- read.csv(file.path(PATH_PROCESSED, "invs_future_pred_data.csv")) #HMWB flags and climate
    his <- readRDS(file.path(PATH_PROCESSED, "invs_historical_pred_data.rds")) #list(basedata, scales)
  } else{
    fut <- read.csv(file.path(PATH_PROCESSED, "fish_future_pred_data.csv")) #HMWB flags and climate
    his <- readRDS(file.path(PATH_PROCESSED, "fish_historical_pred_data.rds")) #list(basedata, scales)
  }
  
  ##Species baseline
  baseline <- readRDS(file.path(PATH_PROCESSED, "baselines.rds"))
  baseline <- list(
    no_metals=baseline$no_metals[baseline$no_metals[,"task_id"]==task_id,-1],
    metals=baseline$metals[baseline$metals[,"task_id"]==task_id,-1]
  )
  
  ##Average historical time varying predictors in the period 2018-2020 - this captures the different baseline years in the water targets
  dat <- his$basedata %>%
    mutate(SITE_ID = as.character(SITE_ID)) %>%
    filter(year %in% c(2018, 2019, 2020)) %>%
    group_by(SITE_ID) %>%
    summarise(
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      across(where(is.character), ~ first(.x)),
      across(where(is.logical), ~ first(.x)),
      across(where(is.factor), ~ first(.x)),
      .groups = "drop"
    ) %>%
    arrange(match(SITE_ID, unique(his$basedata$SITE_ID))) %>%
    left_join(fut %>% mutate(SITE_ID = as.character(SITE_ID)), by = "SITE_ID")
  
  ##Save scalings
  scl <- his$scales
  
  ##Remove unnecessary objects
  rm(fut, his)
  
  ##Function to generate newdata
  make_newdata <- function(i, target_year, model_type){
    if (i==1 | i %% 1000 == 0) {
      message(paste(target_year, model_type, i, sep = ": "))
    }
    x <- pred_args[[model_type]][i,]
    d <- dat
    
    #Simulating changes under scenario i
    if (x$CRI) d$CRI <- ifelse(d$HMWB_A, d$CRI, 0) # CRI: reduce to zero except for HMWBs
    if (x$CAMS) d$CAMS <- ifelse(d$CAMS == "Not Assessed", "Not Assessed", "COMPLIANT") # CAMS: adjust all to COMPLIANT (retain Not Assessed)
    d$CAMS <- factor(d$CAMS, levels = levels(dat$CAMS))
    if (x$sewage) d$sewage <- FALSE # sewage: set to FALSE
    if (x$ASR) d$ASR <- 1 # ASR: set to most favourable category
    d$TIN_12M <- d$TIN_12M * x$TIN_12M
    d$PO4_12M <- d$PO4_12M * x$PO4_12M # TIN/PO4: multiply by proportion
    d$barrier_density <- d$barrier_density * x$barrier_density # barrier_density: multiply by proportion
    d$rainfall <- as.data.frame(d)[,paste("rainfall", target_year, x$RCP, sep="_")] # Rainfall: swap in the appropriate year and RCP scenario column
    d <- d[,-which(substr(colnames(d),1,9)=="rainfall_")] #Remove future rainfall columns no longer needed
    years_ahead <- target_year - 2019 #Base year just to be consistent with the 2018-2020 baselines for targets
    d$wT_12M <- switch(as.character(x$RCP),
                       "26" = d$wT_12M, #no change - already at "mean" 
                       "60" = d$wT_12M + 0.03 * years_ahead, #Increase at long-term rate from the literature
                       "85" = d$wT_12M + (2 * 0.03) * years_ahead #Increase at twice the long-term rate
    )
    if (x$model_type == "metals" && x$metals) {
      # Which rows to apply to
      rows <- if (x$mining) d$mining else rep(TRUE, nrow(d))
      d$Cu_d_12M <- ifelse(rows, pmin(d$Cu_d_12M, 3.76, na.rm = FALSE), d$Cu_d_12M)
      d$Zn_d_12M <- ifelse(rows, pmin(d$Zn_d_12M, 7.8,  na.rm = FALSE), d$Zn_d_12M)
    } # Metals: only if model_type == "metals" AND metals == TRUE
    
    #Scaling
    sc <- d
    for (i in seq_len(nrow(scl))) {
      var_name <- scl$variable[i]
      mu       <- scl$mean[i]
      sigma    <- scl$sd[i]
      if (var_name %in% names(sc)) {
        sc[[var_name]] <- 
          (sc[[var_name]] - mu) / sigma
      }
    }
    sc$SITE_ID <- as.character(sc$SITE_ID)
    sc
  }
  
  ##Function to do predictions
  predict_i <- function(form, pred_data, baseline, fixef_samples, ints, slopes) {
    X <- model.matrix(form, pred_data)
    
    # [n_obs x 150] on count scale
    pred_samples <- exp(X %*% fixef_samples + ints + slopes * pred_data$rainfall)
    
    # Identify valid columns (no NaN or Inf anywhere)
    valid_cols <- which(apply(pred_samples, 2, function(x) !any(is.nan(x) | is.infinite(x))))
    
    # Take first 100 valid samples
    pred_samples <- pred_samples[, valid_cols[1:min(100, length(valid_cols))]]
    
    (colSums(pred_samples) / baseline) * 100
  }
  
  run_preds <- function(task_id){
    mod_file <- file.path(PATH_PROCESSED, "inla_mods", "mods", paste0("mod_", task_id, ".rds"))
    
    preds_no_metals <- NULL
    preds_metals <- NULL
    
    if(file.exists(mod_file)) {
      
      if(mod_args$no_metals==TRUE){
        form <- ~ PC1 + PC2 + CRI + CAMS + ASR + sewage +
          barrier_density + wT_12M + TIN_12M + PO4_12M +
          wT_12M:TIN_12M + wT_12M:PO4_12M + pH_12M + rainfall
        pred_data_dummy <- make_newdata(1, 2030, "no_metals") #i, year doesn't matter
        
        post <- readRDS(file.path(PATH_PROCESSED, "inla_mods", "posts", paste0("post_proc_", task_id, "_no_metals.rds")))
        fixef_samples <- post$fixef_samples
        re_int   <- post$intercept_samples[pred_data_dummy$SITE_ID, ]
        re_slope <- post$slope_samples[pred_data_dummy$SITE_ID, ]
        rm(post)
        
        preds_no_metals <- array(dim=c(nrow(pred_args$no_metals), 100, 2))
        dimnames(preds_no_metals)[[3]] <- c("x2030", "x2042")
        for(i in 1:nrow(pred_args$no_metals)){
          preds_no_metals[i,,"x2030"] <- predict_i(form, make_newdata(i, 2030, "no_metals"), baseline$no_metals, fixef_samples, re_int, re_slope)
          preds_no_metals[i,,"x2042"] <- predict_i(form, make_newdata(i, 2042, "no_metals"), baseline$no_metals, fixef_samples, re_int, re_slope)
        }
      }
      
      if(mod_args$metals==TRUE){
        form <- ~ PC1 + PC2 + CRI + CAMS + ASR + sewage +
          barrier_density + wT_12M + TIN_12M + PO4_12M +
          wT_12M:TIN_12M + wT_12M:PO4_12M + pH_12M +
          Cu_d_12M + Zn_d_12M + pH_12M:Cu_d_12M +
          pH_12M:Zn_d_12M + rainfall
        pred_data_dummy <- na.omit(make_newdata(1, 2030, "metals")) #i, year doesn't matter
        
        post <- readRDS(file.path(PATH_PROCESSED, "inla_mods", "posts", paste0("post_proc_", task_id, "_metals.rds")))
        fixef_samples <- post$fixef_samples
        re_int   <- post$intercept_samples[pred_data_dummy$SITE_ID, ]
        re_slope <- post$slope_samples[pred_data_dummy$SITE_ID, ]
        rm(post)
        
        preds_metals <- array(dim=c(nrow(pred_args$metals), 100, 2))
        dimnames(preds_metals)[[3]] <- c("x2030", "x2042")
        for(i in 1:nrow(pred_args$metals)){
          preds_metals[i,,"x2030"] <- predict_i(form, na.omit(make_newdata(i, 2030, "metals")), baseline$metals, fixef_samples, re_int, re_slope)
          preds_metals[i,,"x2042"] <- predict_i(form, na.omit(make_newdata(i, 2042, "metals")), baseline$metals, fixef_samples, re_int, re_slope)
        }
      }
    }
    
    list(no_metals=preds_no_metals, metals=preds_metals)
  }
  
  preds <- run_preds(task_id)
  
  saveRDS(preds, file.path(PATH_PROCESSED, "inla_future_pred", "preds", paste0(task_id, ".rds")))
}


