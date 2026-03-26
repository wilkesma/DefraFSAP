source(here::here("paths.R"))

####Packages
library(INLA)
library(dplyr)
library(parallel)
library(abind)
library(sf)
library(ggplot2)
library(ggspatial)
library(terra)
library(tidyterra)
library(patchwork)
library(ggeffects)

##Load args
mod_args <- read.csv(file.path(PATH_PROCESSED, "inla_mod_args_filter.csv"), row.names = 1)

pred_args <- list(
  no_metals=expand.grid(
    CRI=c(FALSE, TRUE), #If TRUE, reduce to zero EXCEPT FOR HMWBs - mask these out
    CAMS=c(FALSE, TRUE), #If TRUE, adjust all to COMPLIANT (or retain as Not Assessed)
    sewage=c(FALSE, TRUE), #If TRUE, adjust sewage to 0 (or FALSE in the df?)
    ASR=c(FALSE, TRUE), #If TRUE, adjust ASR to most favourable category (1?)
    TIN_12M=c(1,0.5), #Proportion to multiply long-term mean by
    PO4_12M=c(1,0.5), #Proportion to multiply long-term mean by
    barrier_density=c(1,0.5), #Proportion to multiply barrier density by
    metals=FALSE, #To allow rbind
    mining=FALSE, #To allow rbind
    RCP=c(26, 60, 85) #Rainfall (from CHESS-SCAPE rasters) and wT scenario (long-term mean RCP2.6-SSP1; increase at long-term rate RCP6.0-SSP2; doubling of long-term rate RCP8.5-SSP5)
  ),
  metals=expand.grid(
    CRI=c(FALSE, TRUE), #If TRUE, reduce to zero EXCEPT FOR HMWBs - mask these out
    CAMS=c(FALSE, TRUE), #If TRUE, adjust all to COMPLIANT (or retain as Not Assessed)
    sewage=c(FALSE, TRUE), #If TRUE, adjust sewage to 0 (FALSE in the df)
    ASR=c(FALSE, TRUE), #If TRUE, adjust ASR to most favourable category (1)
    TIN_12M=c(1,0.5), #Proportion to multiply long-term mean by
    PO4_12M=c(1,0.5), #Proportion to multiply long-term mean by
    barrier_density=c(1,0.5), #Proportion to multiply barrier density by
    metals=c(FALSE, TRUE), #If TRUE, adjust both copper (3.76) and zinc (7.8) to standard threshold, or long-term mean (whichever is lower)
    mining=c(FALSE, TRUE), #If TRUE and metals==TRUE, only adjust copper and zinc in areas affected by abandoned mine drainage
    RCP=c(26, 60, 85) #Rainfall (from CHESS-SCAPE rasters) and wT scenario (long-term mean RCP2.6-SSP1; increase at long-term rate RCP6.0-SSP2; doubling of long-term rate RCP8.5-SSP5)
  )
)

pred_args$metals <- pred_args$metals[-which(pred_args$metals$metals==FALSE & pred_args$metals$mining==TRUE),] #Redundant combination

pred_args$no_metals$BAU <- FALSE #TRUE under no change to predictors and RCP85 - Assumed to be equivalent BAU
pred_args$no_metals$BAU[pred_args$no_metals$CRI==FALSE & pred_args$no_metals$CAMS==FALSE & pred_args$no_metals$sewage==FALSE & pred_args$no_metals$ASR==FALSE & pred_args$no_metals$TIN_12M==1 & pred_args$no_metals$PO4_12M==1 & pred_args$no_metals$barrier_density==1] <- TRUE

pred_args$metals$BAU <- FALSE #TRUE under no change to predictors and RCP85 - Assumed to be equivalent BAU
pred_args$metals$BAU[pred_args$metals$CRI==FALSE & pred_args$metals$CAMS==FALSE & pred_args$metals$sewage==FALSE & pred_args$metals$ASR==FALSE & pred_args$metals$TIN_12M==1 & pred_args$metals$PO4_12M==1 & pred_args$metals$barrier_density==1 & pred_args$metals$metals==FALSE & pred_args$metals$mining==FALSE] <- TRUE

##Function to generate predictions
get_preds <- function(task_id){
  mod_args_i <- mod_args[mod_args$task_id==task_id,] #Continue job only if $metals OR $no_metals == TRUE
  
  if(mod_args_i$no_metals==TRUE | mod_args_i$metals==TRUE){
    
    ##Load data for right group
    group <- mod_args_i$group
    
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
    
    ##Baseline year (2003) predictors
    basedata <- his$basedata %>%
      filter(year==2003) %>%
      mutate(SITE_ID = as.character(SITE_ID))
    
    for (i in seq_len(nrow(scl))) {
      var_name <- scl$variable[i]
      mu       <- scl$mean[i]
      sigma    <- scl$sd[i]
      if (var_name %in% names(basedata)) {
        basedata[[var_name]] <- 
          (basedata[[var_name]] - mu) / sigma
      }
    }
    
    ##Remove unnecessary objects
    rm(fut, his)
    
    ##Function to generate newdata
    make_newdata <- function(task_id, i, target_year, model_type){
      if (i==1 | i %% 50 == 0) {
        message(paste(task_id, target_year, model_type, i, sep = ": "))
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
      if (model_type == "metals" && x$metals) {
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
      
      if(!is.null(baseline)){
        zero_zero <- baseline < 0.001 & pred_samples < 0.001 #If both effectively zero then the relative abundance index is effectively 100
        near_zero_baseline <- baseline < 0.001 & pred_samples >= 0.001
        pred_samples <- (pred_samples / baseline) * 100
        pred_samples[zero_zero] <- 100
        pred_samples[near_zero_baseline] <- NA
        pred_samples <- data.frame(
          mean=rowMeans(pred_samples, na.rm = TRUE),
          sd=apply(pred_samples, 1, function(x) sd(x, na.rm = TRUE)),
          lwr=apply(pred_samples, 1, function(x) quantile(x, 0.1, na.rm=TRUE)),
          upr=apply(pred_samples, 1, function(x) quantile(x, 0.9, na.rm=TRUE))
        )
      }
      
      pred_samples
    }
    
    run_preds <- function(task_id){
      mod_file <- file.path(PATH_PROCESSED, "inla_mods", "mods", paste0("mod_", task_id, ".rds"))
      
      preds_no_metals_2030 <- NULL
      preds_no_metals_2042 <- NULL
      preds_metals_2030 <- NULL
      preds_metals_2042 <- NULL
      
      if(file.exists(mod_file)) {
        
        if(mod_args_i$no_metals==TRUE){
          form <- ~ PC1 + PC2 + CRI + CAMS + ASR + sewage +
            barrier_density + wT_12M + TIN_12M + PO4_12M +
            wT_12M:TIN_12M + wT_12M:PO4_12M + pH_12M + rainfall
          pred_data_dummy <- make_newdata(task_id, 1, 2030, "no_metals") #i, year doesn't matter
          
          post <- readRDS(file.path(PATH_PROCESSED, "inla_mods", "posts", paste0("post_proc_", task_id, "_no_metals.rds")))
          fixef_samples <- post$fixef_samples
          re_int   <- post$intercept_samples[pred_data_dummy$SITE_ID, ]
          re_slope <- post$slope_samples[pred_data_dummy$SITE_ID, ]
          rm(post)
          
          preds_base <- predict_i(form, basedata, baseline=NULL, fixef_samples, re_int, re_slope)
          
          preds_no_metals_2030 <- array(dim=c(nrow(preds_base), 4, nrow(pred_args$no_metals)))
          dimnames(preds_no_metals_2030) <- list(basedata$SITE_ID,
                                                 c("mean", "sd", "lwr", "upr"),
                                                 NULL)
          preds_no_metals_2042 <- array(dim=c(nrow(preds_base), 4, nrow(pred_args$no_metals)))
          dimnames(preds_no_metals_2042) <- list(basedata$SITE_ID,
                                                 c("mean", "sd", "lwr", "upr"),
                                                 NULL)
          for(i in 1:nrow(pred_args$no_metals)){
            preds_no_metals_2030[,,i] <- as.matrix(predict_i(form, make_newdata(task_id, i, 2030, "no_metals"), preds_base, fixef_samples, re_int, re_slope)) #Works but I need to sort out the dims - add scenarios to get long df with pred_args row, target year, SITE_ID
            preds_no_metals_2042[,,i] <- as.matrix(predict_i(form, make_newdata(task_id, i, 2042, "no_metals"), preds_base, fixef_samples, re_int, re_slope))
          }
        }
        
        if(mod_args_i$metals==TRUE){
          form <- ~ PC1 + PC2 + CRI + CAMS + ASR + sewage +
            barrier_density + wT_12M + TIN_12M + PO4_12M +
            wT_12M:TIN_12M + wT_12M:PO4_12M + pH_12M +
            Cu_d_12M + Zn_d_12M + pH_12M:Cu_d_12M +
            pH_12M:Zn_d_12M + rainfall
          pred_data_dummy <- na.omit(make_newdata(task_id, 1, 2030, "metals")) #i, year doesn't matter
          
          post <- readRDS(file.path(PATH_PROCESSED, "inla_mods", "posts", paste0("post_proc_", task_id, "_metals.rds")))
          fixef_samples <- post$fixef_samples
          re_int   <- post$intercept_samples[pred_data_dummy$SITE_ID, ]
          re_slope <- post$slope_samples[pred_data_dummy$SITE_ID, ]
          rm(post)
          
          basedata_metals <- na.omit(basedata)
          
          preds_base <- predict_i(form, basedata_metals, baseline=NULL, fixef_samples, re_int, re_slope)
          
          preds_metals_2030 <- array(dim=c(nrow(preds_base), 4, nrow(pred_args$metals)))
          dimnames(preds_metals_2030) <- list(basedata_metals$SITE_ID,
                                              c("mean", "sd", "lwr", "upr"),
                                              NULL)
          preds_metals_2042 <- array(dim=c(nrow(preds_base), 4, nrow(pred_args$metals)))
          dimnames(preds_metals_2042) <- list(basedata_metals$SITE_ID,
                                              c("mean", "sd", "lwr", "upr"),
                                              NULL)
          for(i in 1:nrow(pred_args$metals)){
            preds_metals_2030[,,i] <- as.matrix(predict_i(form, na.omit(make_newdata(task_id, i, 2030, "metals")), preds_base, fixef_samples, re_int, re_slope)) #Works but I need to sort out the dims - add scenarios to get long df with pred_args row, target year, SITE_ID
            preds_metals_2042[,,i] <- as.matrix(predict_i(form, na.omit(make_newdata(task_id, i, 2042, "metals")), preds_base, fixef_samples, re_int, re_slope))
          }
        }
      }
      
      list(no_metals=list(x2030=preds_no_metals_2030, x2042=preds_no_metals_2042),
           metals=list(x2030=preds_metals_2030, x2042=preds_metals_2042))
    }
    
    preds <- run_preds(task_id)
    
    saveRDS(preds, file.path(PATH_PROCESSED, "inla_spatial_pred", "preds", paste0(task_id, ".rds")))
  }
}

mclapply(mod_args$task_id, get_preds, mc.cores=150)

##Function to predict species relative abundance index for 2022
get_preds_2022 <- function(task_id){
  mod_args_i <- mod_args[mod_args$task_id==task_id,] #Continue job only if $metals OR $no_metals == TRUE
  
  if(mod_args_i$no_metals==TRUE | mod_args_i$metals==TRUE){
    
    ##Load data for right group
    group <- mod_args_i$group
    
    if(group=="invs"){
      his <- readRDS(file.path(PATH_PROCESSED, "invs_historical_pred_data.rds")) #list(basedata, scales)
    } else{
      his <- readRDS(file.path(PATH_PROCESSED, "fish_historical_pred_data.rds")) #list(basedata, scales)
    }
    
    ##Species baseline
    baseline <- readRDS(file.path(PATH_PROCESSED, "baselines.rds"))
    baseline <- list(
      no_metals=baseline$no_metals[baseline$no_metals[,"task_id"]==task_id,-1],
      metals=baseline$metals[baseline$metals[,"task_id"]==task_id,-1]
    )
    
    ##Baseline year (2003) predictors
    basedata <- his$basedata %>%
      filter(year==2003) %>%
      mutate(SITE_ID = as.character(SITE_ID))
    
    ##Average historical time varying predictors in the period 2018-2020 - this captures the different baseline years in the water targets
    dat <- his$basedata %>%
      mutate(SITE_ID = as.character(SITE_ID)) %>%
      filter(year %in% c(2022)) %>%
      group_by(SITE_ID) %>%
      arrange(match(SITE_ID, unique(his$basedata$SITE_ID)))
    
    ##Scaling
    scl <- his$scales
    
    for (i in seq_len(nrow(scl))) {
      var_name <- scl$variable[i]
      mu       <- scl$mean[i]
      sigma    <- scl$sd[i]
      if (var_name %in% names(basedata)) {
        basedata[[var_name]] <- 
          (basedata[[var_name]] - mu) / sigma
      }
      if (var_name %in% names(dat)) {
        dat[[var_name]] <- 
          (dat[[var_name]] - mu) / sigma
      }
    }
    
    ##Remove unnecessary objects
    rm(his)
    
    ##Function to do predictions
    predict_i <- function(form, pred_data, baseline, fixef_samples, ints, slopes) {
      X <- model.matrix(form, pred_data)
      
      pred_samples <- exp(X %*% fixef_samples + ints + slopes * pred_data$rainfall)
      
      valid_cols <- which(apply(pred_samples, 2, function(x) all(is.finite(x))))
      pred_samples <- pred_samples[, valid_cols[1:min(100, length(valid_cols))]]
      
      if(!is.null(baseline)){
        zero_zero <- baseline < 0.001 & pred_samples < 0.001 #If both effectively zero then the relative abundance index is effectively 100
        near_zero_baseline <- baseline < 0.001 & pred_samples >= 0.001
        pred_samples <- (pred_samples / baseline) * 100
        pred_samples[zero_zero] <- 100
        pred_samples[near_zero_baseline] <- NA
        pred_samples <- data.frame(
          mean=rowMeans(pred_samples, na.rm = TRUE),
          sd=apply(pred_samples, 1, function(x) sd(x, na.rm = TRUE)),
          lwr=apply(pred_samples, 1, function(x) quantile(x, 0.1, na.rm=TRUE)),
          upr=apply(pred_samples, 1, function(x) quantile(x, 0.9, na.rm=TRUE))
        )
      }
      
      pred_samples
    }
    
    run_preds <- function(task_id){
      mod_file <- file.path(PATH_PROCESSED, "inla_mods", "mods", paste0("mod_", task_id, ".rds"))
      
      preds_no_metals <- NULL
      preds_metals <- NULL
      
      if(file.exists(mod_file)) {
        
        if(mod_args_i$no_metals==TRUE){
          form <- ~ PC1 + PC2 + CRI + CAMS + ASR + sewage +
            barrier_density + wT_12M + TIN_12M + PO4_12M +
            wT_12M:TIN_12M + wT_12M:PO4_12M + pH_12M + rainfall
          
          post <- readRDS(file.path(PATH_PROCESSED, "inla_mods", "posts", paste0("post_proc_", task_id, "_no_metals.rds")))
          fixef_samples <- post$fixef_samples
          re_int   <- post$intercept_samples[dat$SITE_ID, ]
          re_slope <- post$slope_samples[dat$SITE_ID, ]
          rm(post)
          
          preds_base <- predict_i(form, basedata, baseline=NULL, fixef_samples, re_int, re_slope)
          preds_no_metals <- as.matrix(predict_i(form, dat, preds_base, fixef_samples, re_int, re_slope))
          row.names(preds_no_metals) <- basedata$SITE_ID
        }
        
        if(mod_args_i$metals==TRUE){
          form <- ~ PC1 + PC2 + CRI + CAMS + ASR + sewage +
            barrier_density + wT_12M + TIN_12M + PO4_12M +
            wT_12M:TIN_12M + wT_12M:PO4_12M + pH_12M +
            Cu_d_12M + Zn_d_12M + pH_12M:Cu_d_12M +
            pH_12M:Zn_d_12M + rainfall
          
          post <- readRDS(file.path(PATH_PROCESSED, "inla_mods", "posts", paste0("post_proc_", task_id, "_metals.rds")))
          fixef_samples <- post$fixef_samples
          re_int   <- post$intercept_samples[na.omit(dat)$SITE_ID, ]
          re_slope <- post$slope_samples[na.omit(dat)$SITE_ID, ]
          rm(post)
          
          preds_base <- predict_i(form, na.omit(basedata), baseline=NULL, fixef_samples, re_int, re_slope)
          preds_metals <- as.matrix(predict_i(form, na.omit(dat), preds_base, fixef_samples, re_int, re_slope))
          row.names(preds_metals) <- na.omit(basedata)$SITE_ID
        }
      }
      
      list(no_metals=preds_no_metals, metals=preds_metals)
    }
    
    run_preds(task_id)
  }
}

preds_2022 <- mclapply(mod_args$task_id, get_preds_2022, mc.cores=120)

##Get site-level 2022 geometric means
get_geomean_2022 <- function(model_type, group) {
  # Get all non-null matrices for this model_type
  mats <- lapply(preds_2022[mod_args$group==group], function(x) x[[model_type]])
  mats <- abind(mats[!sapply(mats, is.null)], along=3)[,c("mean", "lwr", "upr"),]
  apply(mats, c(1,2), function(x) exp(mean(log(x), na.rm=TRUE)))
}

preds_2022 <- list(
  invs=list(
    no_metals = get_geomean_2022("no_metals", "invs"),
    metals    = get_geomean_2022("metals", "invs")
  ),
  fish=list(
    no_metals = get_geomean_2022("no_metals", "fish"),
    metals    = get_geomean_2022("metals", "fish")
  )
)

##Site-level future geometric means
get_geomean_preds <- function(group, n_cores) {
  task_ids_no_metals <- mod_args$task_id[mod_args$group == group & mod_args$no_metals == TRUE]
  task_ids_metals    <- mod_args$task_id[mod_args$group == group & mod_args$metals == TRUE]
  
  read_preds <- function(model_type, task_ids, n_cores) {
    mats <- mclapply(task_ids, function(id) {
      message(id)
      preds <- readRDS(paste0("./preds/preds_", id, ".rds"))
      list(
        x2030 = preds[[model_type]]$x2030,
        x2042 = preds[[model_type]]$x2042
      )
    }, mc.cores=n_cores)
    list(
      x2030 = exp(rowMeans(log(abind(lapply(mats, `[[`, "x2030"), along=4)), dims=3, na.rm=TRUE)),
      x2042 = exp(rowMeans(log(abind(lapply(mats, `[[`, "x2042"), along=4)), dims=3, na.rm=TRUE))
    )
  }
  
  list(
    no_metals = read_preds("no_metals", task_ids_no_metals, n_cores),
    metals    = read_preds("metals",    task_ids_metals, n_cores)
  )
}

preds <- list(
  invs = get_geomean_preds("invs", 120),
  fish = get_geomean_preds("fish", 120)
)

##Check SITE_IDs with NAs in preds
na_sites <- list()
for(group in c("invs", "fish")) {
  for(model_type in c("no_metals", "metals")) {
    for(target_year in c("x2030", "x2042")) {
      arr <- preds[[group]][[model_type]][[target_year]]
      if(is.null(arr)) next
      na_sites[[paste(group, model_type, target_year, sep="_")]] <- 
        dimnames(arr)[[1]][apply(arr, 1, function(x) any(is.na(x)))]
    }
  }
}
na_sites #Only affects invs

na_sites <- unique(do.call(c, na_sites)) #Remove these from invs spatial preds for all years/model_types

##Calculate site-level probabilities
process_preds <- function(group, model_type, target_year) {
  preds_i      <- preds[[group]][[model_type]][[paste0("x", target_year)]]
  preds_2022_i <- preds_2022[[group]][[model_type]]
  
  which_bau      <- which(pred_args[[model_type]]$BAU == TRUE)
  pred_args_bau  <- pred_args[[model_type]][which_bau, ]
  preds_bau      <- preds_i[,, which_bau, drop=FALSE]
  pred_args_scen <- pred_args[[model_type]][-which_bau, ]
  preds_scen     <- preds_i[,, -which_bau, drop=FALSE]
  
  site_ids <- rownames(preds_2022_i)
  n_sites  <- length(site_ids)
  
  approx_draws <- function(x, n=100){
    sd_x <- (x["upr"] - x["lwr"]) / (2 * qnorm(0.9))
    rnorm(n, mean=x["mean"], sd=sd_x)
  }
  
  draws_2022 <- t(apply(preds_2022_i, 1, approx_draws))
  
  bau_draws <- lapply(1:dim(preds_bau)[3], function(j)
    t(apply(preds_bau[,,j], 1, approx_draws)))
  
  # p_cf per site per RCP
  p_cf_lookup <- do.call(rbind, lapply(1:dim(preds_bau)[3], function(j) {
    data.frame(
      SITE_ID = site_ids,
      RCP     = pred_args_bau$RCP[j],
      p_cf    = rowMeans(bau_draws[[j]] > draws_2022)
    )
  }))
  
  # One row per site per scenario
  do.call(rbind, lapply(1:dim(preds_scen)[3], function(j) {
    draws_scen  <- t(apply(preds_scen[,,j], 1, approx_draws))
    bau_j       <- which(pred_args_bau$RCP == pred_args_scen$RCP[j])
    draws_bau   <- bau_draws[[bau_j]]
    rcp_j       <- pred_args_scen$RCP[j]
    pred_arg_idx <- which(!pred_args[[model_type]]$BAU)[j]  # index in original pred_args
    
    df <- data.frame(
      group        = group,
      model_type   = model_type,
      year         = target_year,
      SITE_ID      = site_ids,
      RCP          = rcp_j,
      pred_arg_idx = pred_arg_idx,
      mean_2022    = preds_2022_i[, "mean"],
      mean_future  = preds_scen[, "mean", j],
      p_scen       = rowMeans(draws_scen > draws_2022),
      p_supports   = rowMeans(draws_scen > draws_bau),
      p_decisive   = rowMeans(draws_bau <= draws_2022 & draws_scen > draws_2022)
    )
    left_join(df, p_cf_lookup, by=c("SITE_ID", "RCP"))
  }))
}

site_coords_invs <- readRDS(file.path(PATH_PROCESSED, "invs_env_data_complete_2003_2023.rds")) %>%
  select(SITE_ID, easting, northing) %>%
  distinct(SITE_ID, .keep_all=TRUE)

site_coords_fish <- readRDS(file.path(PATH_PROCESSED, "fish_env_data_complete_2003_2023.rds")) %>%
  select(SITE_ID, easting, northing) %>%
  distinct(SITE_ID, .keep_all=TRUE)

probs <- list()
for(model_type in c("no_metals", "metals")) {
  probs[[model_type]] <- list()
  for(group in c("invs", "fish")) {
    probs[[model_type]][[group]] <- list()
    coords <- if(group=="invs") site_coords_invs else site_coords_fish
    for(target_year in c(2030, 2042)) {
      key <- as.character(target_year)
      message(paste(model_type, group, target_year))
      df <- process_preds(group, model_type, target_year)
      probs[[model_type]][[group]][[key]] <- left_join(df, coords, by="SITE_ID")
    }
  }
}

##Pull out single interventions for visualising in report
get_single_intervention_idx <- function(model_type) {
  pa <- pred_args[[model_type]]
  pa_idx <- which(!pa$BAU)  # indices in original pred_args
  pa <- pa[pa_idx, ]
  
  interventions <- list(
    CRI             = pa$CRI == TRUE,
    CAMS            = pa$CAMS == TRUE,
    sewage          = pa$sewage == TRUE,
    ASR             = pa$ASR == TRUE,
    TIN_12M         = pa$TIN_12M == 0.5,
    PO4_12M         = pa$PO4_12M == 0.5,
    barrier_density = pa$barrier_density == 0.5
  )
  
  if(model_type == "metals") {
    interventions$metals <- pa$metals == TRUE & pa$mining == FALSE
    interventions$mining <- pa$metals == TRUE & pa$mining == TRUE
  }
  
  lapply(names(interventions), function(name) {
    this_active  <- interventions[[name]]
    other_active <- Reduce(`|`, interventions[names(interventions) != name])
    idx <- which(this_active & !other_active)
    if(model_type == "metals" && !name %in% c("metals", "mining")) {
      idx <- idx[pa$metals[idx] == FALSE & pa$mining[idx] == FALSE]
    }
    pa_idx[idx]  # return indices in original pred_args
  }) |> setNames(names(interventions))
}

single_intervention_idx <- list(
  no_metals = get_single_intervention_idx("no_metals"),
  metals    = get_single_intervention_idx("metals")
)

# Extract single intervention rows from probs
si_probs <- lapply(c("no_metals", "metals"), function(model_type) {
  idxs <- single_intervention_idx[[model_type]]
  lapply(c("invs", "fish"), function(group) {
    lapply(c("2030", "2042"), function(year) {
      df <- probs[[model_type]][[group]][[year]]
      do.call(rbind, lapply(names(idxs), function(intervention) {
        rows <- df$pred_arg_idx %in% idxs[[intervention]]
        df_sub <- df[rows, ]
        df_sub$intervention <- intervention
        df_sub
      }))
    }) |> setNames(c("2030", "2042"))
  }) |> setNames(c("invs", "fish"))
}) |> setNames(c("no_metals", "metals"))

si_probs$no_metals <- rbind(
  do.call(rbind, si_probs$no_metals$invs),
  do.call(rbind, si_probs$no_metals$fish))

si_probs$metals <- rbind(
  do.call(rbind, si_probs$metals$invs),
  do.call(rbind, si_probs$metals$fish))
  
##Visualise
#Maps
england <- st_transform(st_read(PATH_COUNTRIES), 27700)
england <- england[england$ctry17nm=="England",]

rcp_labels <- c("26"="RCP 2.6", "60"="RCP 6.0", "85"="RCP 8.5")

make_raster_stack_rcp <- function(model_type, group, year, variable, res=25000) {
  df <- si_probs[[model_type]][si_probs[[model_type]]$group==group & si_probs[[model_type]]$year==as.numeric(year), ]
  sf_df <- st_as_sf(df, coords=c("easting", "northing"), crs=27700)
  template <- rast(ext(vect(sf_df)), resolution=res, crs="EPSG:27700")
  
  rast_list <- lapply(unique(df$intervention), function(intv) {
    lapply(unique(df$RCP), function(rcp) {
      sub <- sf_df[sf_df$intervention==intv & sf_df$RCP==rcp, ]
      r <- rasterize(vect(sub), template, field=variable, fun="mean")
      names(r) <- paste(intv, rcp, sep="__")
      r
    })
  })
  
  rast(unlist(rast_list, recursive=FALSE))
}

plot_grid <- function(model_type, group, year, variable, subscript, interventions=NULL) {
  r <- make_raster_stack_rcp(model_type, group, year, variable)
  
  df_scen <- as.data.frame(r, xy=TRUE) %>%
    tidyr::pivot_longer(-c(x,y), names_to="lyr", values_to="value") %>%
    tidyr::separate(lyr, into=c("intervention","RCP"), sep="__") %>%
    mutate(RCP=rcp_labels[RCP]) %>%
    filter(!is.na(value))
  
  if(!is.null(interventions)) df_scen <- df_scen[df_scen$intervention %in% interventions, ]
  
  # BAU/counterfactual row - p_cf is same for all interventions so just take first
  r_bau <- make_raster_stack_rcp(model_type, group, year, "p_cf")
  df_bau <- as.data.frame(r_bau, xy=TRUE) %>%
    tidyr::pivot_longer(-c(x,y), names_to="lyr", values_to="value") %>%
    tidyr::separate(lyr, into=c("intervention","RCP"), sep="__") %>%
    mutate(RCP=rcp_labels[RCP]) %>%
    filter(!is.na(value)) %>%
    distinct(x, y, RCP, .keep_all=TRUE) %>%  # one value per cell per RCP
    mutate(intervention="BAU")
  
  df <- rbind(df_scen, df_bau) %>%
    mutate(intervention=factor(intervention, 
                               levels=c("BAU", unique(df_scen$intervention))))
  
  ggplot(df, aes(x=x, y=y, fill=value)) +
    geom_raster() +
    geom_sf(data=england, inherit.aes=FALSE, fill=NA, colour="black", linewidth=0.2) +
    scale_fill_gradient2(low='#ef8a62', mid='#f7f7f7', high='#67a9cf',
                         midpoint=0.5, na.value="darkgrey", 
                         name=bquote(P[.(subscript)]),
                         limits=c(0,1),
                         guide=guide_colorbar(ticks.colour="black")) +
    facet_grid(intervention~RCP) +
    coord_sf(crs=27700, expand=FALSE) +
    theme_void() +
    labs(title=year) +
    theme(
      plot.margin=margin(10, 10, 10, 10, "mm"),
      strip.text.y=element_text(angle=0, hjust=0, margin=margin(l=3, unit="mm")),
      strip.text.x=element_text(margin=margin(b=3, t=2, unit="mm"))
    )
}

ggsave("./maps/invs_no_metals_2030.png", plot_grid("no_metals", "invs", 2030, "p_supports", "supports"),
       width=7, height=8, dpi=300, units="in")
ggsave("./maps/invs_no_metals_2042.png", plot_grid("no_metals", "invs", 2042, "p_supports", "supports"),
       width=7, height=8, dpi=300, units="in")
ggsave("./maps/fish_no_metals_2030.png", plot_grid("no_metals", "fish", 2030, "p_supports", "supports"),
       width=7, height=8, dpi=300, units="in")
ggsave("./maps/fish_no_metals_2042.png", plot_grid("no_metals", "fish", 2042, "p_supports", "supports"),
       width=7, height=8, dpi=300, units="in")


ggsave("./maps/invs_metals_2030.png", plot_grid("metals", "invs", 2030, "p_supports", "supports", c("metals", "mining")),
       width=7, height=3.5, dpi=300, units="in")
ggsave("./maps/invs_metals_2042.png", plot_grid("metals", "invs", 2042, "p_supports", "supports", c("metals", "mining")),
       width=7, height=3.5, dpi=300, units="in")
ggsave("./maps/fish_metals_2030.png", plot_grid("metals", "fish", 2030, "p_supports", "supports", c("metals", "mining")),
       width=7, height=3.5, dpi=300, units="in")
ggsave("./maps/fish_metals_2042.png", plot_grid("metals", "fish", 2042, "p_supports", "supports", c("metals", "mining")),
       width=7, height=3.5, dpi=300, units="in")

#Density plots
plot_density_grid <- function(model_type, group, year, interventions=NULL) {
  df <- si_probs[[model_type]][si_probs[[model_type]]$group==group & 
                                 si_probs[[model_type]]$year==as.numeric(year), ]
  
  if(!is.null(interventions)) df <- df[df$intervention %in% interventions, ]
  
  df_long <- df %>%
    select(SITE_ID, RCP, intervention, p_cf, p_scen, p_supports, p_decisive) %>%
    tidyr::pivot_longer(cols=c(p_cf, p_scen, p_supports, p_decisive),
                        names_to="p_type", values_to="value") %>%
    mutate(
      RCP=rcp_labels[as.character(RCP)],
      p_type=factor(p_type, levels=c("p_cf","p_scen","p_supports","p_decisive"),
                    labels=c(expression(P[cf]), expression(P[scen]),
                             expression(P[supports]), expression(P[decisive])))
    )
  
  ggplot(df_long, aes(x=value, fill=p_type, colour=p_type)) +
    geom_density(alpha=0.5, linewidth=0.3, bounds=c(0,1)) +
    scale_fill_brewer(palette="Set1", name=NULL,
                      labels=c(expression(P[cf]), expression(P[scen]),
                               expression(P[supports]), expression(P[decisive]))) +
    scale_colour_brewer(palette="Set1", name=NULL,
                        labels=c(expression(P[cf]), expression(P[scen]),
                                 expression(P[supports]), expression(P[decisive]))) +
    scale_x_continuous(limits=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1")) +
    geom_vline(xintercept=0.5, linetype=2, linewidth=0.33) +
    facet_grid(intervention~RCP) +
    labs(x="Probability", y="Density", title=paste(group, model_type, year)) +
    theme_bw() +
    theme(
      plot.margin=margin(5, 5, 5, 5, "mm"),
      strip.text.y=element_text(angle=0, hjust=0),
      strip.text.x=element_text(margin=margin(b=2, t=2, unit="mm")),
      legend.position="top"
    )
}

pdf("P densities.pdf", width=7, height=8)
plot_density_grid("no_metals", "invs", 2030)
plot_density_grid("no_metals", "invs", 2042)
plot_density_grid("no_metals", "fish", 2030)
plot_density_grid("no_metals", "fish", 2042)
plot_density_grid("metals", "invs", 2030)
plot_density_grid("metals", "invs", 2042)
plot_density_grid("metals", "fish", 2030)
plot_density_grid("metals", "fish", 2042)
dev.off()

pdf("P densities metals.pdf", width=7, height=4)
plot_density_grid("metals", "invs", 2030, c("metals", "mining"))
plot_density_grid("metals", "invs", 2042, c("metals", "mining"))
plot_density_grid("metals", "fish", 2030, c("metals", "mining"))
plot_density_grid("metals", "fish", 2042, c("metals", "mining"))
dev.off()

##Assess P_supports relationship to urban, rural, PC1 and PC2 values
lc_2020 <- app(rast(PATH_LCM)[[20:21]], sum) #Urban + suburban
get_mod_data <- function(group, model_type){
  lc <- data.frame(unique(si_probs[[model_type]][si_probs[[model_type]]$group==group,c("SITE_ID", "easting", "northing")]), year=2020)
  lc <- st_as_sf(lc, coords=c("easting", "northing"), crs=27700)
  lc <- extract(lc_2020, lc, ID=FALSE)
  
  lc <- data.frame(
    SITE_ID=unique(si_probs[[model_type]][si_probs[[model_type]]$group==group,c("SITE_ID", "easting", "northing")]$SITE_ID),
    urban=lc[,1]
  )
  
  output <- left_join(unique(readRDS(file.path(PATH_PROCESSED, paste0(group, "_historical_pred_data.rds")))$basedata[,c("SITE_ID", "PC1", "PC2")]),
                      si_probs[[model_type]][si_probs[[model_type]]$group=="invs",])
  left_join(output, lc)
}

mod_data <- list()
mod_data$no_metals <- list()
mod_data$metals <- list()
mod_data$no_metals$invs <- get_mod_data("invs", "no_metals")
mod_data$no_metals$fish <- get_mod_data("fish", "no_metals")
mod_data$metals$invs <- na.omit(get_mod_data("invs", "metals"))
mod_data$metals$fish <- na.omit(get_mod_data("fish", "metals"))

fit_mod <- function(group, model_type, year){
  mod <- glm(p_supports ~ (PC1 + PC2 + urban) * intervention * RCP,
        data=mod_data[[model_type]][[group]][mod_data[[model_type]][[group]]$year==year,], family="binomial")
  pred_PC1   <- ggpredict(mod, terms = c("PC1 [all]", "intervention", "RCP"))
  pred_PC2   <- ggpredict(mod, terms = c("PC2 [all]", "intervention", "RCP"))
  pred_urban <- ggpredict(mod, terms = c("urban [all]", "intervention", "RCP"))
  pred_all <- bind_rows(
    as.data.frame(pred_PC1)   |> mutate(predictor = "PC1"),
    as.data.frame(pred_PC2)   |> mutate(predictor = "PC2"),
    as.data.frame(pred_urban) |> mutate(predictor = "Urban")
  )
  if(model_type=="metals"){
    pred_all <- pred_all[pred_all$group %in% c("metals", "mining"),]
  }
  intervention_order <- mod_data[[model_type]][[group]]$intervention |> 
    unique() |> 
    as.character()
  pred_all <- pred_all |>
    mutate(
      group = recode(as.character(group),
                     "TIN_12M" = "TIN",
                     "PO4_12M" = "PO4"),
      group = factor(group, levels = recode(intervention_order,
                                            "TIN_12M" = "TIN",
                                            "PO4_12M" = "PO4")),
      predictor = recode(predictor,
                         "PC1" = "PC1 (altitude)",
                         "PC2" = "PC2 (river size)"),
      RCP_label = recode(as.character(facet),
                         "26" = "RCP2.6",
                         "60" = "RCP6.0",
                         "85" = "RCP8.5")
    )
  p <- ggplot(pred_all, aes(x = x, y = predicted, colour = RCP_label, fill = RCP_label)) +
    geom_hline(yintercept=0.5, linewidth=0.33, linetype=2) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, colour = NA) +
    geom_line() +
    facet_grid(group ~ predictor, scales = "free_x") +
    labs(x = "Predictor value",  y = expression(P[supports]), title=paste(group, model_type, year, sep=": ")) +
    theme_bw() +
    theme(legend.title=element_blank())
  print(p)
}

pdf("Logistic regression preds no_metals.pdf", width=7, height=8)
lapply(c("invs", "fish"), function(group) lapply(c(2030, 2042), function(year) fit_mod(group, "no_metals", year)))
dev.off()

pdf("Logistic regression preds metals.pdf", width=7, height=3)
lapply(c("invs", "fish"), function(group) lapply(c(2030, 2042), function(year) fit_mod(group, "metals", year)))
dev.off()

##Save critical object
saveRDS(probs, file.path(PATH_PROCESSED, "spatial_probs.rds")) #Will need a system in the shiny app to index this
