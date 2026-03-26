source(here::here("paths.R"))

##no_metals models
pred_args <- data.frame(
  model_type="no_metals",
  expand.grid(
    CRI=c(FALSE, TRUE), #If TRUE, reduce to zero EXCEPT FOR HMWBs - mask these out
    CAMS=c(FALSE, TRUE), #If TRUE, adjust all to COMPLIANT (or retain as Not Assessed)
    sewage=c(FALSE, TRUE), #If TRUE, adjust sewage to 0 (or FALSE in the df?)
    ASR=c(FALSE, TRUE), #If TRUE, adjust ASR to most favourable category (1?)
    TIN_12M=seq(from=1, to=0.2, by= -0.05), #Proportion to multiply long-term mean by
    PO4_12M=seq(from=1, to=0.2, by= -0.05), #Proportion to multiply long-term mean by
    barrier_density=seq(from=1, to=0, by= -0.1), #Proportion to multiply barrier density by
    metals=FALSE, #To allow rbind
    mining=FALSE, #To allow rbind
    RCP=c(26, 60, 85) #Rainfall (from CHESS-SCAPE rasters) and wT scenario (long-term mean RCP2.6-SSP1; increase at long-term rate RCP6.0-SSP2; doubling of long-term rate RCP8.5-SSP5)
  )
)
pred_args$BAU <- FALSE #TRUE under no change to predictors and RCP85 - Assumed to be equivalent BAU
pred_args$BAU[pred_args$CRI==FALSE & pred_args$CAMS==FALSE & pred_args$sewage==FALSE & pred_args$ASR==FALSE & pred_args$TIN_12M==1 & pred_args$PO4_12M==1 & pred_args$barrier_density==1] <- TRUE

n_metals <- nrow(pred_args)

##metals models
pred_args_metals <- data.frame(
  model_type="metals",
  expand.grid(
    CRI=c(FALSE, TRUE), #If TRUE, reduce to zero EXCEPT FOR HMWBs - mask these out
    CAMS=c(FALSE, TRUE), #If TRUE, adjust all to COMPLIANT (or retain as Not Assessed)
    sewage=c(FALSE, TRUE), #If TRUE, adjust sewage to 0 (FALSE in the df)
    ASR=c(FALSE, TRUE), #If TRUE, adjust ASR to most favourable category (1)
    TIN_12M=seq(from=1, to=0.2, by= -0.05), #Proportion to multiply long-term mean by
    PO4_12M=seq(from=1, to=0.2, by= -0.05), #Proportion to multiply long-term mean by
    barrier_density=seq(from=1, to=0, by= -0.1), #Proportion to multiply barrier density by
    metals=c(FALSE, TRUE), #If TRUE, adjust both copper (3.76) and zinc (7.8) to standard threshold, or long-term mean (whichever is lower)
    mining=c(FALSE, TRUE), #If TRUE and metals==TRUE, only adjust copper and zinc in areas affected by abandoned mine drainage
    RCP=c(26, 60, 85) #Rainfall (from CHESS-SCAPE rasters) and wT scenario (long-term mean RCP2.6-SSP1; increase at long-term rate RCP6.0-SSP2; doubling of long-term rate RCP8.5-SSP5)
  )
)
pred_args_metals <- pred_args_metals[-which(pred_args_metals$metals==FALSE & pred_args_metals$mining==TRUE),] #Redundant combination
pred_args_metals$BAU <- FALSE
pred_args_metals$BAU[pred_args_metals$CRI==FALSE & pred_args_metals$CAMS==FALSE & pred_args_metals$sewage==FALSE & pred_args_metals$ASR==FALSE & pred_args_metals$TIN_12M==1 & pred_args_metals$PO4_12M==1 & pred_args_metals$barrier_density==1 & pred_args_metals$metals==FALSE] <- TRUE

n_no_metals <- nrow(pred_args_metals)

##Prep final pred_args
pred_args <- rbind(pred_args, pred_args_metals)
rm(pred_args_metals)

write.csv(pred_args, file.path(PATH_PROCESSED, "future_pred_args.csv"), row.names=FALSE)

#One task for each species/model_type (up to 266*2 tasks - taken from mod_args)
