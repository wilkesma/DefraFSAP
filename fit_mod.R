#!/usr/bin/Rscript
task.id <- as.integer(Sys.getenv("SGE_TASK_ID"))

source(here::here("paths.R"))
dir.create(file.path(PATH_PROCESSED, "mods"), recursive = TRUE, showWarnings = FALSE)

####Packages
library(dplyr)
library(ggplot2)
library(mgcv)
library(tictoc)
library(lubridate)
library(tidyr)
library(gratia)
library(ddspWQ)
library(gridExtra)
library(stringr)

mod.args <- read.csv(file.path(PATH_PROCESSED, "modargs.csv"), row.names=1)

wq_matched <- read.csv(file.path(PATH_PROCESSED, "wq_matched.csv"), row.names=1)
wq_matched$area <- substr(wq_matched$sample.samplingPoint.notation,1,2)

post_process_predictions <- function(pred.df) {
  cols = c(Concentration = '.fitted', lower.ci = '.lower_ci', upper.ci = '.upper_ci')
  
  pp <- 
    pred.df |> 
    mutate(across(all_of(cols), ~if_else(transform == 'log', 10^.x, .x)), 
           Date = decimal_date(ymd(paste(year, month, '15', sep = '-'))), 
           season = factor(case_match(month, 3:5 ~ 'Spring', 
                                      6:8 ~ 'Summer', 
                                      9:11 ~ 'Autumn', 
                                      c(12,1,2) ~ 'Winter'), 
                           levels = c('Spring', 'Summer', 'Autumn', 'Winter')))
  
  return(pp)
}

fit.mod <- function(task.id){
  det <- mod.args$det[task.id]
  area <- mod.args$area[task.id]
  start <- Sys.time()
  det.data <- wq_matched[wq_matched$determinand.name==det & wq_matched$area==area,]
  
  sknots <- list(month=c(0.5, 12.5))
  nearly.zero <- 0.000001
  log.zero <- log10(nearly.zero)
  
  samples <- 
    det.data[,] |>
    select(notation = sample.samplingPoint.notation, 
           DATE = date, 
           DETE_SHORT_DESC = determinand.name,
           MEAS_SIGN = resultQualifier.notation, 
           MEAS_RESULT = result,
           sub_id = sub_id) |>
    mutate(DATE = date(DATE), 
           year = year(DATE), 
           month = month(DATE),
           notation = as.factor(notation),
           DATE_DEC = round(decimal_date(DATE),2),
           MEAS_SIGN = replace_na(MEAS_SIGN, '='),
           logMEAS_RESULT = log10(MEAS_RESULT + nearly.zero),
           sqrtMEAS_RESULT = sqrt(MEAS_RESULT),
           meas.cen = if_else(MEAS_SIGN == '<', ifelse(det=="O2_d", 0, log.zero), ifelse(det=="O2_d", sqrtMEAS_RESULT, logMEAS_RESULT)))
  samples$MEAS_SIGN[samples$MEAS_SIGN==""] <- "="
  print(samples |>
          group_by(DETE_SHORT_DESC, MEAS_SIGN) |>
          tally() |>
          pivot_wider(names_from = MEAS_SIGN, values_from = n))
  
  if(det %in% c("pH", "wT")){
    y <- samples$MEAS_RESULT
  } else{
    y <- cbind(samples$meas.cen, samples$logMEAS_RESULT) # for lod values, meas.cen is 0
  }
  samples$sub_id <- as.factor(samples$sub_id)
  
  gc(full = TRUE, verbose = TRUE)
  
  if(det %in% c("pH", "wT")){
    mod <- bam(y ~ s(month, bs = 'cc') +
                 s(sub_id, DATE_DEC, bs = 'fs', k=5) +
                 s(notation, bs = "re", k=5), #s(notation, bs = "re", by = sub_id, k=5, xt=list(sparse=TRUE)) would allow different levels of shrinkage of sites towards their subsegment means but memory demand was impractical
               method = 'fREML', discrete=TRUE, chunk.size=1, use.chol=TRUE, select=TRUE, #chunk.size will reset to 4*p, where p is the number of coefficients (approx 13,000)
               family = gaussian(), knots = sknots, data = samples, nthreads=1)
  } else{
    mod <- bam(y ~ s(month, bs = 'cc') +
                 s(sub_id, DATE_DEC, bs = 'fs', k=5) +
                 s(notation, bs = "re", k=5), #s(notation, bs = "re", by = sub_id, k=5, xt=list(sparse=TRUE)) would allow different levels of shrinkage of sites towards their subsegment means but memory demand was impractical
               method = 'fREML', discrete=TRUE, chunk.size=1, use.chol=TRUE, select=TRUE, #chunk.size will reset to 4*p, where p is the number of coefficients (approx 13,000)
               family = cnorm(), knots=sknots, data = samples, nthreads=1)
  }
  
  gc(full = TRUE, verbose = TRUE)
  
  date.sequence <- 
    expand_grid(year = min(samples$year):max(samples$year), month = 1:12) |> 
    mutate(DATE_DEC = decimal_date(ymd(paste(year, month, '15', sep = '-')))) |> as.data.frame() #NOTE: Predictions only for year range with data
  ds <- expand_grid(date.sequence, notation = unique(samples$notation))
  ds <- left_join(ds, unique(samples[,c("notation", "sub_id")]))
  
  if(det %in% c("pH", "wT")){
    preds <- 
      fitted_values(mod, data = ds, scale = 'response') |>
      mutate(transform = 'none') |>
      post_process_predictions()
  } else{
    preds <- 
      fitted_values(mod, data = ds, scale = 'response') |>
      mutate(transform = 'log') |>
      post_process_predictions()
  }
  preds$date <- as.Date(paste(preds$year, ifelse(nchar(as.character(preds$month))==1, paste0("0", preds$month), as.character(preds$month)), "15", sep="-"))
  
  preds.fixef <- fitted_values(mod, data = ds, scale = 'response')
  preds.fixef <- aggregate(cbind(.fitted, .se)~year+month+DATE_DEC+sub_id, preds.fixef, FUN=mean)
  preds.fixef$.lower_ci <- preds.fixef$.fitted - (1.96*preds.fixef$.se)
  preds.fixef$.upper_ci <- preds.fixef$.fitted + (1.96*preds.fixef$.se)
  preds.fixef$transform <- ifelse(det %in% c("pH", "wT"), "none", "log")
  preds.fixef <- post_process_predictions(preds.fixef)
  preds.fixef$date <- as.Date(paste(preds.fixef$year, ifelse(nchar(as.character(preds.fixef$month))==1, paste0("0", preds.fixef$month), as.character(preds.fixef$month)), "15", sep="-"))
  preds.fixef <- left_join(ds[which(!duplicated(paste0(ds$DATE_DEC, ds$sub_id))),], preds.fixef)
  
  finish <- Sys.time()
  runtime <- round(as.numeric(difftime(finish, start, units = "hours")),3)
  
  save(mod, samples, preds, preds.fixef, runtime,
       file = file.path(PATH_PROCESSED, "mods", paste0(det, "_", area, ".RData")))
  gc()
  NULL
}

fit.mod(task.id)

