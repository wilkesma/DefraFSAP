source(here::here("paths.R"))

####Packages
library(sf)
library(ggplot2)
library(ggspatial)
library(terra)
library(stringr)
library(lubridate)
library(parallel)
library(dplyr)
library(tidyr)
library(readxl)

##Modelled data
invs_env <- readRDS(file.path(PATH_PROCESSED, "invs_env_data_complete_2003_2023.rds"))
fish_env <- readRDS(file.path(PATH_PROCESSED, "fish_env_data_complete_2003_2023.rds"))

##Scalings used for training data
vars_to_scale <- c(
  "PC1","PC2","CRI","ASR",
  "wT_12M","TIN_12M","PO4_12M","pH_12M","Cu_d_12M","Zn_d_12M",
  "wT_9M","TIN_9M","PO4_9M","pH_9M","Cu_d_9M","Zn_d_9M",
  "wT_6M","TIN_6M","PO4_6M","pH_6M","Cu_d_6M","Zn_d_6M",
  "wT_3M","TIN_3M","PO4_3M","pH_3M","Cu_d_3M","Zn_d_3M",
  "rainfall"
)

invs_scale <- invs_env %>%
  summarise(across(all_of(vars_to_scale),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE)))) %>%
  pivot_longer(
    everything(),
    names_to = c("variable", ".value"),
    names_pattern = "^(.*)_(mean|sd)$"
  )

fish_scale <- fish_env %>%
  summarise(across(all_of(vars_to_scale),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE)))) %>%
  pivot_longer(
    everything(),
    names_to = c("variable", ".value"),
    names_pattern = "^(.*)_(mean|sd)$"
  )

##Static variables - get these from training data
vars_static <- c("PC1","PC2","CRI", "CAMS","ASR", "sewage", "barrier_density")
fish_static <- unique(fish_env[,c("SITE_ID", vars_static)])
invs_static <- unique(invs_env[,c("SITE_ID", vars_static)])

##Rainfall
get.clim <- function(var){
  mon.files <- sapply(2002:2023, haduk_annual_path)
  mon.files <- mon.files[str_detect(mon.files, paste0(var, "_hadukgrid_uk_5km_mon_"))]
  long.files <- c(PATH_HADUK_30Y_1961_1990, PATH_HADUK_30Y_1991_2020)
  long.files <- long.files[str_detect(long.files, paste0(var, "_hadukgrid_uk_5km_mon-"))]
  mon.rasts <- lapply(mon.files, rast) #Monthly rasters, 1989-2024 inclusive
  names(mon.rasts) <- paste0("x", sapply(mon.files, function(x) substr(str_split(x, "_")[[1]][6], 1, 4)))
  long.rasts <- lapply(long.files, rast) #Two raster stacks of monthly values covering 1961-2020 inclusive
  long.rasts.agg <- list()
  for(i in 1:12){
    long.rasts.agg[[i]] <- app(c(long.rasts[[1]][[i]], long.rasts[[2]][[i]]), mean)
  }
  long.rasts <- do.call(c, long.rasts.agg) #Now a single raster stack of long term monthly means
  rm(long.rasts.agg)
  anom.rasts <- list()
  for(i in 1:12){
    anom.rasts[[i]] <-  lapply(mon.rasts, function(x) x[[i]]-long.rasts[[i]])
  }
  anom.rasts <- lapply(seq(1:12), function(x) project(do.call(c, lapply(names(mon.rasts), function(y) anom.rasts[[x]][[y]])), "EPSG:27700"))
  for(i in 1:12){
    names(anom.rasts[[i]]) <- names(mon.rasts)
  } #anom.rasts is a list of rasters, each element is a month and each band is a year
  anom.rasts
}

extract.clim.anom <- function(env, x){
  years.x <- 2003:2023
  months.x <- c(seq(6+1, 12, by=1), seq(1:6))
  times <- expand_grid(
    window_year = years.x,
    month = months.x
  ) %>%
    mutate(
      year = if_else(month >= 7, window_year - 1, window_year)
    ) %>%
    select(window_year, year, month) %>%
    as.data.frame()
  
  output <- do.call(rbind, lapply(years.x, function(y) {
    times.x <- times[times$window_year==y,]
    clim.x <- matrix(nrow=12, ncol=1)
    colnames(clim.x) <- "rainfall"
    for(i in 1:12){
      clim.x[i,] <- terra::extract(clim[[times.x$month[i]]][[paste0("x", times.x$year[i])]], st_as_sf(unique(env[env$SITE_ID==x,c("SITE_ID", "easting", "northing")]), coords=c("easting", "northing"), crs=27700))[,2]
    }
    if(any(is.na(clim.x))){
      for(i in 1:12){
        clim.x[i,] <- mean(terra::extract(clim[[times.x$month[i]]][[paste0("x", times.x$year[i])]], st_buffer(st_as_sf(unique(env[env$SITE_ID==x,c("SITE_ID", "easting", "northing")]), coords=c("easting", "northing"), crs=27700), 20000))[,2], na.rm=T) #For coastal meta falling just of the raster, calculate the mean anomalies in a 20 km buffer
      }
    }
    data.frame(SITE_ID=x, year=y, t(as.data.frame(colMeans(clim.x))))
  }))
  row.names(output) <- NULL
  output
}

clim <- get.clim("rainfall")

invs_clim <- do.call(rbind, mclapply(unique(invs_env$SITE_ID), function(x) extract.clim.anom(invs_env, x), mc.cores=50))
fish_clim <- do.call(rbind, mclapply(unique(fish_env$SITE_ID), function(x) extract.clim.anom(fish_env, x), mc.cores=50))

##Water quality
dets <- c("pH", "AmN", "NO2", "NO3", "PO4", "Cu_d", "Zn_d", "wT")
areas <- c("AN", "MD", "NE", "NW", "SO", "SW", "TH")

get.wq <- function(det, area){
  message(det, area)
  load(file.path(PATH_PROCESSED, "mods", paste0(det, "_", area, ".RData"))) #mod, samples, preds, preds.fixef, runtime
  output <- data.frame(preds.fixef[,c("year", "month", "sub_id")], x=preds.fixef$Concentration)
  colnames(output)[4] <- det
  gc()
  output
}

wq_preds <- do.call(rbind, mclapply(areas, function(area) {
  Reduce(function(x, y) left_join(x, y), lapply(dets, function(det) get.wq(det, area)))
}, mc.cores = length(areas)))

invs_meta_match <- read.csv(file.path(PATH_PROCESSED, "invs_meta_matched.csv"), row.names=1)
fish_meta_match <- read.csv(file.path(PATH_PROCESSED, "fish_meta_matched.csv"), row.names=1)

extract.wq <- function(meta_match, x){
  years.x <- 2003:2023
  months.x <- c(seq(6+1, 12, by=1), seq(1:6))
  times <- expand_grid(
    window_year = years.x,
    month = months.x
  ) %>%
    mutate(
      year = if_else(month >= 7, window_year - 1, window_year)
    ) %>%
    select(window_year, year, month) %>%
    as.data.frame()

  preds.x <- wq_preds[wq_preds$sub_id==meta_match$sub_id[meta_match$SITE_ID==x][1],]
  if(nrow(preds.x)>0){
    output <- do.call(rbind, lapply(years.x, function(y) {
      times.x <- times[times$window_year==y,]
      wq.x <- matrix(nrow=12, ncol=length(dets))
      colnames(wq.x) <- dets
      for(i in 1:12){
        wq.x[i,] <- as.numeric(colMeans(preds.x[preds.x$year==y & preds.x$month==times.x$month[i], dets]))
      }
      data.frame(SITE_ID=x, year=y, t(as.data.frame(colMeans(wq.x))))
    }))
    row.names(output) <- NULL
  } else{
    output <- df <- data.frame(SITE_ID = x, year = years.x)
    output[dets] <- NA
  }
  colnames(output) <- c("SITE_ID", "year", paste0(dets, "_12M"))
  output
}

invs_wq <- do.call(rbind, mclapply(unique(invs_env$SITE_ID), function(x) extract.wq(invs_meta_match, x), mc.cores=50))
invs_wq$TIN_12M <- rowSums(invs_wq[,c("AmN_12M", "NO2_12M", "NO3_12M")])
invs_wq <- invs_wq[,-which(substr(colnames(invs_wq), 1, 3) %in% paste0(c("AmN", "NO2", "NO3")))]

fish_wq <- do.call(rbind, mclapply(unique(fish_env$SITE_ID), function(x) extract.wq(fish_meta_match, x), mc.cores=50))
fish_wq$TIN_12M <- rowSums(fish_wq[,c("AmN_12M", "NO2_12M", "NO3_12M")])
fish_wq <- fish_wq[,-which(substr(colnames(fish_wq), 1, 3) %in% paste0(c("AmN", "NO2", "NO3")))]

##Joining
invs_data <- left_join(left_join(invs_static, invs_clim), invs_wq)
fish_data <- left_join(left_join(fish_static, fish_clim), fish_wq)

##Saving
saveRDS(list(basedata=invs_data, scales=invs_scale), file.path(PATH_PROCESSED, "invs_historical_pred_data.rds"))
saveRDS(list(basedata=fish_data, scales=fish_scale), file.path(PATH_PROCESSED, "fish_historical_pred_data.rds"))


