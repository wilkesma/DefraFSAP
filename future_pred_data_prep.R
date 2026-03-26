source(here::here("paths.R"))

####Packages
library(terra)
library(sf)
library(dplyr)
library(parallel)

##Load data
#Training data
invs <- readRDS(file.path(PATH_PROCESSED, "invs_historical_pred_data.rds")) #list(basedata, scales)
fish <- readRDS(file.path(PATH_PROCESSED, "fish_historical_pred_data.rds")) #list(basedata, scales)

#Original data to get waterbody ID and coordinates
invs_meta <- readRDS(file.path(PATH_PROCESSED, "invs_env_data_complete_2003_2023.rds"))
fish_meta <- readRDS(file.path(PATH_PROCESSED, "fish_env_data_complete_2003_2023.rds"))

#List of HMWB and Artificial waterbody IDs
hmwbs <- read.csv(PATH_WFD_CLASS) #https://data.catchmentbasedapproach.org/maps/4cff14f14acc41d2ac7c2a6b1e200c10
hmwbs <- hmwbs[which(hmwbs$Heavily.Modified.Designation %in% c("Artificial", "Heavily Modified")),]

##Flagging HMWB/Artifical SITE_IDs
#Get HMWB/Artificial SITE_IDs
invs_hmwbs <- invs_meta$SITE_ID[which(invs_meta$WFD_WATERBODY_ID %in% hmwbs$EA.Waterbody.ID)]
fish_hmwbs <- fish_meta$SITE_ID[which(fish_meta$GEO_WATERBODY %in% hmwbs$EA.Waterbody.ID)]

length(unique(invs_hmwbs))/length(unique(invs_meta$SITE_ID)) #29% HMWBs
length(unique(fish_hmwbs))/length(unique(fish_meta$SITE_ID)) #37% HMWBs

#Make data.frame to add mining flag and rainfall anomalies to
invs_fut <- data.frame(SITE_ID=unique(invs$basedata$SITE_ID), HMWB_A=unique(invs$basedata$SITE_ID) %in% invs_hmwbs)
sum(invs_fut$HMWB_A)/nrow(invs_fut) #29% HMWBs
fish_fut <- data.frame(SITE_ID=unique(fish$basedata$SITE_ID), HMWB_A=unique(fish$basedata$SITE_ID) %in% fish_hmwbs)
sum(fish_fut$HMWB_A)/nrow(fish_fut) #27% HMWBs

##Flagging areas affected by mining
invs_fut <- left_join(invs_fut, unique(invs_meta[,c("SITE_ID", "mining")]))
fish_fut <- left_join(fish_fut, unique(fish_meta[,c("SITE_ID", "mining")]))

##Rainfall anomalies to June 2030 and 2042 inclusive under three RCPs
#Load rasters
clim <- lapply(c(2030, 2042), function(y) lapply(c(26, 60, 85), function(rcp) rast(file.path(PATH_PROCESSED, paste0("clim_", y, "_", rcp, ".tif")))))
names(clim) <- paste0("x", c(2030, 2042))
for(i in 1:2){ names(clim[[i]]) <- paste0("RCP", c(26, 60, 85)) }

#Convert unique SITE_ID locations to sf format
invs_sf <- st_as_sf(unique(invs_meta[,c("SITE_ID", "easting", "northing")]), coords=c("easting", "northing"), crs=27700)
fish_sf <- st_as_sf(unique(fish_meta[,c("SITE_ID", "easting", "northing")]), coords=c("easting", "northing"), crs=27700)

#Function to extract anomalies from rasters
get_anom <- function(df, y, rcp){
  output <- extract(clim[[paste0("x", y)]][[paste0("RCP", rcp)]], df, ID=FALSE)
  which_na <- which(is.na(output$mean))
  for(i in which_na){
    output[i,] <- mean(extract(clim[[paste0("x", y)]][[paste0("RCP", rcp)]],
                                     st_buffer(df[i,], 20000), ID=FALSE)[,1], na.rm=TRUE) #For coastal meta falling just of the raster, calculate the mean anomalies in a 20 km buffer
  }
  colnames(output) <- paste0("rainfall_", y, "_", rcp)
  output
}

invs_clim <- do.call(cbind, lapply(c(2030, 2042), function(y) do.call(cbind, lapply(c(26, 60, 85), function(rcp) get_anom(invs_sf, y, rcp)))))
fish_clim <- do.call(cbind, lapply(c(2030, 2042), function(y) do.call(cbind, lapply(c(26, 60, 85), function(rcp) get_anom(fish_sf, y, rcp)))))

nrow(invs_clim); nrow(na.omit(invs_clim))
nrow(fish_clim); nrow(na.omit(fish_clim))

invs_clim$SITE_ID <- invs_sf$SITE_ID
fish_clim$SITE_ID <- fish_sf$SITE_ID

##Join HMWBs and clim data together
invs_fut <- left_join(invs_fut, invs_clim)
fish_fut <- left_join(fish_fut, fish_clim)

write.csv(invs_fut, file.path(PATH_PROCESSED, "invs_future_pred_data.csv"), row.names=FALSE)
write.csv(fish_fut, file.path(PATH_PROCESSED, "fish_future_pred_data.csv"), row.names=FALSE)
