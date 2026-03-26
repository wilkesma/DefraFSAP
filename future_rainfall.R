source(here::here("paths.R"))

####Packages
library(terra)
library(tidyterra)
library(ggplot2)
library(gridExtra)
library(stringr)
library(sf)

#Data overview - https://catalogue.ceda.ac.uk/uuid/8194b416cbee482b89e0dfbe17c5786c/
#We used the standard member (01)
#e.g. rcp2.6 - https://data.ceda.ac.uk/badc/deposited2021/chess-scape/data/rcp26_bias-corrected/01/monthly

##First compare the CHESS-SCAPE and HadUKGrid annual means for the study period 2003-2023
england <- st_read(PATH_COUNTRIES)
england <- st_transform(england[england$ctry17nm=="England",], 27700)

get_diff <- function(rcp){
  pr <- rast(chess_scape_path(rcp))
  
  had_files <- sapply(2003:2023, haduk_annual_path)
  
  # Load and combine all HadUKGrid data
  had_all <- lapply(had_files, rast)
  had_all <- rast(had_all)
  
  #Mask
  pr <- crop(mask(pr, england), england)
  had_all <- crop(mask(had_all, england), england)
  
  # Get corresponding CHESS-SCAPE data (2003-2023, all months)
  pr_dates <- as.Date(time(pr))
  pr_years <- as.numeric(format(pr_dates, "%Y"))
  idx_2003_2023 <- which(pr_years >= 2003 & pr_years <= 2023)
  pr_2003_2023 <- pr[[idx_2003_2023]]
  
  # Convert to mm
  pr_2003_2023_mm <- pr_2003_2023 * 30 * 24 * 60 * 60
  units(pr_2003_2023_mm) <- "mm"
  
  # Resample to HadUKGrid
  pr_2003_2023_resampled <- project(pr_2003_2023_mm, had_all[[1]], method = "bilinear")
  
  # Calculate annual totals for both datasets
  # CHESS-SCAPE: 12 months per year
  chess_annual <- tapp(pr_2003_2023_resampled, 
                       index = rep(2003:2023, each = 12), 
                       fun = sum)
  
  # HadUKGrid: group by year
  had_times <- time(had_all)
  had_years_all <- as.numeric(format(as.Date(had_times), "%Y"))
  had_annual <- tapp(had_all, 
                     index = had_years_all, 
                     fun = sum)
  
  # Mean annual difference
  diff_annual <- chess_annual - had_annual
  mean_diff <- mean(diff_annual)
  
  # Overall mean difference across all cells and years
  overall_mean_diff <- global(mean_diff, "mean", na.rm = TRUE)[1,1]
  
  # Plot
  plotdata <- c(mean_diff)
  names(plotdata) <- "mean_annual_diff"
  
  p <- ggplot() +
    geom_spatraster(data = plotdata) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = 0, na.value = "transparent",
                         name = "Difference (mm)") +
    theme_void() +
    labs(title = paste0("RCP", rcp/10, " mean annual difference 2003-2023"),
         subtitle = paste0("Overall mean: ", round(overall_mean_diff, 1), " mm/year"))
  print(p)
  
  cat("\nRCP", rcp/10, "mean offset:", round(overall_mean_diff, 1), "mm/year\n")
  cat("Summary of spatial variation:\n")
  print(summary(values(mean_diff), na.rm = TRUE))
  
  return(overall_mean_diff)
}

#Calculate offsets for each RCP
offset_rcp26 <- get_diff(26)
offset_rcp60 <- get_diff(60)
offset_rcp85 <- get_diff(85) #All offsets are close

##Get the data for 12 months to end of June 2030 and 2042
#Load and wrangle long-term monthly means from HadUK grid data
long.files <- c(PATH_HADUK_30Y_1961_1990, PATH_HADUK_30Y_1991_2020)
long.rasts <- lapply(long.files, rast) #Two raster stacks of monthly values covering 1961-2020 inclusive
long.rasts.agg <- list()
for(i in 1:12){
  long.rasts.agg[[i]] <- app(c(long.rasts[[1]][[i]], long.rasts[[2]][[i]]), mean)
}
long.rasts <- do.call(c, long.rasts.agg) #Now a single raster stack of long term monthly means
long.rasts <- project(long.rasts, "EPSG:27700")
rm(long.rasts.agg)

#Functon to get 2030 and 2042 anomalies for a single RCP
get_anoms <- function(rcp){
  #Get the RCP name for later
  rcp_pt <- sub("(\\d)(\\d)$", "\\1.\\2", as.character(rcp))
  
  #Read in thew corresponding raster
  pr <- rast(chess_scape_path(rcp))
  
  # Get the time information
  pr_times <- time(pr)
  pr_dates <- as.Date(pr_times)
  
  #Find the layer that corresponds to June 2030 and June 2042
  june_2030_idx <- which(format(pr_dates, "%Y-%m") == "2030-06")[1]
  june_2042_idx <- which(format(pr_dates, "%Y-%m") == "2042-06")[1]
  
  #Get 12 months ENDING in June (i.e., previous 11 months + June)
  idx_2030 <- (june_2030_idx - 11):june_2030_idx
  idx_2042 <- (june_2042_idx - 11):june_2042_idx
  
  #Extract the layers
  pr_2030 <- pr[[idx_2030]]
  pr_2042 <- pr[[idx_2042]]
  
  #Convert units
  pr_2030 <- pr_2030 * 30 * 24 * 60 * 60 #Seconds per 30-day month
  units(pr_2030) <- "mm"
  pr_2042 <- pr_2042 * 30 * 24 * 60 * 60 #Seconds per 30-day month
  units(pr_2042) <- "mm"
  
  ##Calculate anomalies
  #Resample CHESS-SCAPE data
  pr_2030 <- project(pr_2030, long.rasts[[1]], method = "bilinear")
  pr_2042 <- project(pr_2042, long.rasts[[1]], method = "bilinear")
  
  #Final calculation
  anom_2030 <- list()
  anom_2042 <- list()
  for(i in 1:12){
    message(i)
    anom_2030[[i]] <-  pr_2030[[i]]-long.rasts[[i]]
    anom_2042[[i]] <-  pr_2042[[i]]-long.rasts[[i]]
  }
  anom_2030 <- do.call(c, anom_2030)
  names(anom_2030) <- month.name[c(7:12, 1:6)]
  anom_2042 <- do.call(c, anom_2042)
  names(anom_2042) <- month.name[c(7:12, 1:6)]
  
  pdf(paste0("Future rainfall monthly anomalies RCP ", rcp_pt, ".pdf"), height=11, width=7)
  grid.arrange(
    ggplot() +
      geom_spatraster(data=crop(anom_2030, st_buffer(england, 25000))) +
      facet_wrap(~lyr) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0, na.value = "transparent",
                           name = "Anomaly (mm)") +
      coord_sf(expand=FALSE) +
      labs(title="2030") +
      theme_dark(),
    ggplot() +
      geom_spatraster(data=crop(anom_2042, st_buffer(england, 25000))) +
      facet_wrap(~lyr) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0, na.value = "transparent",
                           name = "Anomaly (mm)") +
      coord_sf(expand=FALSE) +
      labs(title="2042") +
      theme_dark(),
    ncol=1)
  dev.off()
  
  ##Calculate mean anomalies and plot
  out_2030 <- app(anom_2030, mean)
  out_2042 <- app(anom_2042, mean)
  
  pdf(paste0("Future rainfall annual anomalies RCP ", rcp_pt, ".pdf"), height=4, width=7)
  grid.arrange(
    ggplot() +
      geom_spatraster(data=crop(out_2030, st_buffer(england, 25000))) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0, na.value = "transparent",
                           name = "Anomaly (mm)") +
      coord_sf(expand=FALSE) +
      labs(title="2030") +
      theme_dark(),
    ggplot() +
      geom_spatraster(data=crop(out_2042, st_buffer(england, 25000))) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                           midpoint = 0, na.value = "transparent",
                           name = "Anomaly (mm)") +
      coord_sf(expand=FALSE) +
      labs(title="2042") +
      theme_dark(),
    ncol=2)
  dev.off()
  
  #Save grids
  writeRaster(out_2030, file.path(PATH_PROCESSED, paste0("clim_2030_", rcp, ".tif")), overwrite=TRUE)
  writeRaster(out_2042, file.path(PATH_PROCESSED, paste0("clim_2042_", rcp, ".tif")), overwrite=TRUE)
}

lapply(c(26, 60, 85), get_anoms)  


