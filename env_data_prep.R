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

##Data
rivers <- st_read(PATH_OS_RIVERS)
meta_snap <- st_read(PATH_META_SNAP)

load(file.path(PATH_PROCESSED, "invs_data.RData")) #invs_abun.sc2, invs_info.sc2, invs_meta
rm(invs_abun.sc2, invs_info.sc2) #We only want invs_meta

load(file.path(PATH_PROCESSED, "fish_data.RData")) #fish_abun, fish_meta
rm(fish_abun) #We only want fish_meta

##Sewage (water companies)
stws_snap <- st_read(PATH_STWS_SNAP)
stws_snap <- stws_snap[which(stws_snap$DSI_TYPE_D %in% c("WwTW/Sewage Treatment Works (water company)", "Storm Tank/CSO on Sewerage Network (water company)")),]

meta_join <- st_join(meta_snap, rivers, join=st_nearest_feature)
stws_join <- st_join(stws_snap, rivers, join=st_nearest_feature)

meta_join$sewage <- meta_join$identifier %in% stws_join$identifier

pdf("Sewage linking test.pdf")
ggplot() +
  geom_sf(data=st_intersection(rivers, st_buffer(st_union(meta_join[meta_join$name1=="River Aln",]), 5000)), colour="blue") +
  geom_sf(data=meta_join[which(meta_join$name1=="River Aln"),], aes(colour=sewage)) +
  scale_colour_manual(values=c("darkgrey", "darkgreen")) +
  geom_sf(data=st_intersection(stws_join, st_buffer(st_union(meta_join[meta_join$name1=="River Aln",]), 5000)), colour="red") +
  theme_light() #Tested - works fine
dev.off()

all(meta_snap$SITE_ID==meta_join$SITE_ID) #Fine to just add sewage column to meta_snap
meta_snap$sewage <- meta_join$sewage

##RICT PCA axes
rict <- rast(PATH_RICT)

rict_points <- data.frame(SITE_ID=meta_snap$SITE_ID, terra::extract(rict, meta_snap, ID=FALSE))
nrow(rict_points[is.na(rict_points$PC1),])/nrow(rict_points)*100 #30% are NA

rict_points_na_shp <- meta_snap[which(meta_snap$SITE_ID %in% rict_points$SITE_ID[is.na(rict_points$PC1)]),]
for(i in row.names(rict_points_na_shp)){
  rict_points[i, c("PC1", "PC2")] <- colMeans(terra::extract(rict, st_buffer(rict_points_na_shp[i,], 100), ID=FALSE), na.rm=TRUE)
}

nrow(rict_points[is.na(rict_points$PC1),])/nrow(rict_points)*100 #2.2% are still NA

rict_points_na_shp <- meta_snap[which(meta_snap$SITE_ID %in% rict_points$SITE_ID[is.na(rict_points$PC1)]),]
for(i in row.names(rict_points_na_shp)){
  rict_points[i, c("PC1", "PC2")] <- colMeans(terra::extract(rict, st_buffer(rict_points_na_shp[i,], 250), ID=FALSE), na.rm=TRUE)
}

nrow(rict_points[is.na(rict_points$PC1),])/nrow(rict_points)*100 #0.82% are still NA

rict_points_na_shp <- meta_snap[which(meta_snap$SITE_ID %in% rict_points$SITE_ID[is.na(rict_points$PC1)]),]
for(i in row.names(rict_points_na_shp)){
  rict_points[i, c("PC1", "PC2")] <- colMeans(terra::extract(rict, st_buffer(rict_points_na_shp[i,], 500), ID=FALSE), na.rm=TRUE)
}

nrow(rict_points[is.na(rict_points$PC1),])/nrow(rict_points)*100 #0.34% are still NA - OK

meta_snap$PC1 <- rict_points$PC1
meta_snap$PC2 <- rict_points$PC2

##Rainfall
get.clim <- function(var){
  mon.files <- sapply(1989:2024, haduk_annual_path)
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

clim <- get.clim("rainfall")

get.clim.plotdata <- function(month){
  output <- data.frame(month=month, year=seq(1989, 2024), var="rainfall", t(apply(na.omit(values(clim[[month]])), 2, function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975)))))
  colnames(output)[4:6] <- c("mean", "lwr", "upr")
  output$date <- as.Date(paste("15", as.character(output$month), as.character(output$year), sep="/"), "%d/%m/%Y")
  output
}
clim.plotdata <- do.call(rbind, lapply(seq(1:12), get.clim.plotdata))
clim.plotdata$var <- factor(clim.plotdata$var)
levels(clim.plotdata$var) <-"Monthly rainfall anomaly (mm)"

pdf("Monthly UK rainfall anomalies 1989 to 2024.pdf", width=6, height=3)
ggplot(clim.plotdata, aes(x=date, y=mean)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.35) +
  geom_line() + ylab("Monthly mean rainfall anomaly (mm)") +
  scale_x_date(date_breaks="5 year", date_minor_breaks="12 month", date_labels="%Y", expand =c(0,0)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.title.x=element_blank())
dev.off()

extract.clim.anom <- function(meta, x, id_column, date_column){
  date.x <- as.character(meta[date_column][meta[id_column]==x])
  month.x <- month(round_date(as.Date(date.x, "%d/%m/%Y"), unit="month"))-1 #-1 to give inclusive month
  years.x <- year(as.Date(date.x, "%d/%m/%Y"))
  if(month.x!=12){
    months.x <- c(seq(month.x+1, 12, by=1), seq(1:month.x))
    years.x <- c(rep(years.x-1, length(months.x[1:which(months.x==1)-1])), rep(years.x, month.x))
  } else{
    months.x <- seq(1:12)
  }
  times <- data.frame(month=months.x, year=years.x)
  clim.x <- matrix(nrow=12, ncol=1)
  colnames(clim.x) <- "rainfall"
  for(i in 1:12){
    clim.x[i,] <- terra::extract(clim[[times$month[i]]][[paste0("x", times$year[i])]], st_as_sf(meta[meta[id_column]==x,], coords=c("easting", "northing"), crs=27700))[,2]
  }
  if(any(is.na(clim.x))){
    for(i in 1:12){
      clim.x[i,] <- mean(terra::extract(clim[[times$month[i]]][[paste0("x", times$year[i])]], st_buffer(st_as_sf(meta[meta[id_column]==x,], coords=c("easting", "northing"), crs=27700), 20000))[,2], na.rm=T) #For coastal meta falling just of the raster, calculate the mean anomalies in a 20 km buffer
    }
  }
  output <- data.frame(ANALYSIS_ID=x, t(as.data.frame(colMeans(clim.x))))
  colnames(output)[1] <- id_column
  output
}

invs_clim_points <- do.call(rbind, mclapply(invs_meta$ANALYSIS_ID, function(x) extract.clim.anom(invs_meta, x, "ANALYSIS_ID", "SAMPLE_DATE"), mc.cores=50))
fish_clim_points <- do.call(rbind, mclapply(fish_meta$SURVEY_ID, function(x) extract.clim.anom(fish_meta, x, "SURVEY_ID", "EVENT_DATE"), mc.cores=50))
#NOTE - Would have been more efficient to do this on *_meta_match

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

extract.wq <- function(meta, meta_match, x, id_column, date_column){
  date.x <- as.character(meta[date_column][meta[id_column]==x])
  month.x <- month(round_date(as.Date(date.x, "%d/%m/%Y"), unit="month"))-1 #-1 to give inclusive month
  years.x <- year(as.Date(date.x, "%d/%m/%Y"))
  if(month.x!=12){
    months.x <- c(seq(month.x+1, 12, by=1), seq(1:month.x))
    years.x <- c(rep(years.x-1, length(months.x[1:which(months.x==1)-1])), rep(years.x, month.x))
  } else{
    months.x <- seq(1:12)
  }
  times <- c()
  for(i in 1:12) { times[i] <- paste0(years.x[i], months.x[i]) }
  
  site.x <- meta$SITE_ID[meta[id_column]==x]
  preds.x <- wq_preds[wq_preds$sub_id==meta_match$sub_id[meta_match$SITE_ID==site.x][1],]
  if(nrow(preds.x)>0){
    preds.x$time <- NA
    for(i in 1:nrow(preds.x)) { preds.x$time[i] <- paste0(preds.x$year[i], preds.x$month[i]) }
    preds.x <- preds.x[preds.x$time %in% times, -which(colnames(preds.x)=="time")]
    
    preds.12m <- colMeans(preds.x[,dets])
    names(preds.12m) <- paste(dets, "12M", sep="_")
    preds.9m <- colMeans(preds.x[4:12,dets])
    names(preds.9m) <- paste(dets, "9M", sep="_")
    preds.6m <- colMeans(preds.x[7:12,dets])
    names(preds.6m) <- paste(dets, "6M", sep="_")
    preds.3m <- colMeans(preds.x[10:12,dets])
    names(preds.3m) <- paste(dets, "3M", sep="_")
    
    output <- as.data.frame(t(c(id=x, preds.12m, preds.9m, preds.6m, preds.3m)))
    colnames(output)[1] <- id_column
  } else{
    output <- as.data.frame(t(c(id=x, rep(NA, length(dets)*4))))
    colnames(output) <- c(id_column, c(paste(dets, "12M", sep="_"), paste(dets, "9M", sep="_"), paste(dets, "6M", sep="_"), paste(dets, "3M", sep="_")))
  }
  
  output
}

invs_wq <- do.call(rbind, mclapply(invs_meta_match$ANALYSIS_ID[as.numeric(invs_meta_match$year)>2002 & as.numeric(invs_meta_match$year)<2024], function(x) extract.wq(invs_meta, invs_meta_match, x, "ANALYSIS_ID", "SAMPLE_DATE"), mc.cores=50))
fish_wq <- do.call(rbind, mclapply(fish_meta_match$SURVEY_ID[as.numeric(fish_meta_match$year)>2002 & as.numeric(fish_meta_match$year)<2024], function(x) extract.wq(fish_meta, fish_meta_match, x, "SURVEY_ID", "EVENT_DATE"), mc.cores=50))

invs_wq$TIN_12M <- rowSums(invs_wq[,c("AmN_12M", "NO2_12M", "NO3_12M")])
invs_wq$TIN_9M <- rowSums(invs_wq[,c("AmN_9M", "NO2_9M", "NO3_9M")])
invs_wq$TIN_6M <- rowSums(invs_wq[,c("AmN_6M", "NO2_6M", "NO3_6M")])
invs_wq$TIN_3M <- rowSums(invs_wq[,c("AmN_3M", "NO2_3M", "NO3_3M")])
invs_wq <- invs_wq[,-which(substr(colnames(invs_wq), 1, 3) %in% paste0(c("AmN", "NO2", "NO3")))]

fish_wq$TIN_12M <- rowSums(fish_wq[,c("AmN_12M", "NO2_12M", "NO3_12M")])
fish_wq$TIN_9M <- rowSums(fish_wq[,c("AmN_9M", "NO2_9M", "NO3_9M")])
fish_wq$TIN_6M <- rowSums(fish_wq[,c("AmN_6M", "NO2_6M", "NO3_6M")])
fish_wq$TIN_3M <- rowSums(fish_wq[,c("AmN_3M", "NO2_3M", "NO3_3M")])
fish_wq <- fish_wq[,-which(substr(colnames(fish_wq), 1, 3) %in% paste0(c("AmN", "NO2", "NO3")))]

nrow(invs_wq) #59742
nrow(na.omit(invs_wq)) #21770 complete observations
apply(invs_wq[,-1], 2, function(x) length(x[is.na(x)])/length(x)*100) #From 1% (pH) to 33% (Cu_d) and 62% (Zn_d) missingness
#Other dets ~2-4% missing

nrow(fish_wq) #6187
nrow(na.omit(fish_wq)) #1980 complete observations
apply(fish_wq[,-1], 2, function(x) length(x[is.na(x)])/length(x)*100) #From 1.5% (pH) to 38% (Cu_d) and 65% (Zn_d) missingness
#Other dets ~2-8% missing

##Agricultural sediment risk & Channel Resectioning Index
asr_cri <- read_excel(PATH_ASR_CRI)
asr_cri <- st_as_sf(
  asr_cri,
  coords = c("X", "Y"),
  crs = 27700
)
nearest_id <- st_nearest_feature(meta_snap, asr_cri)

cri.dists <- sapply(1:nrow(meta_snap), function(x) {message(x); st_distance(meta_snap[x,], asr_cri[nearest_id[x],])})
summary(cri.dists) #Mean=122 m; 3rd quartile=172 m; max=3.2 km
quantile(cri.dists, 0.95) #95th percentile=223 m
length(cri.dists[cri.dists>500]) #Only 112 sites are >500 m from the CRI point - this is ok

pdf("CRI snapping distance ECDF.pdf", width=6.7/2, height=2.5)
ggplot(data.frame(x=cri.dists), aes(x=x)) +
  stat_ecdf() +
  scale_x_log10("Snapping distance (m)", breaks=c(1, 10, 100, 1000), labels=c("1", "10", "100", "1,000")) +
  ylab("ECDF")
dev.off()

meta_snap$ASR <- as.integer(asr_cri$ASR[nearest_id])
meta_snap$CRI <- as.integer(asr_cri$CRI[nearest_id])

##AMBER riverine barrier density
catchments <- st_read(PATH_WFD_CATCHMENTS)
barriers <- read.csv(PATH_BARRIERS)

barriers_sf <- st_as_sf(
  barriers,
  coords = c("Longitude_WGS84", "Latitude_WGS84"),
  crs = 4326
) |> 
  st_transform(27700)

barriers_catch <- st_join(barriers_sf, catchments)

barrier_counts <- barriers_catch |>
  st_drop_geometry() |>
  count(wb_id, name = "n_barriers")

catchments$area_km2 <- as.numeric(st_area(catchments)) / 1e6

catchments <- catchments |>
  left_join(barrier_counts, by = "wb_id") |>
  mutate(
    n_barriers = ifelse(is.na(n_barriers), 0, n_barriers),
    barrier_density = n_barriers / area_km2
  )

meta_snap <- st_join(
  meta_snap,
  catchments["barrier_density"]
)

##Locations affected by abandoned mine drainage
mines <- read_excel(
  PATH_MINES,
  sheet = "WAMM_list_Nov_2019",
  skip = 1  # adjust so the first real header row becomes row 1 in R
)

catchments$mining <- catchments$wb_id %in% mines$EA_WB_ID

meta_snap <- st_join(
  meta_snap,
  catchments["mining"]
)
meta_snap$mining[is.na(meta_snap$mining)] <- FALSE #NA values were just not on catchments

##CAMS flow compliance bands
flow_comp <- st_read(PATH_FLOW_COMP)
nearest_idx <- st_nearest_feature(meta_snap, flow_comp)

flow_comp.dists <- sapply(1:nrow(meta_snap), function(x) {message(x); st_distance(meta_snap[x,], flow_comp[nearest_idx[x],])})
summary(flow_comp.dists) #Mean=146 m; 3rd quartile=4.7 Km; max=29.3 km!
quantile(flow_comp.dists, 0.95) #95th percentile=1 Km
length(flow_comp.dists[flow_comp.dists>1000]) #1306 sites are >5 km from the CAMS point
#Larger distances are small streams that were not assessed for CAMS and it wouldn't be safe to assume anything
#So we must consider anything beyond a small distance (e.g. 50 m) as "Not Assessed"

pdf("CAMS snapping distance ECDF.pdf", width=6.7/2, height=2.5)
ggplot(data.frame(x=flow_comp.dists), aes(x=x)) +
  stat_ecdf() +
  scale_x_log10("Snapping distance (m)", breaks=c(0.01, 10, 10000), labels=c("0.01", "10", "10,000")) +
  ylab("ECDF")
dev.off()

meta_snap$CAMS <- flow_comp$compliance[nearest_idx]
meta_snap$CAMS[which(flow_comp.dists>50)] <- "Not Assessed"
table(meta_snap$CAMS)

##Join all above back to raw meta data (invs_meta and fish_meta)
#Static variables (sewage, RICT PCs, ASR, mine drainage, CAMS, CRI, barrier density) - these are on meta_snap
meta_snap_df <- st_drop_geometry(meta_snap)
colnames(meta_snap_df)[2] <- "SITE_ID2" #To match correct column in meta_match for joining

invs_meta_match <- left_join(invs_meta_match, meta_snap_df[which(str_detect(meta_snap_df$SITE_ID2, "invs")),c("SITE_ID2", "PC1", "PC2", "sewage", "ASR", "CRI", "mining", "CAMS", "barrier_density")])
fish_meta_match <- left_join(fish_meta_match, meta_snap_df[which(str_detect(meta_snap_df$SITE_ID2, "fish")),c("SITE_ID2", "PC1", "PC2", "sewage", "ASR", "CRI", "mining", "CAMS", "barrier_density")])

#Rainfall anomalies - these are on raw meta (not matched) so just take those in *_meta_match
invs_meta_match <- left_join(invs_meta_match, invs_clim_points) #Discards data for any ANALYSIS_ID not in invs_meta_match
fish_meta_match <- left_join(fish_meta_match, fish_clim_points) #Same for fish SURVEY_ID

#Water quality - on *_meta_match (only observations between 2003 and 2023 inclusive)
invs_meta_match <- left_join(invs_meta_match, invs_wq)
fish_meta_match <- left_join(fish_meta_match, fish_wq)

##Missingness analysis
vars <- c("PC1", "PC2", "sewage", "rainfall", "ASR", "CRI", "mining", "CAMS", "barrier_density",
          paste0(c("pH", "TIN", "PO4", "Cu_d", "Zn_d", "wT"), "_12M"))

invs_miss <- invs_meta_match %>%
  filter(year >= 2003 & year <= 2023) %>%
  select(year, basin, all_of(vars)) %>%
  mutate(across(all_of(vars), as.character)) %>%  # Convert all to character
  pivot_longer(cols = all_of(vars), 
               names_to = "variable", 
               values_to = "value") %>%
  group_by(year, basin, variable) %>%
  summarise(
    prop_present = mean(!is.na(value)),
    n_total = dplyr::n(),
    n_present = sum(!is.na(value)),
    .groups = "drop"
  )

#Some var x basin x year combinations are missing - add them
all_combinations <- expand.grid(
  year = 2003:2023,
  basin = unique(invs_meta_match$basin),
  variable = vars,
  stringsAsFactors = FALSE
)
missing_combinations <- all_combinations %>%
  anti_join(invs_miss, by = c("year", "basin", "variable")) %>%
  mutate(prop_present = 1, n_total = 0, n_present = 0)
invs_miss <- bind_rows(invs_miss, missing_combinations)

invs_miss %>%
  group_by(variable) %>%
  summarise(
    overall_prop_present = mean(prop_present),
    min_year = min(year[prop_present > 0], na.rm = TRUE),
    max_year = max(year[prop_present > 0], na.rm = TRUE)
  ) %>%
  arrange(desc(overall_prop_present))

vars <- vars[-which(vars %in% c("PC1", "PC2", "sewage", "rainfall", "CAMS", "CRI", "barrier_density", "mining"))] #These have 100% coverage or close to 100% anyway
invs_miss <- invs_miss[invs_miss$variable %in% vars,]
invs_miss$variable <- factor(invs_miss$variable, levels=vars)

vars <- c("PC1", "PC2", "sewage", "rainfall", "ASR", "CRI", "mining", "CAMS", "barrier_density",
          paste0(c("pH", "TIN", "PO4", "Cu_d", "Zn_d", "wT"), "_12M"))

fish_miss <- fish_meta_match %>%
  filter(year >= 2003 & year <= 2023) %>%
  select(year, basin, all_of(vars)) %>%
  mutate(across(all_of(vars), as.character)) %>%  # Convert all to character
  pivot_longer(cols = all_of(vars), 
               names_to = "variable", 
               values_to = "value") %>%
  group_by(year, basin, variable) %>%
  summarise(
    prop_present = mean(!is.na(value)),
    n_total = dplyr::n(),
    n_present = sum(!is.na(value)),
    .groups = "drop"
  )

#Some var x basin x year combinations are missing - add them
all_combinations <- expand.grid(
  year = 2003:2023,
  basin = unique(fish_meta_match$basin),
  variable = vars,
  stringsAsFactors = FALSE
)
missing_combinations <- all_combinations %>%
  anti_join(fish_miss, by = c("year", "basin", "variable")) %>%
  mutate(prop_present = 1, n_total = 0, n_present = 0)
fish_miss <- bind_rows(fish_miss, missing_combinations)

fish_miss %>%
  group_by(variable) %>%
  summarise(
    overall_prop_present = mean(prop_present),
    min_year = min(year[prop_present > 0], na.rm = TRUE),
    max_year = max(year[prop_present > 0], na.rm = TRUE)
  ) %>%
  arrange(desc(overall_prop_present))

vars <- vars[-which(vars %in% c("PC1", "PC2", "sewage", "rainfall", "CAMS", "CRI", "barrier_density", "mining"))] #These have 100% coverage or close to 100% anyway
fish_miss <- fish_miss[fish_miss$variable %in% vars,]
fish_miss$variable <- factor(fish_miss$variable, levels=vars)

pdf("Missingness.pdf", width=10, height=7)
ggplot(invs_miss, aes(x = year, y = prop_present, fill = variable)) +
  geom_bar(stat="identity") +
  facet_grid(basin~variable) +
  scale_y_continuous(labels = scales::percent_format(), 
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = seq(2003, 2023, 5)) +
  labs(
    x = "Year",
    y = "Observations complete (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(face = "bold", size = 11, angle = 0, hjust = 0)
  )
ggplot(fish_miss, aes(x = year, y = prop_present, fill = variable)) +
  geom_bar(stat="identity") +
  facet_grid(basin~variable) +
  scale_y_continuous(labels = scales::percent_format(), 
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = seq(2003, 2023, 5)) +
  labs(
    x = "Year",
    y = "Observations complete (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(face = "bold", size = 11, angle = 0, hjust = 0)
  )
dev.off()

#ASR missing values are where dist from source couldn't be modelled; PSYCIC model maybe fell over

##Save all obs (including NAs - e.g. all WQ before 2003 and after 2023)
saveRDS(invs_meta_match, file.path(PATH_PROCESSED, "invs_env_data_full_1990_2024.rds"))
saveRDS(fish_meta_match, file.path(PATH_PROCESSED, "fish_env_data_full_1990_2024.rds"))

##Prep data frame for modelling - remove observations with any NAs for predictor variables EXCEPT Cu_d and Zn_d
#Vars we want complete cases for (i.e. not copper and zinc)
vars <- c("PC1", "PC2", "sewage", "rainfall", "ASR", "CRI", "mining", "CAMS", "barrier_density",
          paste0(c("pH", "TIN", "PO4", "wT"), "_12M"))

invs_meta_mod <- invs_meta_match[complete.cases(invs_meta_match[, vars]), ]
fish_meta_mod <- fish_meta_match[complete.cases(fish_meta_match[, vars]), ]

nrow(invs_meta_mod); nrow(invs_meta_match) #51239 of 123858 observations
nrow(fish_meta_mod); nrow(fish_meta_match) #4999 of 9338 observations

#Sort out classes
invs_meta_mod$year <- as.factor(invs_meta_mod$year)
invs_meta_mod$SITE_ID <- as.factor(invs_meta_mod$SITE_ID)
invs_meta_mod$CAMS <- factor(invs_meta_mod$CAMS, levels=c("COMPLIANT", "BAND1", "BAND2", "BAND3", "Not Assessed"))

fish_meta_mod$year <- as.factor(fish_meta_mod$year)
fish_meta_mod$SITE_ID <- as.factor(fish_meta_mod$SITE_ID)
fish_meta_mod$CAMS <- factor(fish_meta_mod$CAMS, levels=c("COMPLIANT", "BAND1", "BAND2", "BAND3", "Not Assessed"))

#CV fold membership
invs_site_n <- as.data.frame(table(invs_meta_mod$SITE_ID))
colnames(invs_site_n) <- c("SITE_ID", "n") #No frequency threshold for inclusion as INLA handles this well - will shrink singleton sites towards population mean for random effects

set.seed(123)

K <- 5
invs_site_n <- invs_site_n[order(-invs_site_n$n), ] #Helps ensure a balance of site observation frequency per fold
invs_site_n$fold <- NA
fold_totals <- rep(0, K)
for (i in seq_len(nrow(invs_site_n))) {
  f <- which.min(fold_totals)
  invs_site_n$fold[i] <- f
  fold_totals[f] <- fold_totals[f] + invs_site_n$n[i]
}
table(invs_site_n[,c("n", "fold")]) #Verify balance
invs_meta_mod <- left_join(invs_meta_mod, invs_site_n) #Join back to full data

fish_site_n <- as.data.frame(table(fish_meta_mod$SITE_ID))
colnames(fish_site_n) <- c("SITE_ID", "n") #No frequency threshold for inclusion as INLA handles this well - will shrink singleton sites towards population mean for random effects

fish_site_n <- fish_site_n[order(-fish_site_n$n), ] #Helps ensure a balance of site observation frequency per fold
fish_site_n$fold <- NA
fold_totals <- rep(0, K)
for (i in seq_len(nrow(fish_site_n))) {
  f <- which.min(fold_totals)
  fish_site_n$fold[i] <- f
  fold_totals[f] <- fold_totals[f] + fish_site_n$n[i]
}
table(fish_site_n[,c("n", "fold")])
fish_meta_mod <- left_join(fish_meta_mod, fish_site_n) #Join back to full data

saveRDS(invs_meta_mod, file.path(PATH_PROCESSED, "invs_env_data_complete_2003_2023.rds"))
saveRDS(fish_meta_mod, file.path(PATH_PROCESSED, "fish_env_data_complete_2003_2023.rds"))







