source(here::here("paths.R"))

####Packages
library(sf)
library(dplyr)
library(rnrfa)
library(shiny)
library(ggplot2)
library(lwgeom)
library(ggExtra)
library(gridExtra)
library(reshape2)

##Pre-preared biological data for inverts and fish
load(file.path(PATH_PROCESSED, "invs_data.RData")) #invs_abun.sc2, invs_info.sc2, invs_meta
rm(invs_abun.sc2, invs_info.sc2)
head(invs_meta)

load(file.path(PATH_PROCESSED, "fish_data.RData")) #fish_abun, fish_meta
rm(fish_abun)
head(fish_meta)

##Wrangled wq data - see /home/mw22803/ea_wq_project_v2/wq_data_prep.R
wq <- read.csv(file.path(PATH_PROCESSED, "ea_wq.csv"))
dets <- c("pH", "AmN", "NO2", "NO3", "PO4", "Cu_d", "Zn_d", "wT") #TIN is sum of N forms listed here

##Other data
#OS open rivers
rivers <- st_read(PATH_OS_RIVERS)

#Consented discharges
stws <- read.csv(PATH_CONSENTS)
#NOTE - Not necessarily active discharges during study period; includes all discharges, not just WWTWs

##Snapping
wq_sites <- unique(wq[,c("sample.samplingPoint.notation", "sample.samplingPoint.label", "easting", "northing")])
write.csv(wq_sites, file.path(PATH_PROCESSED, "ea_wq_sites.csv"))
invs_meta_sites <- unique(invs_meta[,c("SITE_ID", "easting", "northing")])
invs_meta_sites$SITE_ID <- paste0(invs_meta_sites$SITE_ID, "_invs")
fish_meta_sites <- unique(fish_meta[,c("SITE_ID", "easting", "northing")])
fish_meta_sites$SITE_ID <- paste0(fish_meta_sites$SITE_ID, "_fish")
meta_sites <- rbind(invs_meta_sites, fish_meta_sites)
write.csv(meta_sites, file.path(PATH_PROCESSED, "meta_sites.csv"))

nrow(meta_sites) #28635
nrow(wq_sites) #20054

#SAGA 'Snap points to lines' function with 100 m search distance - snaps to nearest point (inserting a new vertex where required)
#This left points further than 100 m away from OS rivers unsnapped, so then added buffer to OS open rivers and used this to clip snapped points
#QGIS process:
  #Snap points to lines
  #Clip to buffer
  #Save as shapefile

meta_snap <- st_read(PATH_META_SNAP)
wq_snap <- st_read(PATH_WQ_SNAP)
stws_snap <- st_read(PATH_STWS_SNAP)

nrow(meta_snap) #25116
nrow(invs_meta[which(paste0(invs_meta$SITE_ID, "_invs") %in% meta_snap$SITE_ID),]) #181646 invert samples
nrow(fish_meta[which(paste0(fish_meta$SITE_ID, "_fish") %in% meta_snap$SITE_ID),]) #20682 fish samples
nrow(wq_snap) #15932
nrow(wq[wq$sample.samplingPoint.notation %in% wq_snap$sample.sam,]) #~7.1M observations
nrow(stws_snap) #74277

##Detect river segments with at least one BIOSYS AND at least one WQA site
meta_join <- st_join(meta_snap, rivers, join=st_nearest_feature)
length(unique(meta_join$identifier)) #13758
wq_join <- st_join(wq_snap, rivers, join=st_nearest_feature)
length(unique(wq_join$identifier)) #11423
match_segs <- intersect(unique(wq_join$identifier), unique(meta_join$identifier))
length(match_segs) #7902

##Join stws to enable quick lookup of any with a match_seg
stws_join <- st_join(stws_snap, rivers, join=st_nearest_feature) #Use unique(stws_join$identifier) to lookup matches
length(unique(stws_join$identifier)) #13859

##Subset rivers and stws by match_segs
rivers <- rivers[rivers$identifier %in% match_segs,]
stws_join <- stws_join[stws_join$identifier %in% match_segs,]

##Split rivers at each stw and add a sub_id
split_rivers <- function(x){
  message(paste0(which(match_segs==x), " of ", length(match_segs)))
  line <- rivers[rivers$identifier==x, ]
  if(x %in% stws_join$identifier){
    line <- st_collection_extract(st_split(line, st_buffer(stws_join[stws_join$identifier==x,], 1)), "LINESTRING") #Split line at each discharge
  }
  line$sub_id <- 1:nrow(line) #Assign sub_id for each split segment
  line
}
nrow(rivers) #6481
rivers <- unique(do.call(rbind, lapply(match_segs, split_rivers)))
rivers$sub_id <- paste(rivers$identifier, rivers$sub_id, sep="_")
nrow(rivers) #33256

##Find matching subsegments
meta_join2 <- st_join(meta_join[meta_join$identifier %in% match_segs,], rivers[,"sub_id"], join=st_nearest_feature)
length(unique(meta_join2$sub_id)) #10269
wq_join2 <- st_join(wq_join[wq_join$identifier %in% match_segs,], rivers[,"sub_id"], join=st_nearest_feature)
length(unique(wq_join2$sub_id)) #9293
match_subsegs <- intersect(unique(wq_join2$sub_id), unique(meta_join2$sub_id))
length(match_subsegs) #7870

##Add subsegment id to meta and wq, filtering out observations with no matches
meta_join2$SITE_ID2 <- meta_join2$SITE_ID

invs_meta$SITE_ID2 <- paste0(invs_meta$SITE_ID, "_invs")
invs_meta_matched <- left_join(invs_meta, st_drop_geometry(meta_join2)[,c("SITE_ID2", "sub_id")])
invs_meta_matched <- invs_meta_matched[!is.na(invs_meta_matched$sub_id),]
nrow(invs_meta_matched) #144160
length(unique(invs_meta_matched$SITE_ID)) #13312

fish_meta$SITE_ID2 <- paste0(fish_meta$SITE_ID, "_fish")
fish_meta_matched <- left_join(fish_meta, st_drop_geometry(meta_join2)[,c("SITE_ID2", "sub_id")])
fish_meta_matched <- fish_meta_matched[!is.na(fish_meta_matched$sub_id),]
nrow(fish_meta_matched) #11957
length(unique(fish_meta_matched$SITE_ID)) #3213

colnames(wq_join2)[2] <- "sample.samplingPoint.notation"
wq_matched <- left_join(wq, unique(st_drop_geometry(wq_join2)[,c("sample.samplingPoint.notation", "sub_id")]))
wq_matched <- wq_matched[!is.na(wq_matched$sub_id),]
nrow(wq_matched) #~5.9M
length(unique(wq_matched$sample.samplingPoint.notation)) #11200

length(unique(invs_meta_matched$sub_id)) #9450
length(unique(fish_meta_matched$sub_id)) #2549
length(unique(wq_matched$sub_id)) #9293

invs_meta_matched <- invs_meta_matched[invs_meta_matched$sub_id %in% match_subsegs,]
fish_meta_matched <- fish_meta_matched[fish_meta_matched$sub_id %in% match_subsegs,]
wq_matched <- wq_matched[wq_matched$sub_id %in% match_subsegs,]

length(unique(invs_meta_matched$sub_id)) #7468 subsegments
length(unique(fish_meta_matched$sub_id)) #1907 subsegments
length(unique(wq_matched$sub_id)) #7870 subsegments
nrow(invs_meta_matched) #123858 samples
length(unique(invs_meta_matched$SITE_ID)) #10758 sites
nrow(fish_meta_matched) #9338 samples
length(unique(fish_meta_matched$SITE_ID)) #2481 sites
nrow(wq_matched) #~5.1M observations
length(unique(wq_matched$sample.samplingPoint.notation)) #9627 sites

#Visualise number of sites per subsegment
sub_id.freq <- do.call(rbind, lapply(unique(wq_matched$sub_id),
                                     function(x) data.frame(sub_id=x,
                                                            nsites=length(unique(wq_matched$sample.samplingPoint.notation[wq_matched$sub_id==x])),
                                                            ts_len=difftime(min(wq_matched$date[wq_matched$sub_id==x]), max(wq_matched$date[wq_matched$sub_id==x]), units="days"))))

pdf("Matched WQA site summary per subsegment.pdf", width=7.5, height=2.2)
grid.arrange(
  ggplot(sub_id.freq, aes(x=nsites)) +
    geom_histogram(binwidth=1, fill="lightgrey", colour="black") +
    labs(x="Number of matched WQA sites in reach", y="Frequency"),
  ggplot(sub_id.freq, aes(x= 1+(-ts_len/365.25))) +
    geom_histogram(binwidth=1, fill="lightgrey", colour="black") +
    labs(x="Length of matched WQA timeseries (years)", y="Frequency"),
  ncol=2)
dev.off()

#Visualise distribution of determinand coverage
get_det_coverage <- function(x){
  x_data <- wq_matched[wq_matched$sub_id==x,]
  x_data <- as.data.frame(table(x_data$determinand.name)/length(unique(substr(x_data$date, 1,7)))*100) #% of sampled months with values
  colnames(x_data) <- c("det", "coverage")
  x_data$coverage[x_data$coverage>100] <- 100
  x_data <- left_join(data_frame(det=dets), x_data)
  x_data$coverage[is.na(x_data$coverage)] <- 0
  data.frame(sub_id=x, x_data)
}
det_coverage <- do.call(rbind, lapply(unique(wq_matched$sub_id), get_det_coverage))
det_coverage$det <- factor(as.character(det_coverage$det), levels=dets)

pdf("Determinand coverage per subsegment.pdf", width=7.5, height=8)
ggplot(det_coverage, aes(x=coverage)) +
  geom_histogram(binwidth=10, fill="lightgrey", colour="black") +
  facet_wrap(~det, ncol=2) +
  labs(x="Determinand coverage by subsegment (% of sampled months)", y="Frequency")
dev.off()

write.csv(invs_meta_matched, file.path(PATH_PROCESSED, "invs_meta_matched.csv"))
write.csv(fish_meta_matched, file.path(PATH_PROCESSED, "fish_meta_matched.csv"))
write.csv(wq_matched, file.path(PATH_PROCESSED, "wq_matched.csv"))

