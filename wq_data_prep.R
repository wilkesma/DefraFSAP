source(here::here("paths.R"))

####Packages
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)
library(gridExtra)

##EA WQ archive data
wq <- lapply(seq(2002, 2023), function(x) read.csv(paste0(PATH_EXTERNAL, "/wq_archive/", x, ".csv"))) #WQA website down - can only use 2002-2023 for now
lapply(wq, colnames) #All sets of column names are the same
wq <- do.call(rbind, wq)

#Basic exploration and filtering
head(wq)
table(wq$sample.sampledMaterialType.label)
wq <- wq[wq$sample.sampledMaterialType.label=="RIVER / RUNNING SURFACE WATER",]
table(wq$sample.isComplianceSample)
table(wq$sample.purpose.label)

#Get determinand names used in ChemPop for populating a spreadsheet
detect.det <- function(x){
  table(wq$determinand.definition[wq$determinand.definition %in%
                                    unique(wq$determinand.definition[str_detect(wq$determinand.definition, x) |
                                                                       str_detect(wq$determinand.definition, str_to_upper(x)) |
                                                                       str_detect(wq$determinand.definition, str_to_title(x))])])
}

detect.det("pH")
detect.det("nitrate")
detect.det("nitrite")
detect.det("ammonia")
detect.det("phosphate")
detect.det("zinc")
detect.det("copper")
detect.det("temp") #Temperature of Water

#Wrangle WQ data
wq.vars <- data.frame(
  Variable = c(
    "Zn_d",
    "Cu_d",
    "NO3",
    "NO2",
    "pH",
    "PO4",
    "AmN",
    "wT"
  ),
  Variable_name = c(
    "Zinc, Dissolved",
    "Copper, Dissolved",
    "Nitrate as N",
    "Nitrite as N",
    "pH",
    "Orthophosphate, reactive as P",
    "Ammoniacal Nitrogen as N",
    "Temperature of Water"
  ),
  stringsAsFactors = TRUE
)
wq.vars$determinand.definition <- wq.vars$Variable.name #To ensure link to wq
wq <- wq[wq$determinand.definition %in% wq.vars$determinand.definition,]
wq <- left_join(wq, wq.vars[,c("determinand.definition", "Variable")])
table(wq$Variable)

colnames(wq)
wq <- wq[,c("sample.samplingPoint.notation", "sample.samplingPoint.label",
            "sample.sampleDateTime", "determinand.definition",
            "resultQualifier.notation", "result", "determinand.unit.label",
            "sample.isComplianceSample", "sample.purpose.label",
            "sample.samplingPoint.easting", "sample.samplingPoint.northing", "Variable")]
colnames(wq)[ncol(wq)] <- "determinand.name"

nrow(wq) #7,750,172

#Site coordinates checks
unique_sites <- unique(wq$sample.samplingPoint.notation)
check.coords <- function(x){
  message(paste0("Site ", which(unique_sites==x), " of ", length(unique_sites)))
  dists <- dist(wq[wq$sample.samplingPoint.notation==x,c("sample.samplingPoint.easting", "sample.samplingPoint.northing")])
  data.frame(site=x, mean_dist=mean(dists), max_dist=max(dists))
}
site_dists <- do.call(rbind, mclapply(unique(wq$sample.samplingPoint.notation), check.coords, mc.cores=50))
summary(site_dists)
site_dists[is.na(site_dists)] <- 0
site_dists$max_dist[is.infinite(site_dists$max_dist)] <- 0

pdf("Maximum WQA site movement.pdf", width=6.7/2, height=2)
ggplot(site_dists, aes(x=max_dist/1000)) +
  geom_histogram() +
  scale_x_log10(breaks=c(0.01, 1, 100), labels=c("0.01", "1", "100")) +
  labs(x="Maximum distance variation (km)", y="Count of sites")
dev.off()

site_dists[site_dists$max_dist>1000,] #Furthest realistic move is 1.1 km (SO-Y0004370)

#Deal with sites moving >2 km
sites_moved <- unique(site_dists$site[site_dists$max_dist>2000])
for(i in sites_moved){
  print(unique(wq[wq$sample.samplingPoint.notation==i,c("sample.samplingPoint.notation", "sample.samplingPoint.easting", "sample.samplingPoint.northing")]))
} #4/5 due to same bad coordinates for some samples (5000 x 4); 1/5 due to typo on some coords
#Just exclude them as they are only 5 of 19,767 sites
wq <- wq[-which(wq$sample.samplingPoint.notation %in% sites_moved),]

#For each remaining site, create new easting and northing columns with latest coordinates
wq$date <- as.POSIXct(wq$sample.sampleDateTime)

latest_locations <- wq %>%
  group_by(sample.samplingPoint.notation) %>%
  slice_max(order_by = desc(date), n = 1) %>%
  select(sample.samplingPoint.notation, 
         easting = sample.samplingPoint.easting, 
         northing = sample.samplingPoint.northing)

wq <- wq %>% left_join(unique(latest_locations), by = "sample.samplingPoint.notation")

unique(wq[wq$sample.samplingPoint.notation=="SO-Y0004370", c("sample.samplingPoint.notation", "sample.samplingPoint.easting", "sample.samplingPoint.northing", "easting", "northing")])

colnames(wq)
nrow(wq) #7,747,848

##Exclude sites on mine drains, adits etc., with mean pH <4
#First, show why it's necessary (it was causing issues with modelling WQ in subsegments with multiple WQ sites)
wq_mines <- wq[wq$determinand.name %in% c("pH", "Zn_d"),]
wq_mines$mine <- str_detect(wq_mines$sample.samplingPoint.label, "MINE ")
wq_mines <- aggregate(result~sample.samplingPoint.notation+determinand.name+mine, wq_mines, mean)
table(wq_mines$mine[wq_mines$determinand.name=="pH"])/nrow(wq_mines[wq_mines$determinand.name=="pH",])*100 #0.3% of sites with pH observations are associated with mines
wq_mines_wider <- pivot_wider(wq_mines[,c("sample.samplingPoint.notation", "result", "determinand.name", "mine")],
                              names_from=determinand.name, values_from=result)

pdf("Mine drainage pH and Zn_d mean site distributions.pdf", width=6.7, height=8)
grid.arrange(
  ggplot(wq_mines[wq_mines$determinand.name=="pH",], aes(x=result, fill=mine)) +
    geom_density(alpha=0.4, colour=NA) +
    geom_vline(xintercept=4, linetype=2) +
    scale_x_continuous(breaks=c(2,6,10,14)) +
    labs(x="pH", fill="Mine drainage", y="Probability density") +
    theme(legend.position="top"),
  ggplot(wq_mines[wq_mines$determinand.name=="Zn_d",], aes(x=result, fill=mine)) +
    geom_density(alpha=0.4, colour=NA) +
    scale_x_log10() +
    labs(x=expression("Dissolved zinc (" * mu * "g L"^{-1} * ")"), fill="Mine drainage", y="Probability density") +
    theme(legend.position="top"),
  ggplot(wq_mines_wider, aes(x=pH, y=Zn_d, colour=mine)) +
    geom_point(shape=16, alpha=0.4) +
    geom_point(data=wq_mines_wider[wq_mines_wider$mine==TRUE,], aes(x=pH, y=Zn_d), inherit.aes=FALSE, shape=16, colour="#00BFC4", alpha=0.1) +
    geom_vline(xintercept=4, linetype=2) +
    labs(x="pH", y=expression("Dissolved zinc (" * mu * "g L"^{-1} * ")"), fill="Mine drainage") +
    scale_y_log10() +
    theme(legend.position="none"),
  layout_matrix=rbind(c(1,2), c(3,3), c(3,3)))
dev.off()

#Now, exclude sites that are mine drains, adits, etc
length(unique(wq_mines$sample.samplingPoint.notation))
nrow(unique(wq_mines[,c("sample.samplingPoint.notation", "mine")])) #All sample.samplingPoint.notations have a single mine classification
wq_mines_sites <- unique(wq_mines$sample.samplingPoint.notation[wq_mines$mine==TRUE])
length(wq_mines_sites) #51 sites

wq <- wq[-which(wq$sample.samplingPoint.notation %in% wq_mines_sites),]

#Repeat plots to show difference
wq_mines <- wq[wq$determinand.name %in% c("pH", "Zn_d"),]
wq_mines$mine <- str_detect(wq_mines$sample.samplingPoint.label, "MINE ")
wq_mines <- aggregate(result~sample.samplingPoint.notation+determinand.name+mine, wq_mines, mean)
table(wq_mines$mine[wq_mines$determinand.name=="pH"])/nrow(wq_mines[wq_mines$determinand.name=="pH",])*100 #0.2% of sites with pH observations are associated with mines
wq_mines_wider <- pivot_wider(wq_mines[,c("sample.samplingPoint.notation", "result", "determinand.name", "mine")],
                              names_from=determinand.name, values_from=result)

pdf("Mine drainage pH and Zn_d mean site distributions AFTER EXCLUSION.pdf", width=6.7, height=8)
grid.arrange(
  ggplot(wq_mines[wq_mines$determinand.name=="pH",], aes(x=result, fill=mine)) +
    geom_density(alpha=0.4, colour=NA) +
    geom_vline(xintercept=4, linetype=2) +
    scale_x_continuous(breaks=c(2,6,10,14)) +
    labs(x="pH", fill="Mine drainage", y="Probability density") +
    theme(legend.position="top"),
  ggplot(wq_mines[wq_mines$determinand.name=="Zn_d",], aes(x=result, fill=mine)) +
    geom_density(alpha=0.4, colour=NA) +
    scale_x_log10() +
    labs(x=expression("Dissolved zinc (" * mu * "g L"^{-1} * ")"), fill="Mine drainage", y="Probability density") +
    theme(legend.position="top"),
  ggplot(wq_mines_wider, aes(x=pH, y=Zn_d, colour=mine)) +
    geom_point(shape=16, alpha=0.4) +
    geom_point(data=wq_mines_wider[wq_mines_wider$mine==TRUE,], aes(x=pH, y=Zn_d), inherit.aes=FALSE, shape=16, colour="#00BFC4", alpha=0.1) +
    geom_vline(xintercept=4, linetype=2) +
    labs(x="pH", y=expression("Dissolved zinc (" * mu * "g L"^{-1} * ")"), fill="Mine drainage") +
    scale_y_log10() +
    theme(legend.position="none"),
  layout_matrix=rbind(c(1,2), c(3,3), c(3,3)))
dev.off()

nrow(wq) #7,742,658

#Save data
write.csv(wq, file.path(PATH_PROCESSED, "ea_wq.csv"), row.names=FALSE)
