source(here::here("paths.R"))

####Packages
library(sf)
library(rnrfa)
library(ggplot2)
library(dplyr)
library(taxize)
library(reshape2)
library(stringr)
library(gridExtra)

####Functions
get.season <- function(x){
  x.month <- as.numeric(format(x, "%m"))
  if(x.month==12 | x.month==1 | x.month==2) { season <- "winter" }
  if(x.month==3 | x.month==4 | x.month==5) { season <- "spring" }
  if(x.month==6 | x.month==7 | x.month==8) { season <- "summer" }
  if(x.month==9 | x.month==10 | x.month==11) { season <- "autumn" }
  season
}

get.coords <- function(x){
  if("NGR_10_FIG" %in% colnames(x)){
    osg_parse(x$NGR_10_FIG)
  } else{
    osg_parse(x$SURVEY_RANKED_NGR)
  }
}

get.basin <- function(meta){
  meta.sp <- st_as_sf(meta, coords = c("easting", "northing"),
                      crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
  meta.sp$basin <- NA
  meta.basins <- as.data.frame(st_intersects(meta.sp, basins, sparse=F))
  for(i in 1:length(meta.sp$basin)){
    message(paste0(i, " of ", length(meta.sp$basin)))
    if(sum(as.numeric(meta.basins[i,]))>0){
      meta.sp$basin[i] <- unique(basins$rbd_name)[which.max(meta.basins[i,])]
    } else{ meta.sp$basin[i] <- NA }
  }
  meta.sp$basin
}

get.taxares <- function(year, basin, abun, meta, x.info){
  message(paste(year, basin))
  y.abun <- abun[as.character(meta$ANALYSIS_ID[meta$year==year & meta$basin==basin]),]
  if(!is.null(nrow(y.abun))){
    n <- nrow(y.abun)
    y.abun <- data.frame(TAXON_LIST_ITEM_KEY=colnames(y.abun), freq=as.numeric(colSums(y.abun)/sum(y.abun))*100)
    y.abun <- left_join(y.abun, x.info[,c("TAXON_LIST_ITEM_KEY", "TAXON_RANK")])
    if(any(!is.na(y.abun$freq))){
      y.abun <- aggregate(freq~TAXON_RANK, y.abun, sum)
      y.abun$TAXON_RANK <- as.character(y.abun$TAXON_RANK)
      output <- data.frame(year=year, basin=basin, n=n,  y.abun)
    } else{
      output <- data.frame(year=year, basin=basin, n=0, freq=NA, TAXON_RANK=NA)
    }
  } else{
    output <- data.frame(year=year, basin=basin, n=1, freq=NA, TAXON_RANK=NA) #Returns NA if <2 samples
  }
  output
}

get.taxon.info <- function(x){
  message(x$name[nrow(x)])
  x.class <- x[x$rank %in% c("phylum", "class", "order", "family", "genus", "species", "species group", "species hybrid"), c("name", "rank", "id")]
  x.df <- data.frame(final.name=NA, final.rank=NA, phylum=NA, class=NA, order=NA, family=NA, genus=NA, species=NA)
  x.df$final.name <- x.class$name[nrow(x.class)]
  x.df$final.rank <- x.class$rank[nrow(x.class)]
  if(substr(x.df$final.rank, 1, 8)=="species "){
    x.df$final.rank <- "species"
  }
  if("phylum" %in% x.class$rank){ x.df$phylum <- x.class$name[x.class$rank=="phylum"] }
  if("class" %in% x.class$rank){ x.df$class <- x.class$name[x.class$rank=="class"] }
  if("order" %in% x.class$rank){ x.df$order <- x.class$name[x.class$rank=="order"] }
  if("family" %in% x.class$rank){ x.df$family <- x.class$name[x.class$rank=="family"] }
  if("genus" %in% x.class$rank){ x.df$genus <- x.class$name[x.class$rank=="genus"] }
  if("species" %in% substr(x.class$rank, 1, 7)){ x.df$species <- x.class$name[substr(x.class$rank, 1, 7)=="species"] }
  x.df
}

merge.cols <- function(x, df){
  x.taxa <- as.data.frame(df[,which(colnames(df)==x)])
  rowSums(x.taxa)
}

####Spatial data
basins <- read_sf(PATH_WFD_BASINS)
basins <- st_transform(basins,
                       crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")

####Raw data
##Load EA datasets (includes other BQEs)
england.bio.data <- list(
  INV_OPEN_DATA_METRICS.rds = readRDS(PATH_EA_INV_METRICS),
  INV_OPEN_DATA_SITE.rds    = readRDS(PATH_EA_INV_SITE),
  INV_OPEN_DATA_TAXA.rds    = readRDS(PATH_EA_INV_TAXA),
  OPEN_DATA_TAXON_INFO.rds  = readRDS(PATH_EA_TAXON_INFO),
  FW_Fish_Counts.rds        = readRDS(PATH_EA_FISH_COUNTS)
)

####Invertebrates
##invs_metadata
invs_meta <- merge(england.bio.data$INV_OPEN_DATA_METRICS.rds, england.bio.data$INV_OPEN_DATA_SITE.rds)
invs_meta$easting <- get.coords(invs_meta)$easting
invs_meta$northing <- get.coords(invs_meta)$northing #Duplicates existing columns but better to be consistent
invs_meta$basin <- get.basin(invs_meta)
invs_meta$year <- format(as.Date(as.character(invs_meta[,"SAMPLE_DATE"]), format = "%d/%m/%Y"), format="%Y")
invs_meta$time <- as.numeric(difftime(as.Date(as.character(invs_meta[,"SAMPLE_DATE"]), format = "%d/%m/%Y"), "1990-01-01", units = c("days")))/365.25
invs_meta$julian <- as.POSIXlt(invs_meta[,"SAMPLE_DATE"], format = "%d/%m/%Y")$yday
invs_meta$season <- sapply(as.POSIXlt(invs_meta[,"SAMPLE_DATE"], format = "%d/%m/%Y"), get.season)
invs_meta$julian_season <- NA
invs_meta$julian_season[invs_meta$season=="spring"] <- invs_meta$julian[invs_meta$season=="spring"]-aggregate(julian~season, invs_meta, min)[2,2]
invs_meta$julian_season[invs_meta$season=="autumn"] <- invs_meta$julian[invs_meta$season=="autumn"]-aggregate(julian~season, invs_meta, min)[1,2]
nrow(invs_meta) #276080 samples

##Data filtering and checking
#Metadata
pdf("Invert samples per year 1965-2025.pdf", height=7, width=6.7/2)
ggplot(data.frame(Year=as.numeric(names(table(invs_meta$year))), Samples=as.numeric(table(invs_meta$year))), aes(x=Year, y=Samples))+
  geom_bar(stat="identity") + coord_flip() + scale_y_log10() + scale_x_reverse(breaks=seq(1965, 2025, by=1)) +
  theme(axis.title.y=element_blank())
dev.off()
invs_meta <- invs_meta[invs_meta$year>=1985,] #Earliest year worth considering (all years consistently have >1000 samples from 1985)
invs_meta <- invs_meta[invs_meta$year<2025,] #2025 only partial data

table(invs_meta$WATERBODY_TYPE)/nrow(invs_meta)*100
invs_meta <- invs_meta[invs_meta$WATERBODY_TYPE=="WBRV",] #Rivers only (96.98% of remaining samples)

table(invs_meta$ANALYSIS_METHOD)/nrow(invs_meta)*100
invs_meta <- invs_meta[which(invs_meta$ANALYSIS_METHOD %in% c("ANAA", "ANLA", "ANLE")),] #Only main abundance estimate methods (>99% of remaining samples)

table(invs_meta$SAMPLE_METHOD)/nrow(invs_meta)*100
invs_meta <- invs_meta[which(invs_meta$SAMPLE_METHOD=="S3PO"),] #Just standard 3-minute kick samples (99% of remaining records)

table(invs_meta$season)/nrow(invs_meta)*100
table(invs_meta$season[as.numeric(invs_meta$year)>2001])/nrow(invs_meta[as.numeric(invs_meta$year)>2001,])*100
invs_meta <- invs_meta[which(invs_meta$season %in% c("spring", "autumn")),] #Majority of records in these seasons (79% of remaining records), and the vast majority of records (93%) since 2002

pdf("Invert samples per year and basin 1985-2024.pdf", height=10.2, width=6.7)
ggplot(melt(table(invs_meta[,c("year", "basin")])), aes(x=year, y=value))+
  geom_bar(stat="identity") + coord_flip() + scale_y_log10() + scale_x_reverse(breaks=seq(1985, 2025, by=5)) +
  facet_wrap(~basin, nrow=3) + ylab("Samples") +
  theme(axis.title.y=element_blank())
dev.off()

invs_meta <- invs_meta[invs_meta$year>=1990,] #Earliest year when all basins have data

nrow(invs_meta) #194261 samples remaining

pdf("Invert samples per year and basin 1990-2024.pdf", height=10.2, width=6.7)
ggplot(melt(table(invs_meta[,c("year", "basin")])), aes(x=year, y=value))+
  geom_bar(stat="identity") + coord_flip() + scale_y_log10() + scale_x_reverse(breaks=seq(1985, 2025, by=5)) +
  facet_wrap(~basin, nrow=3) + ylab("Samples") +
  theme(axis.title.y=element_blank())
dev.off()

#Checking consistency of coordinates for SITE_IDs
unique.coords <- do.call(c, lapply(invs_meta$SITE_ID, function(x) nrow(unique(invs_meta[invs_meta$SITE_ID==x, c("FULL_EASTING", "FULL_NORTHING")]))))
summary(unique.coords) #All sites have only a single set of coordinates
rm(unique.coords)

##Abundance data
gc()
invs_abun <- with(england.bio.data$INV_OPEN_DATA_TAXA.rds[which(england.bio.data$INV_OPEN_DATA_TAXA.rds$ANALYSIS_ID %in% invs_meta$ANALYSIS_ID),], tapply(TOTAL_ABUNDANCE, list(ANALYSIS_ID, TAXON_LIST_ITEM_KEY), FUN=sum))
gc()
invs_abun[is.na(invs_abun)] <- 0
invs_abun <- invs_abun[,colSums(invs_abun)>0]
invs_abun <- invs_abun[as.character(invs_meta$ANALYSIS_ID),]
nrow(invs_abun); ncol(invs_abun) #194261 samples, 2315 taxa

##Taxonomic classifications
info <- england.bio.data$OPEN_DATA_TAXON_INFO.rds
invs_info <- info[which(info$TAXON_LIST_ITEM_KEY %in% colnames(invs_abun)),]
row.names(invs_info) <- invs_info$TAXON_LIST_ITEM_KEY
invs_info <- invs_info[colnames(invs_abun),]
nrow(invs_info) #2315

#Basic checking and filtering
table(invs_info$TAXON_TYPE)/nrow(info)
invs_info <- invs_info[invs_info$TAXON_TYPE=="Other Macroinvertebrates",] #89.5% are this type; others are mostly fish - presumably other groups recorded in kicknet samples
nrow(invs_info) #2073 taxa remaining

invs_abun <- invs_abun[,as.character(info$TAXON_LIST_ITEM_KEY)]
nrow(invs_abun); ncol(invs_abun) #194261 samples, 2073 taxa

#Taxonomic resolution over time (assuming info$TAXON_RANK is correct)
taxares <- do.call(rbind, lapply(seq(1990, 2024, by=1), function(year) do.call(rbind, lapply(unique(invs_meta$basin), function(basin) get.taxares(year, basin, invs_abun, invs_meta)))))
taxares$TAXON_RANK[which(taxares$TAXON_RANK %in% c("Species sensu lato", "Subspecies", "Form"))] <- "Species"
taxares <- aggregate(freq~year+basin+TAXON_RANK, taxares, sum)
taxares$TAXON_RANK <- factor(taxares$TAXON_RANK, levels= c("Phylum", "Class", "Order", "Infraorder", "Family", "Subfamily", "Genus", "Subgenus", "Tribe", "Species group", "Species aggregate", "Species"))

pdf("Invs taxonomic resolution over time.pdf", height=5, width=6.7)
ggplot(taxares, aes(x=year, y=freq, fill=TAXON_RANK)) +
  geom_area(alpha=0.65) + facet_wrap(~basin) +
  xlim(c(1990,2024)) +
  ylab("Detected invs_abundance (%)") +
  scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#fb9a99','#e31a1c')) +
  theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#Subset info to schedule 2 species
sc2 <- read.csv(PATH_SCHEDULE2)
sc2$level <- "species"
sc2$level[which(str_detect(sc2$species, " spp."))] <- "genus"
sc2$species <- str_replace(sc2$species, " spp.", "")
colnames(sc2)[1] <- "taxon"

invs_sc2 <- sc2[sc2$group=="Freshwater invertebrates", c("taxon", "level")]

#Deal with synonyms/groups/spp. aggs
invs_sc2$taxon.alt1 <- NA
invs_sc2$taxon.alt2 <- NA
invs_sc2$taxon.alt3 <- NA
invs_sc2$taxon.alt4 <- NA

invs_sc2$taxon.alt1[invs_sc2$taxon=="Alboglossiphonia heteroclita"] <- "Glossiphonia heteroclita"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Ampullaceana balthica"] <- "Radix balthica"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Baetis muticus"] <- "Alainites muticus"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Baetis scambus group"] <- "Baetis scambus"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Baetis scambus group"] <- "Baetis scambus/fuscatus"
invs_sc2$taxon.alt3[invs_sc2$taxon=="Baetis scambus group"] <- "Baetis fuscatus"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Baetis vernus"] <- "Baetis tenax"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Bithynia leachii"] <- "Bithynia leachi"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Caenis luctuosa group"] <- "Caenis luctuosa"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Caenis luctuosa group"] <- "Caenis luctuosa/macrura"
invs_sc2$taxon.alt3[invs_sc2$taxon=="Caenis luctuosa group"] <- "Caenis macrura"
invs_sc2$taxon.alt4[invs_sc2$taxon=="Caenis luctuosa group"] <- "Caenis moesta"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Schmidtea polychroa group"] <- "Schmidtea polychroa"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Schmidtea polychroa group"] <- "Schmidtea lugubris or polychroa"
invs_sc2$taxon.alt3[invs_sc2$taxon=="Schmidtea polychroa group"] <- "Schmidtea lugubris"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Gammarus pulex group"] <- "Gammarus pulex"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Gammarus pulex group"] <- "Gammarus pulex/fossarum agg."
invs_sc2$taxon.alt3[invs_sc2$taxon=="Gammarus pulex group"] <- "Gammarus fossarum"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Haliplus ruficollis group"] <- "Haliplus ruficollis"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Lepidostoma basale"] <- "Lasiocephala basalis"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Nebrioporus depressus group"] <- "Nebrioporus depressus"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Nebrioporus depressus group"] <- "Nebrioporus elegans"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Nemurella pictetii"] <- "Nemurella picteti"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Perlodes mortoni"] <- "Perlodes microcephala"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Simulium angustitarse group"] <- "Simulium angustitarse"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Simulium argyreatum group"] <- "Simulium argyreatum"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Simulium argyreatum group"] <- "Simulium argyreatum/variegatum"
invs_sc2$taxon.alt3[invs_sc2$taxon=="Simulium argyreatum group"] <- "Simulium variegatum"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Simulium aureum group"] <- "Simulium aureum"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Simulium cryophilum-vernum group"] <- "Simulium cryophilum"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Simulium cryophilum-vernum group"] <- "Simulium vernum"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Simulium ornatum group"] <- "Simulium ornatum"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Simulium ornatum group"] <- "Simulium ornatum/intermedium/trifasciatum"
invs_sc2$taxon.alt3[invs_sc2$taxon=="Simulium ornatum group"] <- "Simulium intermedium"
invs_sc2$taxon.alt4[invs_sc2$taxon=="Simulium ornatum group"] <- "Simulium trifasciatum"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Simulium tuberosum complex"] <- "Simulium tuberosum"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Siphonoperla torrentium"] <- "Chloroperla torrentium"
invs_sc2$taxon.alt1[invs_sc2$taxon=="Rhyacophila fasciata"] <- "Rhyacophila dorsalis/septentrionis"
invs_sc2$taxon.alt2[invs_sc2$taxon=="Rhyacophila fasciata"] <- "Rhyacophila septentrionis"
invs_sc2$taxon.alt3[invs_sc2$taxon=="Rhyacophila fasciata"] <- "Rhyacophila septentrionis/fasciata"

info <- info[-which(info$PREFERRED_TAXON_NAME=="Musculium transversum"),] #Was causing trouble later - remove as the target shouldn't include this single INNS
info$PARENT_TAXON_NAME[which(info$TAXON_NAME=="Baetis niger")] <- "Nigrobaetis" #To keep consistent with Schedule 2 taxonomy

link.invs_sc2 <- function(i){
  taxon <- do.call(c, (invs_sc2[i,-which(colnames(invs_sc2)=="level")]))
  taxon <- taxon[!is.na(taxon)]
  message(taxon[1])
  level <- invs_sc2$level[i]
  if(level=="species"){
    output <- data.frame(taxon=taxon[1], level=level, invs_info[invs_info$PREFERRED_TAXON_NAME %in% taxon,])
  } else{
    output <- data.frame(taxon=taxon, level=level, rbind(invs_info[invs_info$PREFERRED_TAXON_NAME %in% taxon,], invs_info[invs_info$PARENT_TAXON_NAME %in% taxon,]))
  }
  output
}
invs_info.sc2 <- do.call(rbind, lapply(1:nrow(invs_sc2), link.invs_sc2))

##Filter invs_abun df, rename columns to species names, aggregate synonymous species
invs_abun.sc2 <- invs_abun[,as.character(invs_info.sc2$TAXON_LIST_ITEM_KEY)]
ncol(invs_abun.sc2) #555

colnames(invs_abun.sc2) <- invs_info.sc2$taxon
invs_abun.sc2 <- do.call(cbind, lapply(colnames(invs_abun.sc2), function(x) merge.cols(x, invs_abun.sc2)))
colnames(invs_abun.sc2) <- invs_info.sc2$taxon
invs_abun.sc2 <- invs_abun.sc2[,!duplicated(colnames(invs_abun.sc2))] 
ncol(invs_abun.sc2) #235 non-duplicated species
summary(colSums(invs_abun.sc2))

#Order and check consistency of data frames
invs_abun.sc2 <- invs_abun.sc2[,order(colnames(invs_abun.sc2))]
invs_abun.sc2 <- invs_abun.sc2[order(as.numeric(row.names(invs_abun.sc2))),]
invs_info.sc2 <- invs_info.sc2[order(invs_info.sc2$taxon),]
invs_meta <- invs_meta[which(as.character(invs_meta$ANALYSIS_ID) %in% row.names(invs_abun.sc2)), ]
invs_meta <- invs_meta[order(as.numeric(as.character(invs_meta$ANALYSIS_ID))),]

nrow(invs_abun.sc2); nrow(invs_meta)
ncol(invs_abun.sc2); nrow(invs_info.sc2) #Doesn't need to match
all(row.names(invs_abun.sc2)==as.character(invs_meta$ANALYSIS_ID))

##Saving
save(invs_abun.sc2, invs_info.sc2, invs_meta, file=file.path(PATH_PROCESSED, "invs_data.RData"))

####Fish
##Raw abundance data
fish_abun <- with(england.bio.data$FW_Fish_Counts.rds, tapply(SPCSNO, list(SURVEY_ID, LATIN_NAME), FUN=sum)) #Population estimate (Carle and Strub MWL) for multi-pass strategies only
fish_abun <- fish_abun[,-1] #The first column doesn't have a name
nrow(fish_abun); ncol(fish_abun) #75444 samples, 90 taxa

#Metadata, including getting basic additional data
fish_meta <- england.bio.data$FW_Fish_Counts.rds[!duplicated(england.bio.data$FW_Fish_Counts.rds$SURVEY_ID), c(1:27, 47)]
fish_meta$easting <- get.coords(fish_meta)$easting
fish_meta$northing <- get.coords(fish_meta)$northing #Duplicates existing columns but better to be consistent
fish_meta$basin <- get.basin(fish_meta)
fish_meta[sample(1:nrow(fish_meta), 40),c("REGION", "basin")] #Check basin assignment - good
fish_meta$year <- format(as.Date(as.character(fish_meta[,"EVENT_DATE"]), format = "%d/%m/%Y"), format="%Y")
fish_meta$time <- as.numeric(difftime(as.Date(as.character(fish_meta[,"EVENT_DATE"]), format = "%d/%m/%Y"), "2002-01-01", units = c("days")))/365.25
fish_meta$julian <- as.POSIXlt(fish_meta[,"EVENT_DATE"], format = "%d/%m/%Y")$yday
fish_meta$season <- sapply(as.POSIXlt(fish_meta[,"EVENT_DATE"], format = "%d/%m/%Y"), get.season)
nrow(fish_meta) #75444 samples

##Data filtering and checking
#Metadata
table(fish_meta$year)

pdf("Fish samples per year 1973-2025.pdf", height=7, width=6.7/2)
ggplot(data.frame(Year=as.numeric(names(table(fish_meta$year))), Samples=as.numeric(table(fish_meta$year))), aes(x=Year, y=Samples))+
  geom_bar(stat="identity") + coord_flip() + scale_y_log10() + scale_x_reverse(breaks=seq(1970, 2025, by=1)) +
  theme(axis.title.y=element_blank())
dev.off()
fish_meta <- fish_meta[fish_meta$year>=1990,] #Earliest year worth considering (all years consistently have close to 1000 samples from 1990, and this aligns with the invertebrate data)
fish_meta <- fish_meta[fish_meta$year<2025,] #2025 only partial data

fish_meta <- fish_meta[fish_meta$SURVEY_METHOD %in% c("AC ELECTRIC FISHING", "DC ELECTRIC FISHING", "ELECTRIC FISHING", "PDC ELECTRIC FISHING"),]
fish_meta <- fish_meta[fish_meta$SURVEY_STRATEGY %in% c("CATCH DEPLETION SAMPLE", "CATCH DEPLETION SAMPLE (PART WIDTH)"),] #Only catch depletion surveys
fish_meta <- fish_meta[!is.na(fish_meta$FISHED_AREA), ]
fish_meta <- fish_meta[fish_meta$FISHED_AREA>0, ]

table(fish_meta$season)/nrow(fish_meta)*100 #   autumn, spring and summer roughly equal proportions; winter only 3.7% of samples
fish_meta <- fish_meta[which(fish_meta$season %in% c("spring", "summer", "autumn")),] #Removing winter in case of bias

pdf("Fish samples per year and basin 1990-2024.pdf", height=10.2, width=6.7)
ggplot(melt(table(fish_meta[,c("year", "basin")])), aes(x=year, y=value))+
  geom_bar(stat="identity") + coord_flip() + scale_y_log10() + scale_x_reverse(breaks=seq(1990, 2025, by=5)) +
  facet_wrap(~basin, nrow=3) + ylab("Samples") +
  theme(axis.title.y=element_blank())
dev.off()

nrow(fish_meta) #194261 samples remaining

#Checking consistency of coordinates for SITE_IDs
unique.coords <- do.call(c, lapply(fish_meta$SITE_ID, function(x) nrow(unique(fish_meta[fish_meta$SITE_ID==x, c("easting", "northing")]))))
summary(unique.coords) #Many sites have multiple coordinates
problem_sites <- unique(fish_meta$SITE_ID[unique.coords>1])
length(problem_sites)/length(unique(fish_meta$SITE_ID)) #Affects 29% of sites

check_site <- function(x){
  x.obs <- st_as_sf(fish_meta[fish_meta$SITE_ID==x,], coords=c("easting", "northing"), crs=27700)
  max(as.numeric(st_distance(x.obs))) #Maximum distance between observations in m
}
site_dists <- sapply(problem_sites, check_site)
summary(site_dists)

pdf("Fish site distance movement.pdf", width=6, height=4)
ggplot(data.frame(dist=site_dists[site_dists<10000]), aes(x=dist)) +
  geom_histogram() +
  labs(x="Maximum distance between observations with the same SITE_ID (m)", y="Count of sites")
dev.off()

#1 km is a pragmatic choice - remove any sites with maximum distance >1 km, then assign the latest set of coordinates to every observation from each site
problem_sites <- problem_sites[site_dists>1000]
fish_meta <- fish_meta[-which(fish_meta$SITE_ID %in% problem_sites),]
nrow(fish_meta) #21260 of 22237 samples remaining

for(i in unique(fish_meta$SITE_ID)){
  message(paste0(which(unique(fish_meta$SITE_ID)==i), " of ", length(unique(fish_meta$SITE_ID))))
  fish_meta[fish_meta$SITE_ID==i,c("easting", "northing")] <- fish_meta[fish_meta$SITE_ID==i,c("easting", "northing")][nrow(fish_meta[fish_meta$SITE_ID==i,c("easting", "northing")]),]
}

unique.coords <- do.call(c, lapply(fish_meta$SITE_ID, function(x) nrow(unique(fish_meta[fish_meta$SITE_ID==x, c("easting", "northing")]))))
summary(unique.coords) #All observations within sites now have the same coordinates

#Subset taxa data based on meta subsetting
fish_abun <- fish_abun[as.character(fish_meta$SURVEY_ID),] #Only catch depletion samples
fish_abun[is.na(fish_abun)] <- 0
fish_abun <- fish_abun[,colSums(fish_abun)>0]
nrow(fish_abun); ncol(fish_abun) #21260 samples, 67 taxa

#Remove hybrids (Wlkes etal. 2025, Ecography, showed this is appropriate)
fish_abun <- fish_abun[,which(str_detect(colnames(fish_abun), " x ")==FALSE)] #Hybrids removed
ncol(fish_abun) #52 taxa

#Deal with synonyms
colnames(fish_abun)[which(colnames(fish_abun)=="Abramis bjoerkna")] <- "Blicca bjoerkna"
colnames(fish_abun)[which(colnames(fish_abun)=="Leuciscus cephalus")] <- "Squalius cephalus"

#Subset to schedule 2 species
fish_abun <- fish_abun[,which(colnames(fish_abun) %in% sc2$taxon[sc2$group=="Fish"])]

#Order and check consistency of data frames
fish_abun <- fish_abun[,order(colnames(fish_abun))]
fish_abun <- fish_abun[order(as.numeric(row.names(fish_abun))),]
fish_meta <- fish_meta[which(as.character(fish_meta$SURVEY_ID) %in% row.names(fish_abun)), ]
fish_meta <- fish_meta[order(as.numeric(as.character(fish_meta$SURVEY_ID))),]

nrow(fish_abun); nrow(fish_meta)
all(row.names(fish_abun)==as.character(fish_meta$SURVEY_ID))

##Saving
save(fish_abun, fish_meta, file=file.path(PATH_PROCESSED, "fish_data.RData"))

