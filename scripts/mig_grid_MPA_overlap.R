## How many of high route density cells overlap with MPAs? ## -----------------

pacman::p_load(dplyr, sf, ggplot2, sp, stringr, rasterVis)

datatypes <- c("raw_filtered","interpolated")
# datatypes <- c("raw_filtered","interpolated", "sensitivity_analysis")

# for(y in seq_along(datatypes)){
datatype <- datatypes[y]
if(datatype == "interpolated"){
  migrid <- readRDS("data/analysis/mig_grid/interpolated_mig_grid_10x10_22.rds")
} else if(datatype == "raw_filtered"){
  migrid <- readRDS("data/analysis/mig_grid/raw_filtered_mig_grid_10x10_22.rds")
} else if(datatype == "sensitivity_analysis"){
  # folder <- "data/sensitivity_analysis/raw_filtered_gpsids/" # data for sensitivity analysis (only GPS data for 2019/2020 data)
}

# threshold of # of mig routes in a cell to consider it 'high density' 
thresh <- quantile(na.omit(values(migrid)), .95)

hidense <- migrid>thresh  

mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp")

cellover <- extract(migrid, mpas, weights=T, normalizeWeight=F)    
cellover <- as.data.frame(cellover[1][[1]])  

# percentage of high route-density cells which overlap MPAs
percover <- nrow(subset(cellover, cellover$value>=thresh)) / nrow(subset(cellover, cellover$value>=0)) * 100  

percover

## How many MPAs does each turtle pass through during migration? ##


datatype <- datatypes[y]
if(datatype == "interpolated"){
  folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods
} else if(datatype == "raw_filtered"){
  folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods
} else if(datatype == "sensitivity_analysis"){
  # folder <- "data/sensitivity_analysis/raw_filtered_gpsids/" # data for sensitivity analysis (only GPS data for 2019/2020 data)
}
## Data input ~~~~~~~~~~~~~~~~~~
# choose whether to include sub-6 satellite data or naw
# satfilt <- "satfilt4"
satfilt <- "satfilt6"

period <- "migration"

tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))

if("id" %in% colnames(tracks)){
  tracks <- formatFields(tracks,
                         fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date") 
}

## filter to only turtles leaving Bijagos ##
tracks <- tracks %>% filter(destination != "Bijagos")

## filter out points w/in Bijos ## 
tracks <- tracks %>% filter(Latitude > 11.6)

n <- n_distinct(tracks$ID) # one individuals drops out (182454)?

# project tracks to data-centered projection
TD <- projectTracks(tracks, projType = "azim", custom=T) # equal area azimuthal proj

## mpas
# mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp")
mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021.shp")

mpas_prj <- spTransform(mpas, CRSobj=crs(TD))

TD_mpa_u <- over(TD, mpas_prj)

TD$mpa_yn <- ifelse(is.na(TD_mpa_u$WDPAID), F, 
                    ifelse(TD_mpa_u$NAME=="Banc d'Arguin National Park", F, 
                           ifelse(TD_mpa_u$NAME=="Ilhas Formosa, Nago & Tchedia (Urok)", F, T)))

ind_summ <- TD@data %>% group_by(ID) %>% summarise(
  n_pnts_mpa = sum(mpa_yn),
  visit_mpa  = ifelse(n_pnts_mpa > 0, T, F)
)
## % of individuals which do not pass through an MPA on their way to the foraging grounds
sum(ind_summ$visit_mpa == F)/n_distinct(ind_summ$ID) * 100

mapview(subset(TD, TD$ID %in% c("182459", "205282", "205277", "205288", "60886", "60890", "60893", "60898")))
