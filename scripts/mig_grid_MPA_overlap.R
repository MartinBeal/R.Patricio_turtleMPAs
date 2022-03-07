## How many of high route density cells overlap with MPAs? ## -----------------

pacman::p_load(dplyr, sf, ggplot2, sp, stringr, rasterVis)

# migrid <- readRDS("data/analysis/mig_grid/rawdata_dbbmm_mig_grid_10x10_UD99_18_perc.rds") ## using 99% UD
migrid <- readRDS("data/analysis/mig_grid/rawdata_dbbmm_mig_grid_10x10_UD95_18_perc.rds") ## using 95% UD
# migrid <- readRDS("data/analysis/mig_grid/rawdata_dbbmm_mig_grid_10x10_UD95_22_perc.rds") ## using 95% UD

# % of migration routes passing through a cell to consider it 'high density' 
thresh <- 25 #%
# thresh <- 20 #%

hidense <- migrid>thresh  
mapview(hidense)

hi_tf <- hidense
hi_tf[values(hi_tf)<=0] <- NA

## just MPAs ## ---------------------------------------------------
# mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp") #w/out goree
mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_Nov2021_dissolve.shp") #w/ goree

cellover <- extract(migrid, mpas, weights=T, normalizeWeight=F)    
cellover <- as.data.frame(cellover[1][[1]])  

range(na.omit(cellover$value)) # range of route density in MPAs 
(range(na.omit(cellover$value)) / 100) * 18 # range of turtles in MPAs

## get hi dense cells that fall w/in MPA polygons
cellover_hi <- extract(hi_tf, mpas, weights=T, normalizeWeight=F)    
cellover_hi <- as.data.frame(cellover_hi[1][[1]])  

## total number of hi density cells # 
n_cells_hi <- sum(na.omit(values(hi_tf) == TRUE))

## number of hi density cells W/IN MPAs #
n_mpa <- sum(na.omit(cellover_hi$value))

## percentage of hi density cells in MPAs ## 
perccover_mpa <- n_mpa / n_cells_hi * 100
perccover_mpa

## BABR included ## ---------------------------------------------------
babr <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/BABR_polygon.shp")

cons_areas <- aggregate(rbind(mpas[,-c(31:32)], babr))

## get hi dense cells that fall w/in MPA polygons
cellover_hi_ca <- extract(hi_tf, cons_areas, weights=T, normalizeWeight=F)    
cellover_hi_ca <- as.data.frame(cellover_hi_ca[1][[1]])  

## number of hi density cells W/IN MPAs #
n_babr <- sum(na.omit(cellover_hi_ca$value))

## percentage of hi density cells in MPAs ## 
perccover_babr <- n_babr / n_cells_hi * 100
perccover_babr


# -----------------------------------------------------------------------------
## How many MPAs does each turtle pass through during migration? ## -----------
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

# mapview(subset(TD, TD$ID %in% c("182459", "205282", "205277", "205288", "60886", "60890", "60893", "60898")))
