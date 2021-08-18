#------------------------------------------------------------------------------#
## Calculate % of raw data points w/in PAs irrespective of error radius ## 

## raw UNfiltered
## Data input ~~~~~~~~~~~~~~~~~~
folder <- "data/analysis/" # repository w/ datasets split into periods
folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods

tracks <- readRDS(list.files(folder, pattern=fixed(period), full.names = T))

# project tracks to data-centered projection
TD <- projectTracks(tracks, projType = "azim", custom=T) # equal area azimuthal proj

### remove on-beach points for inter-nesting period data ### 
if(period == "internest"){
  poilao <- data.frame(label="Poilão", "Longitude" = -15.725273, "Latitude" = 10.870530)
  poilao <- st_as_sf(poilao, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")
  
  poilao_buff <- poilao %>% st_transform(crs = TD@proj4string) %>% st_buffer(dist = 700) # 700m buffer from island center
  
  idxs <- TD %>% st_as_sf() %>% st_within(poilao_buff)
  inside <- apply(idxs, 1, any)
  
  TD <- st_as_sf(TD)[which(inside==F),] %>% as_Spatial()
  
  # only keep IDs which have had enough data to estimate KDE for 
  KDE_ids <- names(readRDS(paste0("data/analysis/UDs/ind_areas/individual_UDs_h0.89_c1_internest.rds")))
  
  TD <- subset(TD, ID %in% KDE_ids)
  
} else if(period == "foraging"){
  # only keep IDs which have had enough data to estimate KDE for 
  KDE_ids <- names(KDE <- readRDS(paste0("data/analysis/UDs/ind_areas/individual_UDs_h2.08_c1_foraging.rds")))
}

TD_ll <- spTransform(TD, CRS("+init=epsg:4326")) # back to lat/lon

### estimate proportion of positions within and outside MPAs ##----------------
folder <- "data/geodata/WA_conservation area_shp/"
files <- str_subset(list.files("data/geodata/WA_conservation area_shp/", full.names = T),  pattern = fixed(".shp"))

## load MPA shapefiles ## 
shp_list <- lapply(files, function(x) shapefile(x) ) # list of shps

## bind all areas of conservation interest from bijagos to mauritania ##
cons_areas <- bind(shp_list)
# mapview(cons_areas)

# filter out biosphere reserve for overlaying points
mpas <- subset(cons_areas, NAME != "Bijagos Archipelago Biosphere Reserve")

mpas_prj <- spTransform(mpas, CRSobj=crs(TD))

## extract info of MPA which each point falls over (NA if none) ##
mpas_prj_u   <- rgeos::gUnaryUnion(mpas_prj)

TD$over_mpa  <- over(TD, mpas_prj_u)
TD$over_mpa <- ifelse(is.na(TD$over_mpa), F, TD$over_mpa)
# pnts falling w/in MPAs
# TD_inmpa <- subset(TD, TD$over_mpa == 1)
# dist to MPA border (for pnts inside)
pnt2pol_dists <- rgeos::gDistance(TD, 
                                 as(mpas_prj_u, "SpatialLines"), 
                                      byid = TRUE)
TD$edgedist <- unname(pnt2pol_dists[1,])
# set error radius distance from purported position
TD$error_band <- ifelse(TD$Argos.Location.Quality=="3", 250, 
                           ifelse(TD$Argos.Location.Quality=="2", 500,
                                  ifelse(TD$Argos.Location.Quality=="1", 1500,
                                         ifelse(TD$Argos.Location.Quality=="0", 3000, 100))))
# is point definitely locatin in MPA? (i.e. is entire error band located therein)
TD$in_mpa_werror <- ifelse(TD$over_mpa == T & (TD$error_band < TD$edgedist), T, F)

mapview(TD) + mapview(as(mpas_prj_u, "SpatialLines"))
