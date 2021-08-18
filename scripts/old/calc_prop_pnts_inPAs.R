### estimate proportion of positions within and outside MPAs ##----------------

pacman::p_load(track2KBA, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter, data.table)

datatypes <- c("raw_filtered","interpolated")
# datatypes <- c("raw_filtered","interpolated", "sensitivity_analysis")

for(y in seq_along(datatypes)){
  datatype <- datatypes[y]
  if(datatype == "interpolated"){
    folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods
  } else if(datatype == "raw_filtered"){
    folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods
  } else if(datatype == "sensitivity_analysis"){
    # folder <- "data/sensitivity_analysis/raw_filtered_gpsids/" # data for sensitivity analysis (only GPS data for 2019/2020 data)
  }
}
## Data input ~~~~~~~~~~~~~~~~~~
# folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods
# folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods
# folder <- "data/sensitivity_analysis/raw_filtered_gpsids/" # data for sensitivity analysis (only GPS data for 2019/2020 data)

# choose whether to include sub-6 satellite data or naw
satfilt <- "satfilt6"

periods <- c("internest", "foraging")

TD_ll_list <- list()

for(i in seq_along(periods)){
  print(i)
  ## analyze inter-nesting, foraging or migration?
  period <- periods[i]
  # period <- "internest"
  # period <- "foraging"
  
  tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
  
  if("id" %in% colnames(tracks)){
  tracks <- formatFields(tracks,
                         fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date") 
  }
  
  # project tracks to data-centered projection
  TD <- projectTracks(tracks, projType = "azim", custom=T) # equal area azimuthal proj
  
  ### remove on-beach points for inter-nesting period data ### 
  if(period == "internest"){
    TD$period <- rep(period)
    poilao <- data.frame(label="Poilão", "Longitude" = -15.725273, "Latitude" = 10.870530)
    poilao <- st_as_sf(poilao, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")
    
    poilao_buff <- poilao %>% st_transform(crs = TD@proj4string) %>% st_buffer(dist = 700) # 700m buffer from island center
    
    idxs <- TD %>% st_as_sf() %>% st_within(poilao_buff)
    inside <- apply(idxs, 1, any)
    
    TD <- st_as_sf(TD)[which(inside==F),] %>% as_Spatial()
    
    # only keep IDs which have had enough data to estimate KDE for 
    # KDE_ids <- names(readRDS(paste0("data/analysis/UDs/ind_areas/individual_UDs_h0.89_c1_internest.rds")))
    
    TD <- subset(TD, ID %in% KDE_ids)
  
  } else if(period == "foraging"){
    TD$period <- rep(period)
    # only keep IDs which have had enough data to estimate KDE for 
    # KDE_ids <- names(KDE <- readRDS(paste0("data/analysis/UDs/ind_areas/individual_UDs_h2.08_c1_foraging.rds")))
  }
  
  TD_ll <- spTransform(TD, CRS("+init=epsg:4326")) # back to lat/lon
  
  ### estimate proportion of positions within and outside MPAs ##----------------
  # shpfolder <- "data/geodata/WA_conservation area_shp/"
  files <- str_subset(list.files("data/geodata/WA_conservation area_shp/", full.names = T),  pattern = fixed(".shp"))
  
  ## load MPA shapefiles ## 
  shp_list <- lapply(files, function(x) shapefile(x) ) # list of shps
  
  ## **** NEED TO UPDATE THESE **** bind all areas of conservation interest from bijagos to mauritania ##
  cons_areas <- bind(shp_list)
  # mapview(cons_areas)
  
  # filter out biosphere reserve for overlaying points
  mpas <- subset(cons_areas, NAME != "Bijagos Archipelago Biosphere Reserve")
  # mapview(mpas)
  
  ## extract info of MPA which each point falls over (NA if none) ##
  over_mpas <- over(TD_ll, mpas)
  
  # add info on which PA each point overlaps, and if none call is 'no PA'
  TD_ll$WDPA_PID <- ifelse(is.na(over_mpas$WDPA_PID), "no PA", over_mpas$WDPA_PID)
  TD_ll$PA_NAME <- ifelse(is.na(over_mpas$NAME), "no PA", over_mpas$NAME)
  # table(TD_ll$PA_NAME)  # how many points fall within EACH PA?
  TD_ll$over_PA <- ifelse(TD_ll$WDPA_PID == "no PA", F, T)
  table(TD_ll$over_PA)    # how many points fall in ANY PA?
  
  ## including BABR ##
  over_cons <- over(TD_ll, cons_areas)

  # add info on which CA each point overlaps, and if none call is 'no PA'
  TD_ll$WDPA_PID_CA <- ifelse(is.na(over_cons$WDPA_PID), "no CA", over_cons$WDPA_PID)
  TD_ll$CA_NAME <- ifelse(is.na(over_cons$NAME), "no CA", over_cons$NAME)
  table(TD_ll$CA_NAME)  # how many points fall within EACH Conservation Area?
  TD_ll$over_CA <- ifelse(TD_ll$WDPA_PID_CA == "no CA", F, T)
  table(TD_ll$over_CA)    # how many points fall in ANY PA?


  TD_ll_list[[i]] <- TD_ll

}

TD_all <- rbind.data.frame(TD_ll_list[[1]]@data, TD_ll_list[[2]]@data)

## summarize total points, and proportion of points ##
gen_summ <- TD_all %>% group_by(period) %>% summarise(
  n_ind   = n_distinct(ID),
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
  prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
)

## Calculate by foraging area (for foraging period only) ## --------------------
dest_summ <- TD_all %>% group_by(period, destination) %>% summarise(
  n_ind   = n_distinct(ID),
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
  prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
)


## save ## ~~~
filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_gen.csv")
fwrite(gen_summ, filename)

filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_bydest.csv")
fwrite(dest_summ, filename)

## quick check of points and MPAs for specific region
# tracks_sf <- TD_ll %>% st_as_sf(crs = 4326, agr = "constant") %>% filter(destination == "Bijagos")
# mapview(mpas) + mapview::mapview(tracks_sf, col.regions="red") 

## Averages per individual ## 
## summarize total points, and proportion of points ##
id_summ <- TD_all %>% group_by(period, ID, destination) %>% summarise(
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total w/in any CONSERVATION AREA
  prop_in_ca = pnts_in_ca/n_pnts # prop. of locations w/in CONSERVATION AREA
) 

gen_id_summ <- id_summ %>% group_by(period) %>% summarise(
  n_ind    = n_distinct(ID),
  tot_pnts = sum(n_pnts),
  m_pnts   = mean(n_pnts),
  sd_pnts  = sd(n_pnts),
  m_pnts_in_pa  = mean(pnts_in_pa), # Protected areas only
  sd_pnts_in_pa = sd(prop_in_pa),
  m_prop_in_pa  = mean(prop_in_pa),
  sd_prop_in_pa = sd(prop_in_pa),
  m_pnts_in_ca  = mean(pnts_in_ca), # all conservation areas
  sd_pnts_in_ca = sd(prop_in_ca),   # including BABR
  m_prop_in_ca  = mean(prop_in_ca),
  sd_prop_in_ca = sd(prop_in_ca)
)

dest_id_summ <- id_summ %>% group_by(period, destination) %>% summarise(
  n_ind    = n_distinct(ID),
  tot_pnts  = sum(n_pnts),
  m_pnts  = mean(n_pnts),
  sd_pnts  = sd(n_pnts),
  m_pnts_in_pa = mean(pnts_in_pa), # Protected areas only
  sd_pnts_in_pa = sd(prop_in_pa),
  m_prop_in_pa = mean(prop_in_pa),
  sd_prop_in_pa = sd(prop_in_pa),
  m_pnts_in_ca  = mean(pnts_in_ca), # all conservation areas
  sd_pnts_in_ca = sd(prop_in_ca),   # including BABR
  m_prop_in_ca  = mean(prop_in_ca),
  sd_prop_in_ca = sd(prop_in_ca)
)

## save ## ~~~
filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_byiddest.csv")
# fwrite(dest_id_summ, "data/analysis/summaries/pnts_in_pa_byiddest.csv")

filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_byidgen.csv")
# fwrite(gen_id_summ,  "data/analysis/summaries/pnts_in_pa_byidgen.csv")

filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_gen.csv")
# fwrite(gen_id_summ, "data/analysis/summaries/pnts_in_pa_byiddest_internest.csv")


#------------------------------------------------------------------------------#
## Sensitivity analysis comparing PTT only data vs GPS only data ##

gpsids <- unique(TD_all[TD_all$sensor=="GPS",]$ID)
pttids <- unique(TD_all[!TD_all$sensor=="GPS",]$ID)

## Calculate by device type ## --------------------
gen_summ_s <- TD_ll@data %>% group_by(period, sensor) %>% summarise(
  n_ind   = n_distinct(ID),
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
  prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
)
gen_summ_s

dest_summ_s <- TD_ll@data %>% group_by(period,destination, sensor) %>% summarise(
  n_ind   = n_distinct(ID),
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
  prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
)
dest_summ_s

## summarize total points, and proportion of points ##
id_summ_s <- TD_ll@data %>% group_by(period, ID, destination, sensor) %>% summarise(
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total w/in any CONSERVATION AREA
  prop_in_ca = pnts_in_ca/n_pnts # prop. of locations w/in CONSERVATION AREA
) 

gen_id_summ_s <- id_summ_s %>% group_by(period, sensor) %>% summarise(
  n_ind    = n_distinct(ID),
  tot_pnts = sum(n_pnts),
  m_pnts   = mean(n_pnts),
  sd_pnts  = sd(n_pnts),
  m_pnts_in_pa  = mean(pnts_in_pa), # Protected areas only
  sd_pnts_in_pa = sd(prop_in_pa),
  m_prop_in_pa  = mean(prop_in_pa),
  sd_prop_in_pa = sd(prop_in_pa),
  m_pnts_in_ca  = mean(pnts_in_ca), # all conservation areas
  sd_pnts_in_ca = sd(prop_in_ca),   # including BABR
  m_prop_in_ca  = mean(prop_in_ca),
  sd_prop_in_ca = sd(prop_in_ca)
)
gen_id_summ_s

dest_id_summ_s <- id_summ_s %>% group_by(period, destination, sensor) %>% summarise(
  n_ind    = n_distinct(ID),
  tot_pnts  = sum(n_pnts),
  m_pnts  = mean(n_pnts),
  sd_pnts  = sd(n_pnts),
  m_pnts_in_pa = mean(pnts_in_pa), # Protected areas only
  sd_pnts_in_pa = sd(prop_in_pa),
  m_prop_in_pa = mean(prop_in_pa),
  sd_prop_in_pa = sd(prop_in_pa),
  m_pnts_in_ca  = mean(pnts_in_ca), # all conservation areas
  sd_pnts_in_ca = sd(prop_in_ca),   # including BABR
  m_prop_in_ca  = mean(prop_in_ca),
  sd_prop_in_ca = sd(prop_in_ca)
)
dest_id_summ_s

#------------------------------------------------------------------------------#
## Sensitivity analysis (2) comparing PTT pnts to GPS pnts for IDs w/ both ##

TD_gpsids <- TD_all %>% filter(ID %in% gpsids)

## Calculate by device type ## --------------------
gen_summ_s2 <- TD_gpsids %>% group_by(period,sensor) %>% summarise(
  n_ind   = n_distinct(ID),
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
  prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
)
gen_summ_s2

dest_summ_s2 <- TD_gpsids %>% group_by(period,destination, sensor) %>% summarise(
  n_ind   = n_distinct(ID),
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
  prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
)
dest_summ_s2
## by ID ## 
id_summ_s2 <- TD_gpsids %>% group_by(period,ID, destination, sensor) %>% summarise(
  n_pnts  = n(),              # total number of tracking locations
  pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
  prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
  pnts_in_ca = sum(over_CA),     # total w/in any CONSERVATION AREA
  prop_in_ca = pnts_in_ca/n_pnts # prop. of locations w/in CONSERVATION AREA
) 
dest_id_summ_s2 <- id_summ_s2 %>% group_by(period,destination, sensor) %>% summarise(
  n_ind    = n_distinct(ID),
  tot_pnts  = sum(n_pnts),
  m_pnts  = mean(n_pnts),
  sd_pnts  = sd(n_pnts),
  m_pnts_in_pa = mean(pnts_in_pa), # Protected areas only
  sd_pnts_in_pa = sd(prop_in_pa),
  m_prop_in_pa = mean(prop_in_pa),
  sd_prop_in_pa = sd(prop_in_pa),
  m_pnts_in_ca  = mean(pnts_in_ca), # all conservation areas
  sd_pnts_in_ca = sd(prop_in_ca),   # including BABR
  m_prop_in_ca  = mean(prop_in_ca),
  sd_prop_in_ca = sd(prop_in_ca)
)
dest_id_summ_s2
