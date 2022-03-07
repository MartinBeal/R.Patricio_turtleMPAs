### Produce gridded density raster of proportion of migrating individuals passing through certain areas ### 
# proportion estimated from simple overlap of raw or interpolated points

### estimate proportion of positions within and outside MPAs ##----------------

pacman::p_load(track2KBA, dplyr, sp, ggplot2, raster, mapview)

datatypes <- c("raw_filtered","interpolated")

for(y in seq_along(datatypes)){
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
  
  prj <- TD@proj4string
  
  xtnt <- extent(TD)
  
  r <- raster(xtnt, crs=prj)
  res(r) <- c(10000, 10000)
  
  TD$cellid <- cellFromXY(r,TD)
  ovrsumm <- TD@data %>% 
    group_by(cellid) %>% 
    summarise(
      n_IDs  = n_distinct(ID),
      n_pnts = n()) %>% arrange(cellid)
  
  r2 <- r
  r2[ovrsumm$cellid] <- ovrsumm$n_IDs
  
  mapview(r2)
  
  
  ### SAVE ### ------------
  saveRDS(r2, paste0("data/analysis/mig_grid/", datatype, "_mig_grid_10x10_", n, ".rds"))
}
