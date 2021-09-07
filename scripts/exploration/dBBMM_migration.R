## estimate migration routes using dynamic BBMM ## ----------------------------

pacman::p_load(track2KBA, dplyr, sp, ggplot2, raster, mapview)

source("C:\\Users\\Martim Bill\\Documents\\R\\source_scripts\\UD_fxns.R")

tictoc::tic()

  folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods
  
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
  
  ## filter out individuals w/ <10 points
  tracks <- tracks %>% filter(!ID %in% c("182459", "60890", "60893"))
  
  tracks <- tracks[order(tracks$ID, tracks$DateTime), ] ## order date time stamps within each individual
  
  ## estimate horizontal error in meters ##
  
  tracks$error_m <- ifelse(tracks$sensor == "Argos doppler shift" & tracks$Argos.Location.Quality == "0", 1500, 
         ifelse(tracks$sensor == "Argos doppler shift" & tracks$Argos.Location.Quality == "1", 750, 
                ifelse(tracks$sensor == "Argos doppler shift" & tracks$Argos.Location.Quality == "2", 375, 
                       ifelse(tracks$sensor == "Argos doppler shift" & tracks$Argos.Location.Quality == "3", 125, tracks$GPS.HDOP * 15) ## NEED TO CHECK REPORTED AVG ERROR FOR FASTLOC GPS
                              )))
  
  # move object
  tracks_m <- move::move(
    x      = tracks$Longitude, 
    y      = tracks$Latitude, 
    time   = tracks$DateTime, 
    proj   = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),  # common project. for all individuals, custom set by tripSplit
    animal = tracks$ID
  )
  
  tracks_m@data$error_m <- tracks$error_m
  
  tracks_m <- sp::spTransform(tracks_m, center=T)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Estimate dBBMM
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # 'odRaster' is the raster surface from simulate_GLS_tracks, this is used here just for comparisons sake
  # aMask <- raster(extent(tracks_m), crs=crs(tracks_m), ncols = 50, nrows = 100)
  aMask <- raster(extent(tracks_m), crs=crs(tracks_m), resolution = 10000)
  
  aMask[] <- 1  # UDLe
  aMask_ext <- extend(aMask,c(50,50)) # need to extend raster to fit all of dBBMM estimate inside
  
  vardbbmm <- move::brownian.motion.variance.dyn(tracks_m,  # calculate variance across each step in the data 
                                                 location.error=tracks_m$error_m, 
                                                 margin=3, 
                                                 window.size=7
  )
  str(vardbbmm)
  
  # now we can manipulate
  xx <- unlist(lapply(move::timeLag(tracks_m,"hours"), function(x){
    c(NA, x)
  }))
  
  vardbbmm@interest[xx > 24] <- FALSE # 50 (h) = time lag are ignored in dBBM calc.
  
  GLSdbb <- move::brownian.bridge.dyn(vardbbmm, #GLSm
                                      raster=aMask_ext, 
                                      # ext=5, 
                                      location.error=tracks_m$error_m, 
                                      margin=3,
                                      window.size=7,
                                      verbose=T
  )
  plot(GLSdbb[[3]])
  GLSdbb_ud <- move::getVolumeUD(GLSdbb) 
  plot(GLSdbb_ud)
  
  GLSdbb_ud95 <- GLSdbb_ud
  GLSdbb_ud95[GLSdbb_ud95 > .99] <- NA
  # GLSdbb_ud95 <- trim(GLSdbb_ud95) # move NAs
  
  GLSdbb_ud95 <- crop(GLSdbb_ud95, aMask)
  plot(GLSdbb_ud95)
  
  mapview::mapview(GLSdbb_ud95[[4]])

  ## average together indvidual mig UDs ---------------------------------------
  comb_dbb <- raster::calc(GLSdbb, fun = mean)
  ## convert to cumulative density function
  comb_cdf <- comb_dbb
  rank <- (1:length(values(comb_cdf)))[rank(values(comb_cdf))]
  values(comb_cdf) <- 1 - cumsum(sort(values(comb_cdf)))[rank]
  ## isopleth levels
  comb_cdf99 <- comb_cdf
  values(comb_cdf99) <- ifelse(values(comb_cdf99) > .99, NA, values(comb_cdf99))
  plot(comb_cdf99)
  
  comb_cdf50 <- comb_cdf
  values(comb_cdf50) <- ifelse(values(comb_cdf50) > .5, NA, values(comb_cdf50))
  plot(comb_cdf50)
  
  mapview(comb_cdf50)
  
  ## count number of overlapping individuals' isopleth areas ------------------
  
  ind_cdf <- split(GLSdbb)
  n <- length(ind_cdf)
  UDlevel <- 95
  
  Noverlaps <- lapply(split(GLSdbb), function(x) {
    
    rank <- (1:length(values(x)))[rank(values(x))]
    values(x) <- 1 - cumsum(sort(values(x)))[rank]
    values(x) <- ifelse(
      values(x) > UDlevel/100, 
      0, 1)
    
    return(x)
    
  })
  
  comb_cnt <- raster::calc(stack(Noverlaps), fun=sum)
  values(comb_cnt) <- ifelse(values(comb_cnt) == 0, NA, values(comb_cnt))
    
  mapview(comb_cnt)
  
  saveRDS(comb_cnt, 
          paste0("data/analysis/mig_grid/rawdata_dbbmm_mig_grid_10x10_", 
                 "UD", UDlevel, "_", n, "_nIDs.rds"))
  
  values(comb_cnt) <- (values(comb_cnt)/n)*100
  
  saveRDS(comb_cnt, 
          paste0("data/analysis/mig_grid/rawdata_dbbmm_mig_grid_10x10_", 
                 "UD", UDlevel, "_", n, "_perc.rds"))
  
# }
  
tictoc::toc()
  