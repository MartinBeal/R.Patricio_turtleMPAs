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
      # TD <- subset(TD, ID %in% KDE_ids)
      
    } else if(period == "foraging"){
      TD$period <- rep(period)
      # only keep IDs which have had enough data to estimate KDE for 
      # KDE_ids <- names(KDE <- readRDS(paste0("data/analysis/UDs/ind_areas/individual_UDs_h2.08_c1_foraging.rds")))
    }
    
    TD_ll <- spTransform(TD, CRS("+init=epsg:4326")) # back to lat/lon
    
    ### estimate proportion of positions within and outside MPAs ##---------------
    mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp")
    babr <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/BABR_polygon.shp")
  
    # combine and union BABR with MPAs dataset
    cons_areas <- aggregate(rbind(mpas[,-c(31:32)], babr))
    # mapview(cons_areas)
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
    TD_ll$over_CA <- ifelse(is.na(over_cons), F, T)
    table(TD_ll$over_CA)    # how many points fall in ANY PA?
    
    
    TD_ll_list[[i]] <- TD_ll
    
  }
  
  TD_all <- rbind.data.frame(TD_ll_list[[1]]@data, TD_ll_list[[2]]@data)
  
  ## remove PA overlap of two individuals which barely overlap due to error
  TD_all$over_PA <- ifelse(TD_all$ID %in% c("197138", "182453"), 
                           FALSE, TD_all$over_PA)
  TD_all$over_CA <- ifelse(TD_all$ID %in% c("197138", "182453"), 
                           FALSE, TD_all$over_CA)
  ### Calculating MPA coverage of tracking points in various ways -------------
  
  ## summarize total points, and proportion of points (including IDs which don't visit at all ##
  gen_summ <- TD_all %>% group_by(period) %>% summarise(
    n_ind   = n_distinct(ID),
    n_pnts  = n(),                  # total number of tracking locations
    pnts_in_pa = sum(over_PA),      # total number of locations falling w/in any PA
    prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
    pnts_in_ca = sum(over_CA),      # total number of locations falling w/in any PA
    prop_in_ca = pnts_in_ca/n_pnts  # proportion of locations w/in PAs
  )
  
  ## Calculate by foraging area (for foraging period only) ## --------------------
  dest_summ <- TD_all %>% group_by(period, destination) %>% summarise(
    n_ind   = n_distinct(ID),
    n_pnts  = n(),                  # total number of tracking locations
    pnts_in_pa = sum(over_PA),      # total number of locations falling w/in any PA
    prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
    pnts_in_ca = sum(over_CA),      # total number of locations falling w/in any PA
    prop_in_ca = pnts_in_ca/n_pnts  # proportion of locations w/in PAs
  ) %>% filter(period == "foraging")
  
  ## quick check of points and MPAs for specific region
  # tracks_sf <- TD_ll %>% st_as_sf(crs = 4326, agr = "constant") %>% filter(destination == "Bijagos")
  # mapview(mpas) + mapview::mapview(tracks_sf, col.regions="red") 
  
  ## Averages per individual ## ----------------------------------------------
  ## summarize total points, and proportion of points ##
  id_summ <- TD_all %>% group_by(period, ID, destination) %>% summarise(
    n_pnts  = n(),              # total number of tracking locations
    pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
    prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
    pnts_in_ca = sum(over_CA),     # total w/in any CONSERVATION AREA
    prop_in_ca = pnts_in_ca/n_pnts # prop. of locations w/in CONSERVATION AREA
  ) 
  
  ## averages for all IDs
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
  ## averages for only IDs which visited PAs at all 
  gen_id_summ2 <- id_summ %>% filter(pnts_in_pa>0) %>% group_by(period) %>% summarise(
    n_ind    = n_distinct(ID),
    tot_pnts = sum(n_pnts),
    m_pnts   = mean(n_pnts),
    sd_pnts  = sd(n_pnts),
    m_pnts_in_pa  = mean(pnts_in_pa), # Protected areas only
    sd_pnts_in_pa = sd(prop_in_pa),
    m_prop_in_pa  = mean(prop_in_pa),
    max_prop_in_pa = max(prop_in_pa),
    min_prop_in_pa = min(prop_in_pa),
    sd_prop_in_pa = sd(prop_in_pa),
    m_pnts_in_ca  = mean(pnts_in_ca), # all conservation areas
    sd_pnts_in_ca = sd(prop_in_ca),   # including BABR
    m_prop_in_ca  = mean(prop_in_ca),
    sd_prop_in_ca = sd(prop_in_ca),
    max_prop_in_ca = max(prop_in_ca),
    min_prop_in_ca = min(prop_in_ca),
  )
  
  gen_id_summ_jn <- left_join(gen_id_summ, gen_id_summ2, by=c("period"))
  
  ## averages for all IDs
  dest_id_summ <- id_summ %>% group_by(period, destination) %>% summarise(
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
  ) %>% filter(period == "foraging")
  
  ## averages for only IDs which visited PAs at all 
  dest_id_summ2 <- id_summ %>% filter(pnts_in_pa>0) %>% group_by(period, destination) %>% 
    summarise(
      n_ind    = n_distinct(ID),
      tot_pnts = sum(n_pnts),
      m_pnts   = mean(n_pnts),
      sd_pnts  = sd(n_pnts),
      m_pnts_in_pa = mean(pnts_in_pa), # Protected areas only
      sd_pnts_in_pa = sd(prop_in_pa),
      m_prop_in_pa = mean(prop_in_pa),
      sd_prop_in_pa = sd(prop_in_pa),
      m_pnts_in_ca  = mean(pnts_in_ca), # all conservation areas
      sd_pnts_in_ca = sd(prop_in_ca),   # including BABR
      m_prop_in_ca  = mean(prop_in_ca),
      sd_prop_in_ca = sd(prop_in_ca)
    ) %>% filter(period == "foraging")
  
  dest_id_summ_jn <- left_join(dest_id_summ, dest_id_summ2, by=c("period", "destination"))
  
  ## save ## ~~~
  filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_id.csv")
  fwrite(id_summ, filename)
  
  filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_byiddest.csv")
  fwrite(dest_id_summ_jn, filename)
  
  filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_byidgen.csv")
  fwrite(gen_id_summ_jn,  filename)
  
  
  ## ----------------------------------------------
  ## Averages ignoring individual, first removing individuals not using PAs
  gen_summ2 <- TD_all %>% left_join(id_summ) %>% filter(pnts_in_pa>0) %>% 
    group_by(period) %>% summarise(
      n_ind   = n_distinct(ID),
      n_pnts  = n(),              # total number of tracking locations
      pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
      prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
      pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
      prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
  )
  
  gen_summ_jn <- left_join(gen_summ, gen_summ2, by=c("period"))
  
  dest_summ2 <- TD_all %>% left_join(id_summ) %>% filter(pnts_in_pa>0) %>% 
    group_by(period, destination) %>% summarise(
      n_ind   = n_distinct(ID),
      n_pnts  = n(),              # total number of tracking locations
      pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
      prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
      pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
      prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
  ) %>% filter(period == "foraging")
  
  dest_summ_jn <- left_join(dest_summ, dest_summ2, by=c("period", "destination"))
  
  ## save ## ~~~
  filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_gen.csv")
  fwrite(gen_summ_jn, filename)
  
  filename <- paste0("data/analysis/summaries/", datatype, "/pnts_in_pa_bydest.csv")
  fwrite(dest_summ_jn, filename)
  
  #------------------------------------------------------------------------------#
  ## Sensitivity analysis (2) comparing PTT pnts to GPS pnts for IDs w/ both ##
  if(datatype=="raw_filtered"){
    
    gpsids <- unique(TD_all[TD_all$sensor=="GPS",]$ID)
    # pttids <- unique(TD_all[!TD_all$sensor=="GPS",]$ID)
    
    TD_gpsids <- TD_all %>% filter(ID %in% gpsids)
    
    ## Calculate by device type ## --------------------
    gen_summ_s <- TD_gpsids %>% group_by(period,sensor) %>% summarise(
      n_ind   = n_distinct(ID),
      n_pnts  = n(),              # total number of tracking locations
      pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
      prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
      pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
      prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
    )
    gen_summ_s
    
    dest_summ_s <- TD_gpsids %>% group_by(period,destination, sensor) %>% summarise(
      n_ind   = n_distinct(ID),
      n_pnts  = n(),              # total number of tracking locations
      pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
      prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
      pnts_in_ca = sum(over_CA),     # total number of locations falling w/in any PA
      prop_in_ca = pnts_in_ca/n_pnts # proportion of locations w/in PAs
    ) %>% filter(period == "foraging")
    dest_summ_s
    
    ## save ## ~~~
    filename <- paste0("data/analysis/summaries/sensitivity_analysis/", datatype, "_pnts_in_pa_gen.csv")
    fwrite(gen_summ_s, filename)
    
    filename <- paste0("data/analysis/summaries/sensitivity_analysis/", datatype, "_pnts_in_pa_bydest.csv")
    fwrite(dest_summ_s, filename)
    
    ## by ID ## 
    id_summ_s <- TD_gpsids %>% group_by(period,ID, destination, sensor) %>% summarise(
      n_pnts  = n(),              # total number of tracking locations
      pnts_in_pa = sum(over_PA),     # total number of locations falling w/in any PA
      prop_in_pa = pnts_in_pa/n_pnts, # proportion of locations w/in PAs
      pnts_in_ca = sum(over_CA),     # total w/in any CONSERVATION AREA
      prop_in_ca = pnts_in_ca/n_pnts # prop. of locations w/in CONSERVATION AREA
    ) 
    dest_id_summ_s <- id_summ_s %>% group_by(period,destination, sensor) %>% summarise(
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
    ) %>% filter(period == "foraging")
    dest_id_summ_s2 <- id_summ_s %>% filter(pnts_in_pa>0) %>% 
      group_by(period,destination, sensor) %>% 
      summarise(
        n_ind    = n_distinct(ID),
        tot_pnts = sum(n_pnts),
        m_pnts   = mean(n_pnts),
        sd_pnts  = sd(n_pnts),
        m_pnts_in_pa = mean(pnts_in_pa), # Protected areas only
        sd_pnts_in_pa = sd(prop_in_pa),
        m_prop_in_pa = mean(prop_in_pa),
        sd_prop_in_pa = sd(prop_in_pa),
        m_pnts_in_ca  = mean(pnts_in_ca), # all conservation areas
        sd_pnts_in_ca = sd(prop_in_ca),   # including BABR
        m_prop_in_ca  = mean(prop_in_ca),
        sd_prop_in_ca = sd(prop_in_ca)
    ) %>% filter(period == "foraging")
    
    dest_id_summ_s_jn <- left_join(dest_id_summ_s, dest_id_summ_s2, by=c("period", "destination", "sensor"))
    
    ## Save ## 
    filename <- paste0("data/analysis/summaries/sensitivity_analysis/", datatype, "_pnts_in_pa_byiddest.csv")
    fwrite(dest_id_summ_s_jn, filename)
  
  }
}

