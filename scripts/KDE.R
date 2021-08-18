## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Analysis of green turtle movements post-laying ## 

pacman::p_load(track2KBA, amt, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter)

datatypes <- c("raw_filtered","interpolated")

## Data input ~~~~~~~~~~~~~~~~~~
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
  
  for(i in seq_along(periods)){
      print(i)
      ## analyze inter-nesting, foraging or migration?
      period <- periods[i]
  
      tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
      
      if("id" %in% colnames(tracks)){
        tracks <- formatFields(tracks,
                               fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date") 
      }
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ## track2KBA analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
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
      }
      
      # estimate smoothing parameter candidates (output in km)
      HVALS <- findScale(TD, scaleARS = F)
      HVALS
      
      # choose one (reference parameter used here)
      h <- HVALS$href
      h
      
      cres <- 1 ## cell resoluation (kmXkm)
      # calculate UDs for each individual
      KDE  <- estSpaceUse(TD, scale=h, polyOut=F, res = cres)
      UD50 <- estSpaceUse(TD, scale=h, polyOut=T, levelUD = 50, res = cres)
      UD95 <- estSpaceUse(TD, scale=h, polyOut=T, levelUD = 95, res = cres)
      mapKDE(UD50$UDPolygons)
      mapKDE(UD95$UDPolygons)
      
      
      # represent <- repAssess(TD, KDE, iteration=100, nCores=2, levelUD = 50, avgMethod = "mean", bootTable = F)
      
      ## Save ## ------------------------------------------------------------------
      saveRDS(KDE,  
              paste0( paste0("data/analysis/UDs/", datatype, "/ind_areas", "/individual_UDs_h"),
                      h, "_c",  cres, "_", period, ".rds") )
      saveRDS(UD50,  
              paste0( paste0("data/analysis/UDs/", datatype, "/ind_areas", "/individual_UD50_h"),
                      h, "_c",  cres, "_", period, ".rds") )
      saveRDS(UD95,  
              paste0( paste0("data/analysis/UDs/", datatype, "/ind_areas", "/individual_UD95_h"),
                      h, "_c",  cres, "_", period, ".rds") )
      
      ## average together individual UDs weighted by n_pnts ## 
      
      ## convert UDs from estUDm to rasterStack object ##
      KDEraster <- raster::stack(lapply(KDE, function(x) {
        raster::raster(as(x, "SpatialPixelsDataFrame"), values=TRUE)
      } ))
      
      # check that N KDEs == N individuals (if not, remove individuals from tracks b4 next step)
      nlayers(KDEraster)
      n_distinct(tracks$ID)
      
      id_w_kde <- unique(TD$ID)[validNames(unique(TD$ID)) %in% names(KDEraster)]
      
      # weighted mean - number of points per ID -------------------------------
      weights <- tracks %>% 
        dplyr::filter(ID %in% id_w_kde) %>%
        group_by(ID) %>% summarise(n_pnts = n()) %>% pull(n_pnts)
      KDEcmbnd_w <- raster::weighted.mean(KDEraster, w=weights)   
      # arithmetic mean - all individuals equally weighted --------------------
      KDEcmbnd_a <- raster::calc(KDEraster, mean) # arithmetic mean
      
      ## compare ##
      # dev.new()
      plot(KDEcmbnd_w)
      # dev.new()
      plot(KDEcmbnd_a)
      
      ## Save ## --------------------------------------------------------------
      saveRDS(KDEcmbnd_w, 
              paste0( paste0("data/analysis/UDs/", datatype, "/PMF", "/w_groupUD_h"),
                      h, "_c",  cres, "_", period, ".rds"))
      saveRDS(KDEcmbnd_a, 
              paste0( paste0("data/analysis/UDs/", datatype, "/PMF", "/a_groupUD_h"),
                      h, "_c",  cres, "_", period, ".rds"))
      
      
      ## Derive %UD from full UDs ## -----------------------------------------------
      # KDEcmbnd_w <- readRDS(paste0("data/analysis/UDs/w_groupUD_h", h, "_c",  cres, "_", period, ".rds"))
      # KDEcmbnd_a <- readRDS(paste0("data/analysis/UDs/a_groupUD_h", h, "_c",  cres, "_", period, ".rds"))
      
      ### Convert UD in PMF form to CMF (i.e. to % UD) ### -------------------------
      CUD <- KDEcmbnd_w
      # CUD <- KDEcmbnd_a
      
      pixArea <- raster::res(CUD)[1]
      levelUD <- c(50, 95)
      
      percUDs <- lapply(1:2, function(x) {
        lvl <- levelUD[x]
        df <- data.frame(UD = raster::getValues(CUD)) %>%
          mutate(
            rowname = seq_len(length(raster::getValues(CUD))),
            usage = .data$UD * (pixArea^2)
          ) %>%
          arrange(desc(.data$usage)) %>%
          mutate(cumulUD = cumsum(.data$usage),
                 INSIDE = ifelse(.data$cumulUD < (lvl/100), 1, NA)
          ) %>%
          arrange(.data$rowname) %>%
          dplyr::select(.data$INSIDE)
        CUD[] <- df$INSIDE
        
        return(CUD)
      })
      
      names(percUDs) <- levelUD
      
      plot(percUDs[[1]])
      plot(percUDs[[2]])
      
      
      ## Save ## --------------------------------------------------------------------
      saveRDS(percUDs, 
              paste0(paste0("data/analysis/UDs/", datatype, "/CDF", "/w_groupCDFs_h"), 
                     h, "_c",  cres, "_", period, ".rds"))
      
      saveRDS(percUDs, 
              paste0(paste0("data/analysis/UDs/", datatype, "/CDF", "/a_groupCDFs_h"), 
                            h, "_c",  cres, "_", period, ".rds"))
  }
}
