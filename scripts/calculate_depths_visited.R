## Calculate depths visited by turtles ## 

pacman::p_load(
  track2KBA, amt, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter, raster, data.table)

datatypes <- c("raw_filtered","interpolated")
# datatypes <- c("raw_filtered","interpolated", "sensitivity_analysis")

for(y in seq_along(datatypes)){
  datatype <- datatypes[y]
  period <- "foraging"
  
  if(datatype == "interpolated"){
    TDfolder <- "data/analysis/interpolated/" # repository w/ datasets split into periods
  } else if(datatype == "raw_filtered"){
    TDfolder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods
  }
  # choose whether to include sub-6 satellite data or naw
  satfilt <- "satfilt6"
  
  tracks <- readRDS(paste0(TDfolder, period, "_", satfilt, ".rds"))
  
  if("id" %in% colnames(tracks)){
    tracks <- formatFields(tracks,
                           fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date") 
  }
  
  TD_ll <- st_as_sf(x = tracks, 
                    coords = c("Longitude", "Latitude"),
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% as_Spatial()
  
  TD_prj <- projectTracks(tracks, "azim", custom = T)
  prj    <- crs(TD_prj)
  
  # get UD areas
  UDfolder <- paste0("data/analysis/UDs/", datatype, "/ind_areas")
  UD50 <- readRDS(
    list.files(
      UDfolder, full.names=T)[str_detect(
        list.files(UDfolder), pattern=paste(c(period, "UD50"),collapse = '|'))][1])   
  UD95 <- readRDS(
    list.files(
      UDfolder, full.names=T)[str_detect(
        list.files(UDfolder), pattern=paste(c(period, "UD95"),collapse = '|'))][1])   
  
  ## load bathymetry data ##
  # depth <- raster("data/geodata/GEBCO_2020_W Africa/gebco_2020_n20.9_s10.1_w-18.9_e-14.8.tif")
  depth     <- raster("data/geodata/ETOPO1.tiff") # original ETOPO
  # depth     <- raster("data/geodata/ETOPO1_patchedw_bathyrastPNBA.tif") # ETOPO + soundins in MAU
  depth_prj <- projectRaster(depth, crs=prj)
  # nodepthzone <- raster::shapefile("data/geodata/Bathy_BA/no_depth_info_zone.shp")
  
  # plot(depth)
  # plot(nodepthzone, add=T)
  
  ## create depth raster where cells falling within unsampled (untrustworthy) zone are changed to NA ##
  # depth2 <- mask(depth, nodepthzone, inverse=T)
  
  ## extract depth for each tracking point
  TD_ll$depth <- extract(depth, TD_ll)
  hist(TD_ll$depth)
  # mapview(TD_ll)
  TD_prj$depth <- extract(depth_prj, TD_prj)
  
  ## extract depth from cells falling within each individual's UD (both 50 and 95%) 
  # UD50_ll <- UD50$UDPolygons %>% as_Spatial() %>% spTransform(CRS("+init=epsg:4326"))
  # UD95_ll <- UD95$UDPolygons %>% as_Spatial() %>% spTransform(CRS("+init=epsg:4326"))
  
  # UD50_depths <- extract(depth, UD50_ll)
  # UD95_depths <- extract(depth, UD95_ll)
  
  UD50_prj <- UD50$UDPolygons %>% as_Spatial() %>% spTransform(prj)
  UD95_prj <- UD95$UDPolygons %>% as_Spatial() %>% spTransform(prj)
  
  UD50_depths <- extract(depth, UD50_prj)
  UD95_depths <- extract(depth, UD95_prj)
  ## Turn list into a dataframe, where each row is summary of depths used by each individual ##
  # NOTE that 'min' refers to the deepest part (i.e. most negative)
  # UD 50%
  depths_50 <- rbindlist(lapply(UD50_depths, function(x){ 
    data.frame(
      min = min(x),
      mn  = mean(x),
      sd  = sd(x),
      md  = median(x)
    )
  })) %>% bind_cols(UD50_prj@data)
  
  # UD 95%
  depths_95 <- rbindlist(lapply(UD95_depths, function(x){ 
    data.frame(
      min = min(x),
      mn  = mean(x),
      sd  = sd(x),
      md  = median(x)
    )
  })) %>% bind_cols(UD95_prj@data)
  
  ## summarize across individuals - 'min-min' == deepest ##
  summ_50 <- summarise(depths_50, 
                       mn_mn  = mean(mn),                        # avg avg depth 1
                       sd_mn = sd(mn), min_mn = min(mn),         # sd and min
                       mn_md  = mean(md),                        # avg avg depth 2
                       sd_md = sd(md), min_md = min(md),         # sd and min
                       mn_min = mean(min), sd_min = sd(min),     # average min depth 1
                       md_min = median(min),  min_min = min(min), # average min depth 2
                       IQR_min = IQR(min)
  )
  summ_50
  summ_95 <- summarise(depths_95, 
                       mn_mn  = mean(mn),                        # avg avg depth 1
                       sd_mn = sd(mn), min_mn = min(mn),         # sd and min
                       mn_md  = mean(md),                        # avg avg depth 2
                       sd_md = sd(md), min_md = min(md),         # sd and min
                       mn_min = mean(min), sd_min = sd(min),     # average min depth 1
                       md_min = median(min),  min_min = min(min), # average min depth 2
                       IQR_min = IQR(min)
  )
  summ_95
  ## compare with depth-use calculated directly from tracking locations
  # summarise(TD_ll@data, 
  #           mn_min = mean(depth),
  #           sd_min = sd(depth),
  #           min_min = min(depth))
  # by_ind <- TD_ll@data %>% group_by(ID) %>%
  by_ind <- TD_prj@data %>% group_by(ID) %>%
    summarise( 
      mn = mean(depth),   sd = sd(depth),
      md = median(depth), min = min(depth)) 
  by_ind
  crss_ind <- by_ind %>% summarise( 
    mn_mn  = mean(mn),                        # avg avg depth 1
    sd_mn = sd(mn), min_mn = min(mn),         # sd and min
    mn_md  = mean(md),                        # avg avg depth 2
    sd_md = sd(md), min_md = min(md),         # sd and min
    mn_min = mean(min), sd_min = sd(min),     # average min depth 1
    md_min = median(min),  min_min = min(min), # average min depth 2
    IQR_min = IQR(min)
  )
  crss_ind
  
  ## SAVE ##
  filename <- paste0("data/analysis/summaries/", datatype, "/ind_depths_used_fromTD.csv")
  fwrite(by_ind, filename)
  filename <- paste0("data/analysis/summaries/", datatype, "/depths_used_avgfromTD.csv")
  fwrite(crss_ind, filename)
  filename <- paste0("data/analysis/summaries/", datatype, "/depths_used_avgfromUD50.csv")
  fwrite(summ_50,  filename)
  filename <- paste0("data/analysis/summaries/", datatype, "/depths_used_avgfromUD95.csv")
  fwrite(summ_95,  filename)
  
  
  ### Calculate the amount of available depth area inside and outside MPAs ### ---
  mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp")
  # mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021.shp")
  babr <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/BABR_polygon.shp")
  
  mpas <- st_as_sf(mpas) %>% st_transform(crs=proj4string(TD_prj))
  babr <- st_as_sf(babr) %>% st_transform(crs=proj4string(TD_prj))
  
  ## take minimum avg across individuals as maximum depth
  # maxdepth <- crss_ind$min_mn # points
  maxdepth <- crss_ind$min_min # points
  # maxdepth <- summ_95$min_min    # 2UD 95%
  # maxdepth <- summ_95$min_mn    # UD 95%
  
  maxdepth_rndnd <- round(maxdepth) # rounded value
  
  hab <- depth
  hab <- depth_prj
  hab <- setValues(hab, 
                   ifelse( (values(hab) < maxdepth)|(values(hab) >= 0), 
                           NA, values(hab)))
  
  inpas <- raster::intersect(hab, mpas)
  plot(inpas)

  mpas_r <- rasterize(mpas, hab, field=rep(1, nrow(mpas)))
  # plot(mpas_r)
  
  ## habitat within PAs
  inpas <- setValues(hab, 
                     ifelse( (values(mpas_r) == 0), 
                             NA, values(hab)))
  plot(inpas)
  plot(mpas$geometry, add=T)
  ## habitat outside PAs 
  outpas <- setValues(hab, 
                      ifelse( !is.na(values(mpas_r) ), 
                              NA, values(hab)))
  plot(outpas)
  plot(mpas$geometry, add=T)
  
  ## calculate areas ## --------------------------------------------------------
  
  ## first, entire suitable depth habitat ##
  csizes_hab   <- raster::area(hab, na.rm=TRUE, weights=FALSE)
  csizes_inpa  <- raster::area(inpas, na.rm=TRUE, weights=FALSE)
  csizes_outpa <- raster::area(outpas, na.rm=TRUE, weights=FALSE)
  length(na.omit(values(hab))) * res(hab)[1]*res(hab)[2]
  length(na.omit(values(inpas))) * res(inpas)[1]*res(inpas)[2]
  length(na.omit(values(outpas))) * res(outpas)[1]*res(outpas)[2]
  
  res(hab)[1]*res(hab)[2]
  
  #delete NAs from vector of all raster cells
  csizes_hab   <- csizes_hab[!is.na(csizes_hab)]
  csizes_inpa  <- csizes_inpa[!is.na(csizes_inpa)]
  csizes_outpa <- csizes_outpa[!is.na(csizes_outpa)]
  
  #compute area [km2] of all cells in geo_raster
  # hab_area   <- length(csizes_hab) * median(csizes_hab)
  hab_area   <- length(na.omit(values(hab))) * res(hab)[1]*res(hab)[2]
  # inpa_area  <- length(csizes_inpa) * median(csizes_inpa)
  inpa_area  <- length(na.omit(values(inpas))) * res(inpas)[1]*res(inpas)[2]
  # outpa_area <- length(csizes_outpa) * median(csizes_outpa)
  outpa_area <- length(na.omit(values(outpas))) * res(outpas)[1]*res(outpas)[2]

  
  areas <- data.frame(
    type = c("total", "in_pa", "out_pa"),
    area_sqkm = c(hab_area/1000000, inpa_area/1000000, outpa_area/1000000)
  )
  
  areas <- areas %>% mutate(perc = area_sqkm/(hab_area/1000000)*100)
  areas
  
  ## SAVE ##
  filename <- paste0("data/analysis/summaries/", datatype, "/depth_areas_maxdep", abs(maxdepth_rndnd), ".csv")
  fwrite(areas, filename)
}
