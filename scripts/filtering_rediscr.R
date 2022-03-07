## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Second round of processing steps of green turtle data before analysis ##
# i.e., filtering and re-discretizing data #

pacman::p_load(track2KBA,
  data.table, amt, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter, adehabitatLT, mapview)

tracksum <- read.csv("data/tracks/tracking_summary.csv", stringsAsFactors = F) # deployment summary 

## Data input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
folder <- "data/analysis" # repository w/ datasets split into periods

periods <- c("internest", "migration", "foraging")

# how many satellites minimum to retain data? (4 or 6)
satnum_filter <- 6

pnts_id_list <- list()
n_locs_list <- list()

for(i in seq_along(periods)){
  # analyze foraging or migration?
  period <- periods[i]
  print(period)
  ## load datasets split into migration and foraging periods (post-migration) ##
  tracks <- readRDS(list.files(folder, pattern=fixed(period), full.names = T))
  
  ## rename some columns ## 
  tracks <- tracks %>% 
    rename(Argos.Location.Quality = "Location.Quality", GPS.HDOP=HDOP) %>% 
    mutate(
      sensor = ifelse(Argos.Location.Quality == "", "GPS", "Argos doppler shift")
    ) 
  
  n_locs_raw <- tracks %>% group_by(ID, sensor) %>% summarise(
    period = period,
    n_locs_raw = n()
  )
  
  ## Remove PTT points of lowest quality  ~~~~~~~~~~~~~~~~~~~~~
  table(tracks$Argos.Location.Quality)
  b4nrow  <- nrow(tracks)
  pnts_id <- tracks %>% group_by(ID) %>% summarise(n_pnts_r=n(), stage=period)
  tracks  <- tracks[!tracks$Argos.Location.Quality %in% c("Z","B","A"),]
  rmvd    <- 1-(nrow(tracks)/b4nrow)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Speed and inner angle filter for GPS data with SDLfilter ~~~~~~~~~~~~~~~~~~~~
  gps <- tracks[tracks$sensor == "GPS", ]
  
  # 'qi' is the location quality index for use in the SDLfilter functions (satellite number or HDOP)
  # gps <- gps %>% rename(id = ID, lat=Latitude, lon=Longitude) %>% mutate(qi=GPS.HDOP)
  gps <- gps %>% rename(id = ID, lat=Latitude, lon=Longitude) %>% mutate(qi=SatNum)

  # vmaxlp fxn requires at least 10 pnts per ID (remove any ids w/ fewer)
  # fewpnt_ids <- names(table(gps$id)[table(gps$id)<10])
  # if( length(fewpnt_ids )>0){
  #   gps_fltrd <- filter(gps, !id %in% fewpnt_ids)
  # }
  ## SDLfilter steps ## 
  vmax    <- vmax(gps, qi=6, prob=.95)    # maximum linear speed for filter
  vmaxlp  <- vmaxlp(gps, qi=5, prob=0.99) # maximum 'loop-speed'
  
  if(period == "migration"){
    gps_fltrd <- ddfilter(gps, vmax=vmax, vmaxlp = 8 , qi=satnum_filter, method=1)
  } else {
    gps_fltrd <- ddfilter(gps, vmax=vmax, vmaxlp = vmaxlp , qi=satnum_filter, ia=10, method=1)
  }
  
  ## Re-organize ~~
  gps_fltrd <- formatFields(
    gps_fltrd, fieldID = "id", 
    fieldLat="lat", fieldLon="lon", 
    fieldDateTime="DateTime") %>% dplyr::select(-qi) 
  
  ## Speed filter for PTT data with sdafilter ~~~~~~~~~~~~~~~~~~~~
  ptt <- tracks[tracks$sensor != "GPS", ]
  
  ptt_fltrd <- data.table::rbindlist(
    lapply(seq_along(unique(ptt$ID)), function(x){
      print(x)
      one <- subset(ptt, ID == unique(ptt$ID)[x])
      if(nrow(one)<5){
        return(one)
      } else {
        result <- argosfilter::sdafilter(
          lat   = one$Latitude,
          lon   = one$Longitude,
          dtime = one$DateTime,
          lc    = one$Argos.Location.Quality,
          vmax = vmax
        )
        
        one <- one[which(!result == "removed"),]
        return(one)
      }
    })
  )
  
  ## Recombine filtered GPS data with PTT data
  tracks_f <- ptt_fltrd %>% 
    full_join(gps_fltrd) %>% arrange(ID, DateTime)
  
  ## N satellite filter ##
  tracks_f <- filter(tracks_f, SatNum >= satnum_filter | is.na(SatNum))
  
  ## summarise filtering 
  n_locs <- tracks_f %>% group_by(ID, sensor) %>% summarise(
    period = period,
    n_locs_filter = n()
  ) %>% right_join(n_locs_raw) %>% 
    tidyr::pivot_wider(
      names_from = sensor, values_from = c(n_locs_filter, n_locs_raw)
  )
  
  colnames(n_locs)[3:6] <- c("n_filter_argos", "n_filter_gps", "n_raw_argos", "n_raw_gps")
  
  n_locs_list[[i]] <- n_locs

  ## check which points were filtered ##
  # tracks_sf <- st_as_sf(tracks, coords = c("Longitude", "Latitude"),
  #                       crs = 4326, agr = "constant")
  # tracks_f_sf <- st_as_sf(tracks_f, coords = c("Longitude", "Latitude"),
  #                         crs = 4326, agr = "constant")
  # mapview::mapview(tracks_sf) + mapview::mapview(tracks_f_sf, col.regions="red") # map data
  
  # left   <- 1-(nrow(tracks_f)/b4nrow)
  
  ### Filter out individuals with <x? points (pre-interpolation) ### ------------
  n_distinct(tracks_f$ID)   # total n
  ids_wpnts <- tracks_f %>% 
    group_by(ID) %>% summarise(n_pnts=n()) %>% filter(n_pnts > 9)
  n_distinct(ids_wpnts$ID) # n meeting threshold (5 is bare minimum for KDE)
  
  if(period != "migration"){
    tracks_f <- tracks_f %>% filter(ID %in% ids_wpnts$ID) # filter out small n-inds
  }
  
  pnts_id <- tracks_f %>% group_by(ID) %>% summarise(n_pnts_f=n(), stage=period) %>% left_join(pnts_id)
  
  
  ## summarise sampling rate ##--------------------------------------------------
  ## convert to amt 'tracks_xyt' format ##
  tracks_amt <- tracks_f %>% 
    make_track(.x=Longitude, .y=Latitude, .t=DateTime, 
               id = ID, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"))
  
  ## calculate time step (in days) between each location ##
  tracks_amt <- do.call(cbind.data.frame,
                        tracks_amt %>% 
    nest(data = c(x_, y_, t_)) %>% 
    mutate( ts = map(data, function(x){
      c(NA, summarize_sampling_rate(x, time_unit="day", summarize=F))
    })  ) %>% 
    unnest(cols = c(ts, data))
  )
  
  ts_summ <- summary(tracks_amt$ts) # summ stats of time step intervals
  ts_summ * 24 # in hours
  
  ## identify big gaps ## 
  # how often was the time interval more than a specified amount of time? 
  gap_threshold <- 1 # let's say, x days
  tracks_amt$over_thresh <- tracks_amt$ts > gap_threshold
  tracks_amt$over_thresh <- ifelse(is.na(tracks_amt$ts), TRUE, tracks_amt$over_thresh)
  
  table(tracks_amt$over_thresh) # TRUE is number of steps that are 'big gap'
  
  ## designate unique ids for each 'burst' (i.e. non gap period) ##
  ids <- unique(tracks_amt$id)
  
  tracks_amt <- tracks_amt %>% arrange(id, t_)
  
  tracks_amt <- rbindlist(
    lapply(seq_along(ids), function(x){
      one <- tracks_amt[tracks_amt$id == ids[x], ]
      # make ids for periods of non-gappy data
      one$over_thresh <- ifelse(is.na(one$over_thresh), F, one$over_thresh)
      # if(any(is.na(one$id))) print(x)
      one$burst <- paste(one$id,  cumsum(one$over_thresh), sep="_")
      return(one)
    })
  ) %>% dplyr::select(x_, y_, t_, id, ts, burst)
  
  str(tracks_amt)
  
  ## Interpolate ## 
  # convert to 'ltraj' object
  tracks_lt <- as.ltraj(tracks_amt, date=tracks_amt$t_, id=tracks_amt$burst, burst=tracks_amt$burst, typeII = T) 
  
  h <- 3 # how many hours should there be between interpolated points?
  tr <- redisltraj(tracks_lt, 3600*h, type="time", nnew = 2) # re-discretize
  
  # convert back to data frame
  id <- stringr::word(unlist(lapply(tr, function(x) attr(x, "id"))), 1, sep = "\\_")
  burst <- unlist(lapply(tr, function(x) attr(x, "burst")))
  tracks_re <- rbindlist(lapply(seq_along(tr), function(x) {
    one <- tr[[x]]
    one$id <- id[x]
    one$burst <- burst[x]
    return(one)
  }))
  
  pnts_id <- tracks_re %>% group_by(id) %>% summarise(n_pnts_i=n(), stage=period) %>% left_join(pnts_id, by=c("id"="ID", "stage"))
  
  ## convert to sf 
  # tracks_amt_sf <- st_as_sf(tracks_amt, coords = c("x_", "y_"),
  #                           crs = 4326, agr = "constant")
  # tracks_re_sf <- st_as_sf(tracks_re, coords = c("x", "y"),
  #                          crs = 4326, agr = "constant")
  
  # Compare real data to interpolated data for a single individual ##
  # anid <- unique(tracks_re$id)[5]
  # #
  # mapview(filter(tracks_re_sf, id == anid)) +
  # mapview(filter(tracks_amt_sf, id == anid), col.regions="red")
  
  ## add foraging destination classification to tracks ## ---------------------
  tracksum$PTT <- as.character(tracksum$PTT)
  tracks_re <- tracks_re %>% left_join(tracksum[,c("PTT", "destination")], by=c("id"="PTT"))
  
  tracksum$PTT <- as.character(tracksum$PTT)
  tracks_f <- tracks_f %>% left_join(tracksum[,c("PTT", "destination")], by=c("ID"="PTT"))
  
  ## Save ## --------------------------------------------------------------------
  
  ## raw, filtered data ##
  filename <- paste0(folder, "/raw_filtered/", period, "_", "satfilt", satnum_filter, ".rds")
  filename
  saveRDS(tracks_f, filename)
  
  ## interpolated data ##
  filename <- paste0(folder, "/interpolated/", period, "_", "satfilt", satnum_filter,".rds")
  filename
  saveRDS(tracks_re, filename)
  
  ## summarise filtering ## 
  
  # left_internest <- left
  # left_migration <- left
  # left_forage    <- left
  # 
  # first_forage <- 
  # 
  # c(left_internest, left_migration, left_forage)
  
  pnts_id_list[[i]] <- pnts_id
  
  #-------------------------------------------------------------------------------
  ## create dataset w/ only GPS data for those IDs that have it and PTT for rest ##
  gpsids <- unique(gps_fltrd$ID)
  
  ## Recombine filtered GPS data with PTT data, then exclude PTT data for IDs w/ GPS data
  tracks2 <- tracks_f %>% 
    filter(sensor == "Argos doppler shift") %>% 
    full_join(gps_fltrd) %>% arrange(ID, DateTime) %>% 
    filter((sensor == "GPS") | (!ID %in% gpsids))
  
  tracksum$PTT <- as.character(tracksum$PTT)
  tracks2 <- tracks2 %>% left_join(tracksum[,c("PTT", "destination")], by=c("ID"="PTT"))
  
  ## raw, filtered data (only GPS for IDs from 2019/20) ##
  filename <- paste0("data/sensitivity_analysis/raw_filtered_gpsids/", period, "_", "satfilt", satnum_filter, ".rds")
  filename
  saveRDS(tracks2, filename)
}

n_locs_df_period <- rbindlist(n_locs_list) %>% 
  mutate(
    n_filter_argos = ifelse(period == "migration", n_raw_argos, n_filter_argos),
    n_filter_gps = ifelse(period == "migration", n_raw_gps, n_filter_gps)
  ) 

n_locs_df <- n_locs_df_period %>% group_by(ID) %>% 
  summarise(
    n_raw_argos = sum(na.omit(n_raw_argos)),
    n_raw_gps = sum(na.omit(n_raw_gps)),
    n_filter_argos = sum(na.omit(n_filter_argos)),
    n_filter_gps = sum(na.omit(n_filter_gps))
  )

## Save 
fwrite(n_locs_df_period, "data/analysis/summaries/filtering_summary_by_period.csv")
fwrite(n_locs_df, "data/analysis/summaries/filtering_summary.csv")


# pnts_id_internest <- pnts_id
# pnts_id_migrate   <- pnts_id
# pnts_id_forage    <- pnts_id

## combine points per ID for each stage ##

pnts_id_all_s <- bind_rows(pnts_id_list[[1]], pnts_id_list[[2]], pnts_id_list[[3]]) %>% 
  mutate(id = as.factor(id)) %>% group_by(id) %>% 
  summarise(
    n_pnts_r = sum( n_pnts_r ),
    n_pnts_f = sum( n_pnts_f ),
    n_pnts_i = sum( n_pnts_i )
  )

sum(pnts_id_all_s$n_pnts_f)
sum(pnts_id_all_s$n_pnts_r)
sum(pnts_id_all_s$n_pnts_f)/sum(pnts_id_all_s$n_pnts_r)*100 # % pnts retained thru filters

## Save ## 
# pnts summary by individual

fwrite(pnts_id_all_s, paste0("data/analysis/summaries/n_pnts_byid", "_", "satfilt", satnum_filter,".csv"))
