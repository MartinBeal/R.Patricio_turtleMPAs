## Plotting tracking extents per individual and stage ##

pacman::p_load(track2KBA,
               data.table, amt, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter, adehabitatLT, mapview)

datatypes <- c("raw_filtered","interpolated")
satfilt <- "satfilt6"

for(y in seq_along(datatypes)){
  
  datatype <- datatypes[y]
  periods <- c("internest", "migration", "foraging")
  
  if(datatype=="raw_filtered"){
    # Raw data input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rawfolder <- "data/analysis/raw_filtered" # raw data after speed and quality filtering
    
    rawfiles <- list.files(rawfolder, pattern=satfilt, full.names = T)
    
    rawtracks_l <- lapply(seq_along(rawfiles), function(x){
      if(x==1){stage <- "internesting"} else if(x==2){
        stage <- "migration" } else if(x==3){stage <- "foraging"}
      one <-   readRDS(str_subset(rawfiles, pattern = fixed(periods[x])))
      one$stage <- rep(stage)
      return(one)
    })
    ## COMBINE STAGE TD per ID ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tracks <- rbindlist(rawtracks_l)
    tracks_amt <- tracks %>% ## raw data
      make_track(.x=Longitude, .y=Latitude, .t=DateTime,
                 id = ID, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"))
  } else if( datatype == "interpolated"){
    ## Interpolated data input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    intfolder <- "data/analysis/interpolated/" # interpolated data
    
    # choose whether to include sub-6 satellite data or naw
    satfilt <- "satfilt6"
    
    intfiles <- list.files(intfolder, pattern=satfilt, full.names = T)
    
    inttracks_l <- lapply(seq_along(intfiles), function(x){
      if(x==1){stage <- "internesting"} else if(x==2){
        stage <- "migration" } else if(x==3){stage <- "foraging"}
      one <- readRDS(str_subset(intfiles, pattern = fixed(periods[x])))
      one$stage <- rep(stage)
      return(one)
    })
    tracks <- rbindlist(inttracks_l)
    tracks_amt <- tracks %>% ## interpolated data
      make_track(.x=x, .y=y, .t=date,
                 id = id, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), all_cols = F)
  }
  
  ## COMBINE STAGE TD per ID ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## filter out individuals not appearing in other analyses ## 
  # i.e., old tracks which don't have enough high quality data to analyse
  # tracks <- tracks %>% dplyr::filter(!ID %in% c("15242", "15243", "182454"))
  # tracks <- tracks %>% dplyr::filter(!id %in% c("15242", "15243", "182454"))
  
  ## Identify gaps ##--------------------------------------------------
  ## convert to amt 'tracks_xyt' format ##
  # tracks_amt <- tracks %>% ## raw data
  #   make_track(.x=Longitude, .y=Latitude, .t=DateTime,
  #              id = ID, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"))
  # tracks_amt <- tracks %>% ## interpolated data
  #   make_track(.x=x, .y=y, .t=date,
  #              id = id, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"), all_cols = F)
  ## calculate time step (in days) between each location ##
  tracks_amt <- do.call(rbind, 
                  lapply(split(tracks_amt, tracks_amt$id), function(x){
                    xy <- x %>% mutate( 
                      ts =  c(NA, summarize_sampling_rate(x, time_unit="day", summarize=F))  
                    )
                    return(xy)
                  })
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
      over_thresh <- one$ts > gap_threshold
      over_thresh <- ifelse(is.na(over_thresh), F, over_thresh)
      doy  <- yday(one$t_)
      doy2 <- ifelse(doy>182, doy - 182, doy + 182)
      
      one <- one %>% 
        mutate(cross_year=as.numeric(doy2==1 & !duplicated(doy2)))
      
      addid <- cumsum(one$cross_year) + cumsum(over_thresh)
      # if(any(is.na(one$id))) print(x)
      one$burst <- paste(one$id,  addid, sep="_")
      return(one)
    })
  ) %>% dplyr::select(x_, y_, t_, id, ts, burst)
  
  str(tracks_amt)
  
  ## raw data ##
  if(datatype=="raw_filtered"){
  tracks_amt <- left_join(tracks_amt, tracks[, c("Longitude", "Latitude", "DateTime", "destination", "stage")], by=c("x_"="Longitude", "y_"="Latitude", "t_"="DateTime")) ## raw data
  } else if( datatype == "interpolated"){
  ## interpolated data ##
  tracks_amt <- left_join(tracks_amt, tracks[, c("x", "y", "date", "destination", "stage")], by=c("x_"="x", "y_"="y", "t_"="date"))
  }
  tracks <- formatFields(tracks_amt,
                         fieldID = "id", fieldLat="y_", fieldLon="x_", fieldDateTime="t_")
  
  # one <- tracks %>% filter(ID == unique(tracks$ID)[1])
  # one <- tracks %>% filter(ID == "182461")
  one <- tracks
  
  # one$gr <- with(one, paste0(variable,time))
  one$ID <- as.factor(one$ID)
  one$stage <- factor(one$stage, levels=c("internesting", "migration", "foraging"))
  one$doy <- yday(one$DateTime)
  one$doy2 <- ifelse(one$doy > 182, one$doy - 182, one$doy + 182)
  
  one <- one %>% group_by(ID) %>% arrange(doy2)
  
  ## filter out two bad (time) points for 60893
  one$year <- year(one$DateTime)
  oneid <- filter(one, ID == "60893" & year != "2019")
  one <- filter(one, ID != "60893")
  one <- one %>% bind_rows(oneid) %>% 
    arrange(ID, DateTime)
  ## and for 182461 
  oneid <- filter(one, ID == "182461" & year != "2020")
  one <- filter(one, ID != "182461")
  one <- one %>% bind_rows(oneid) %>% 
    arrange(ID, DateTime)
  ## and for 182461 
  oneid <- filter(one, ID == "205286" & year != "2021")
  one <- filter(one, ID != "205286")
  one <- one %>% bind_rows(oneid) %>% 
    arrange(ID, DateTime)
  
  # dates <- c('March', 'June', 'September', 'December') # for display on X axis
  dates <- c('September', 'December', 'March', 'June') # shift year
  
  ## Demarcatae start and ends
  startends <- one %>% group_by(ID, destination) %>% arrange(DateTime) %>% summarise(
    first_DT = first(DateTime),
    last_DT = last(DateTime),
    start = yday(first(DateTime)),
    end   = yday(last(DateTime))
    # start = first(doy),
    # end   = last(doy)
  ) %>% tidyr::pivot_longer(c(start, end), names_to = "type", values_to = "doy")
  
  ## shift yr to start on July 1st
  startends$doy2 <- ifelse(startends$doy>182, startends$doy - 182, startends$doy + 182)
  
  ## add point where there is only data from one day to make sure it shows up on plot
  onex <- do.call(rbind,
                  lapply(split(one, one$burst), function(x){
                    if(n_distinct(x$doy2)==1){
                      xx <- rbind(x[1,], x)
                      xx$doy2 <- c(x$doy2[1]-1, x$doy2)
                    } else {xx <- x}
                    return(xx)
                  })
  )
  
  
  trckxtnts <- ggplot() +
    geom_point(data=startends, aes(x=doy2, y=ID, fill=type), color="white", size=6, pch=21) +
    geom_line(data=onex, aes(x=doy2, y=ID, group=burst, color=stage), size=5) +
    scale_fill_manual(values = c("red", "blue"), guide = guide_legend(reverse = TRUE)) +
    xlim(c(1,365)) + 
    theme_bw() + xlab("") +
    scale_x_continuous(
      breaks=c(62, 153, 242, 335),
      minor_breaks = c(1, 31, 92, 123, 183, 215, 273, 303), labels=dates,
      expand = c(0, 0))
  ## split by destination 
  trckxtnts_facet <- trckxtnts + 
    facet_wrap(~destination, scales = "free_y") 
  trckxtnts_facet
  
  ## SAVE ## ---------------------------------------
  filename <- paste0("figures/ind_track_extents/ind_track_extents_", satfilt, "_", datatype, "X.png")
  ggsave(filename, plot=trckxtnts, width=10, height=7)
  # ggsave(paste0("figures/ind_track_extents_", satfilt, "_", "rawX.png"), plot=trckxtnts, width=10, height=7)
  
  filename <- paste0("figures/ind_track_extents/ind_track_extents_", satfilt, "_", datatype, "_facettedX.png")
  ggsave(filename, plot=trckxtnts_facet, width=10, height=7)
  # ggsave(paste0("figures/ind_track_extents_", satfilt, "_", "raw_facettedX.png"), plot=trckxtnts_facet, width=10, height=7)

}
