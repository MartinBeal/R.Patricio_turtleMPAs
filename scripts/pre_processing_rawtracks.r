## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Pre-processing steps for green turtle PTT/Fast-loc GPS data ## 

pacman::p_load(track2KBA, ctmm, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter, data.table)

## Metadata ##
tracksum <- read.csv("data/tracks/tracking_summary.csv", stringsAsFactors = F) # deployment summary 

## W Africa bounding box (for filtering) ## 
bbox <- raster::shapefile("C:/Users/Martim Bill/Documents/mIBA_package/data/green_turtles/west_africa_bbox.shp")
bbox@proj4string <- CRS(SRS_string = "EPSG:4326")

## load raw tracks ## 
files <- list.files("data/tracks/combined_tracks", full.names = T)

rawdata <- rbindlist(
  lapply(seq_along(files), function(x) {
    print(x)
    one <- fread(files[x]) 
    
    if("[ErrorRadius (m), SemiMajor (m), SemiMinor (m), Orientation ]" %in% colnames(one)){
      xx <- str_split(
        str_remove(one$`[ErrorRadius (m), SemiMajor (m), SemiMinor (m), Orientation ]`, 
                   pattern = "]") %>% 
          str_remove(pattern = fixed("[")), ", ")
      
      one$argos_error_radius <- as.numeric(do.call(rbind, lapply(xx, function(x) x[1])))
    } else {
      one$argos_error_radius <- rep(NA) 
    }

    one <- one %>% 
      dplyr::select(
        Tag_ID, Latitude, Longitude, 
        UTC_Date, UTC_Time, 
        "Location Quality", HDOP, "SatNum", "argos_error_radius") %>% 
      mutate(
        UTC_Date = as.Date(UTC_Date, tryFormats = c("%Y-%m-%d", "%d/%m/%Y")))
    
    return(one)
  })
)

## Identifying stationary vs. migratory periods ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tracks <- formatFields(
  rawdata, 
  fieldID = "Tag_ID", 
  fieldLat="Latitude", fieldLon="Longitude", 
  fieldDate="UTC_Date", fieldTime="UTC_Time", formatDT = "ymd HMS")

tracks$ID <- do.call(rbind,  stringr::str_split(tracks$ID, pattern = stringr::fixed(":")))[,2] # remove colon and number before it 

tracks <- filter(tracks, Latitude != 0 | !is.na(Latitude))

# remove points from NZ (calibration)
tracksSP <- SpatialPointsDataFrame(
  SpatialPoints(data.frame(tracks$Longitude, tracks$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")),
  data=tracks)
tracksSP <- tracksSP[bbox,] ## Keep only points in bbox (a polygon around West Africa)

# mapview::mapview(tracksSP)

tracks <- data.frame(tracksSP) %>% dplyr::select(-tracks.Longitude, -tracks.Latitude)
tracks <- tracks[order(tracks$ID, tracks$DateTime), ] ## order date time stamps within each individual

# remove Leatherbacks and Hawksbills
tracks <- tracks[!tracks$ID %in% c("197219","197220","182460","60861","60862","60895","197139","197140","197141"), ]
# remove data sets with no relevant information for IN, MI, or FG periods (e.g. PTT data 2001)
tracks <- tracks[!tracks$ID %in% c("15240","15241","15242","15243","15244"), ]

# "60861","60862","60895" # IDs to remove?

# remove data sets with no post-migration data
# tracks <- tracks[!tracks$ID %in% c("182455","182456","60900","60866","60888", "182461"), ]

# Remove rows with duplicate time stamps #
tracks_list <- vector(mode="list", dplyr::n_distinct(tracks$ID))
for(i in 1:dplyr::n_distinct(tracks$ID)){
  one <- tracks[tracks$ID==unique(tracks$ID)[i], ]
  one$duplicated <- duplicated(one$DateTime)
  tracks_list[[i]] <- one
}
tracks <- do.call(rbind, tracks_list)

tracks <- tracks[tracks$duplicated == FALSE, ] 


#-------------------------------------------------------------------------------
## Split data into inter-nesting, migratory, and post-migratory periods
# i.e. compare dates in tracking data to dates in metadata spreadsheet
ID_list <- list()
internest_list <- list()
forage_list    <- list()
migrate_list   <- list()

## split tracks into stages based on dates in tracking summary table ##
for(i in 1:n_distinct(tracks$ID)){
  print(i)
  one <- tracks[tracks$ID == unique(tracks$ID)[i], ]
  ID <- unique(tracks$ID)[i]
  dest <- tracksum[tracksum$PTT==ID,]$destination
    
  start.mig <- tryCatch(
    expr = {
      as.Date(tracksum[tracksum$PTT==ID,]$date.leaving.Poilão, tryFormats = c("%d/%m/%Y"))
    },
    error = function(e) NA
  )
  end.mig <- tryCatch(
    expr = {
      as.Date(tracksum[tracksum$PTT==ID,]$first.date.at.FG, tryFormats = c("%d/%m/%Y", "%Y-%m-%d"))    },
    error = function(e) NA
  )
  
  if(is.na(start.mig) & is.na(end.mig)){ # only inter-nesting data
    internest <- one
    migration <- NA
    foraging  <- NA
  } else if(is.na(end.mig)){  # no foraging data 
    internest <- one[  one$DateTime < start.mig, ]
    migration <- one[ (one$DateTime > start.mig & one$DateTime < end.mig), ]
    foraging  <- NA
  } else if(is.na(start.mig) & (!is.na(end.mig))){ # no migration data 
    internest <- one[  one$DateTime < start.mig, ]
    migration <- NA
    foraging  <- one[ (one$DateTime > end.mig), ]
  } else {  # data from all periods 
    internest <- one[  one$DateTime < start.mig, ]
    migration <- one[ (one$DateTime > start.mig & one$DateTime < end.mig), ]
    foraging  <- one[ (one$DateTime > end.mig), ]
  }
  
  if(dest %in% c("Bijagos", "unknown")){ # if no mig. data or destination is Bijagos ignore 'migration'
    migration <- NA
  }
  
  ID_list[[i]]        <- ID
  internest_list[[i]] <- internest
  if(!is.na(foraging)) {
    forage_list[[i]]  <- foraging
    if(any(is.na(foraging$Latitude))) print("foraging NAs!!")
  }
  if(!is.na(migration)) {
    migrate_list[[i]] <- migration
    if(any(is.na(migration$Latitude))) print("foraging NAs!!")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # if(is.na(start.mig)) next
  # end.mig <- as.Date(tracksum[tracksum$PTT==ID,]$first.date.at.FG, tryFormats = c("%d/%m/%Y", "%Y-%m-%d"))
  # 
  # if(!is.na(start.mig)){
  #   internest <- one[  one$DateTime < start.mig, ]
  #   migration <- one[ (one$DateTime > start.mig & one$DateTime < end.mig), ]
  #   if(any(is.na(migration$Latitude))) print("migration NAs!!")
  # } else {
  #   internest <- one
  #   migration <- NA
  #   }
  # 
  # foraging  <- one[ (one$DateTime > end.mig), ]
  # if(any(is.na(foraging$Latitude))) print("foraging NAs!!")
  # 
  # ID_list[[i]]        <- ID 
  # internest_list[[i]] <- internest
  # forage_list[[i]]    <- foraging
  # migrate_list[[i]]   <- migration

}

# summarise number of points of post-migr. and migration for each individ.
# summary <- data.frame(
#   ID=do.call(rbind, lapply(ID_list, function(x) x[[1]])), 
#   n_internest=do.call(rbind, lapply(internest_list, function(x) nrow(x))),
#   n_migrate=do.call(rbind, lapply(migrate_list, function(x) nrow(x))),
#   n_forage=do.call(rbind, lapply(forage_list, function(x) nrow(x)))
# )

tracks_i <- do.call(rbind, internest_list) %>% filter(!is.na(Latitude)) # internesting period
tracks_m <- do.call(rbind, migrate_list)   %>% filter(!is.na(Latitude)) # migration period
tracks_f <- do.call(rbind, forage_list )   %>% filter(!is.na(Latitude)) # foraing period 

# convert to spatial data to visualize
# fSP <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks_f$Longitude, tracks_f$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data=tracks_f) # foraging
# mSP <- SpatialPointsDataFrame(SpatialPoints(data.frame(tracks_m$Longitude, tracks_m$Latitude), proj4string=CRS("+proj=longlat +datum=WGS84")), data=tracks_m) # migration
# visualize in interactive map
# mapview::mapview(fSP)
# mapview::mapview(mSP)

# raw data sample sizes for each period
n_i <- n_distinct(tracks_i$ID) 
n_m <- n_distinct(tracks_m$ID) 
n_f <- n_distinct(tracks_f$ID)
c(n_i, n_m, n_f)

## SAVE ## --------------------------------------------------------------------
saveRDS(tracks_i,
        paste0("data/analysis/internest_unfiltered_GPS_PTT_", n_i, ".rds"))
saveRDS(tracks_f,
        paste0("data/analysis/foraging_unfiltered_GPS_PTT_", n_f, ".rds"))
saveRDS(tracks_m,
        paste0("data/analysis/migration_unfiltered_GPS_PTT_", n_m, ".rds"))
