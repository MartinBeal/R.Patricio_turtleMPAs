## Prep 4 PTT tracks from 2001-02 ##

files <- list.files("data/tracks/fourmoreturtlesgoingtopnba", full.names = T)

rawdata <- rbindlist(
  lapply(seq_along(files), function(x) {
    one <- fread(files[x], check.names = T)
    one <- one %>% mutate(
        HDOP = rep(NA),
        SatNum = rep(NA)
      ) %>% rename(
        Tag_ID=tag_id,"Location Quality"=lc, Latitude=lat1, Longitude=lon1
      )
    
    one$UTC_Date <- do.call(rbind, str_split(one$utc, pattern=" "))[,1]
    one$UTC_Time <- do.call(rbind, str_split(one$utc, pattern=" "))[,2]
    
    one <- dplyr::select(one,
      Tag_ID, Latitude, Longitude, UTC_Date, UTC_Time, "Location Quality", HDOP, "SatNum")
    
    return(one)
  })
)

rawdata$Tag_ID <- str_remove(rawdata$Tag_ID, pattern = fixed("a"))

fwrite(rawdata, "data/tracks/combined_tracks/green_turtles_PTT_2001_2002.csv")

