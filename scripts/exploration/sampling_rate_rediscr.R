### Filter out individuals with <10? points ###

n_distinct(tracks_f$ID)
ids_wpnts <- tracks_f %>% 
  group_by(ID) %>% summarise(n_pnts=n()) %>% filter(n_pnts > 9)

pnts_filt <- tracks_f %>% filter(ID %in% ids_wpnts$ID)
n_distinct(tracks_f$ID) - n_distinct(pnts_filt$ID) # number of inds removed

## summarise sampling rate ##
library(amt)

## convert to amt 'tracks_xyt' format ##
tracks_amt <- pnts_filt %>% make_track(.x=Longitude, .y=Latitude, DateTime, id = ID,
                        crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_d"))

## calculate time step between each location ##
tracks_amt <- tracks_amt %>% nest(
  data = c(x_, y_, t_)) %>% mutate(
    ts = map(data, function(x){
  # summarize_sampling_rate(x, time_unit="day", summarize=F)
  c(NA, summarize_sampling_rate(x, time_unit="day", summarize=F))
})) %>% 
  dplyr::select(ts) %>% unnest(cols = ts) %>% bind_cols(tracks_amt)

## identify big gaps ## 
# how often was the time interval more than a opecified amount of time? 
gap_threshold <- 2 # let's say, 5 days
tracks_amt$over_thresh <- tracks_amt$ts > gap_threshold

table(tracks_amt$over_thresh)

## designate unique ids for each 'burst' ##
ids <- unique(tracks_amt$id)

tracks_amt <- rbindlist(
  lapply(seq_along(ids), function(x){
    one <- tracks_amt[tracks_amt$id == ids[x], ]
    # make ids for periods of non-gappy data
    one$over_thresh <- ifelse(is.na(one$over_thresh), F, one$over_thresh)
    one$burst <- paste(one$id,  cumsum(one$over_thresh), sep="_")
    return(one)
  })
)

str(tracks_amt)

## designate unique ids for each 'burst' ##

## sub-sample/interpolate within bursts ## 



## summarise sampling rate per individual ##
sampl_summ <- tracks_amt %>% nest(data = c(x_, y_, t_)) %>% mutate(ts = map(data, function(x){
  summarize_sampling_rate(x, time_unit="day")
})) %>% 
  select(id, ts) %>% unnest(cols = ts)

# distribution of average sampling rates across individuals - any outliers?
hist(sampl_summ$mean, 20)
hist(sampl_summ$median, 20)

# which individuals (if any) are outliers? (e.g. average sampling rate greater than 1 day)
outlie_ids <- sampl_summ$id[sampl_summ$median > 1]

# remove individuals with outlying sampling rates (too far from others for interpolationg/subsampling)
pnts_filt <- pnts_filt %>% filter(!ID %in% outlie_ids)

# how many and which ids were filtered out?
n_distinct(pnts_filt$ID)
unique(tracks_f$ID)[!unique(tracks_f$ID) %in% unique(pnts_filt$ID)]

# find suitable average sampling rate for remaining individuals
sampl_summ <- sampl_summ %>% filter(id %in% unique(pnts_filt$ID)) 

sampl_summ %>% summarise(
  unit = "hour",
  mean_mean = mean(mean)*24,
  median_median = median(median)*24
)

