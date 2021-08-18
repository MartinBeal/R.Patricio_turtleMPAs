## Calculate depths visited by turtles ## 

pacman::p_load(
  track2KBA, amt, dplyr, sp, sf, move, ggplot2, stringr, SDLfilter, raster, data.table)

## Data input ~~~~~~~~~~~~~~~~~~
folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods

# analyze foraging or migration?
# period <- "internest"
period <- "foraging"
# period <- "migration"

tracks <- readRDS(paste0(folder, period, ".rds"))

tracks <- formatFields(tracks,
                       fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date")

TD_ll <- st_as_sf(x = tracks, 
                        coords = c("Longitude", "Latitude"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>% as_Spatial()

## load bathymetry data ##
# depth <- raster("data/geodata/GEBCO_2020_W Africa/gebco_2020_n20.9_s10.1_w-18.9_e-14.8.tif")
depth <- raster("data/geodata/ETOPO1.tiff")

# nodepthzone <- raster::shapefile("data/geodata/Bathy_BA/no_depth_info_zone.shp")

# plot(depth)
# plot(nodepthzone, add=T)

## create depth raster where cells falling within unsampled (untrustworthy) zone are changed to NA ##
# depth2 <- mask(depth, nodepthzone, inverse=T)

## extract depth for each tracking point
TD_ll$depth <- extract(depth, TD_ll)
hist(TD_ll$depth)
# mapview(TD_ll)

## Load KDE UDs ##
# UDs <- readRDS("data/analysis/UDs/CDF/a_groupCDFs_h2.08_c1_foraging.rds")
# UD50 <- UDs[[1]]
# UD95 <- UDs[[2]]
UD50 <- readRDS(paste0("data/analysis/UDs/ind_areas/individual_UD50_h2.08_c1_foraging.rds"))
UD95 <- readRDS(paste0("data/analysis/UDs/ind_areas/individual_UD95_h2.08_c1_foraging.rds"))

## extract depth from cells falling within each individual's UD (both 50 and 95%) 
UD50_ll <- UD50$UDPolygons %>% as_Spatial() %>% spTransform(CRS("+init=epsg:4326"))
UD95_ll <- UD95$UDPolygons %>% as_Spatial() %>% spTransform(CRS("+init=epsg:4326"))

UD50_depths <- extract(depth, UD50_ll)
UD95_depths <- extract(depth, UD95_ll)

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
})) %>% bind_cols(UD50_ll@data)

# UD 95%
depths_95 <- rbindlist(lapply(UD95_depths, function(x){ 
  data.frame(
    min = min(x),
    mn  = mean(x),
    sd  = sd(x),
    md  = median(x)
  )
})) %>% bind_cols(UD95_ll@data)

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
by_ind <- TD_ll@data %>% group_by(ID) %>% 
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
# fwrite(crss_ind, "data/analysis/summaries/depths_used_avgfromTD.csv")
# fwrite(summ_50, "data/analysis/summaries/depths_used_avgfromUD50.csv")
# fwrite(summ_95, "data/analysis/summaries/depths_used_avgfromUD95.csv")


### Calculate the amount of available depth area inside and outside MPAs ### ---
folder <- "data/geodata/shp/"
files <- str_subset(list.files("data/geodata/WA_conservation area_shp/", full.names = T),  pattern = fixed(".shp"))

## load MPA shapefiles ## 
shp_list <- lapply(files, function(x) shapefile(x) ) # list of shps

## bind all areas of conservation interest from bijagos to mauritania ##
cons_areas <- bind(shp_list)
# mapview(cons_areas)

# filter out biosphere reserve for overlaying points
mpas <- subset(cons_areas, NAME != "Bijagos Archipelago Biosphere Reserve")
# maxdepth <- crss_ind$md_min # points
# maxdepth <- summ_95$min_min    # UD 95%
maxdepth <- summ_95$mn_min    # UD 95%


hab <- depth
hab <- setValues(hab, 
                 ifelse( (values(hab) < maxdepth)|(values(hab) >= 0), 
                        NA, values(hab)))
inpas  
outpas 

inpas <- raster::intersect(hab, mpas)
plot(inpas)

mpas_r <- rasterize(mpas, hab, field=rep(1, nrow(mpas)))
plot(mpas_r)

## habitat within PAs
inpas <- setValues(hab, 
                   ifelse( (values(mpas_r) == 0), 
                           NA, values(hab)))
plot(inpas)
## habitat outside PAs 
outpas <- setValues(hab, 
                    ifelse( !is.na(values(mpas_r) ), 
                            NA, values(hab)))
plot(outpas)
## calculate areas ## --------------------------------------------------------

## first, entire suitable depth habitat ## 
csizes_hab   <- raster::area(hab, na.rm=TRUE, weights=FALSE)
csizes_inpa  <- raster::area(inpas, na.rm=TRUE, weights=FALSE)
csizes_outpa <- raster::area(outpas, na.rm=TRUE, weights=FALSE)

#delete NAs from vector of all raster cells
csizes_hab   <- csizes_hab[!is.na(csizes_hab)]
csizes_inpa  <- csizes_inpa[!is.na(csizes_inpa)]
csizes_outpa <- csizes_outpa[!is.na(csizes_outpa)]

#compute area [km2] of all cells in geo_raster
hab_area   <- length(csizes_hab) * median(csizes_hab)
inpa_area  <- length(csizes_inpa) * median(csizes_inpa)
outpa_area <- length(csizes_outpa) * median(csizes_outpa)

areas <- data.frame(
  type = c("total", "in_pa", "out_pa"),
  area_sqkm = c(hab_area, inpa_area, outpa_area)
)

areas <- areas %>% mutate(perc = area_sqkm/hab_area*100)

## SAVE ##
fwrite(areas, "data/analysis/summaries/depth_areas_maxdep12.csv")
