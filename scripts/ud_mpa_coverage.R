## Overlap UD areas w/ MPAs ## 

pacman::p_load(sp, raster, mapview, sf)

periods <- c("internest", "foraging")


# for()
period <- periods[i]

if(period == "foraging"){
  cdfs <- readRDS("data/analysis/UDs/raw_filtered/CDF/a_groupCDFs_h1.97_c1_foraging.rds")
} else {
  cdfs <- readRDS("data/analysis/UDs/raw_filtered/CDF/a_groupCDFs_h0.75_c1_internest.rds")
}

ud50 <- cdfs[[1]]
ud95 <- cdfs[[2]]

mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp")

mpas <- spTransform(mpas, ud95@crs)

## create empty raster for MPAs (based on UD raster dimensions)
r <- ud95
values(r) <- NA

mpas_rast <- rasterize(mpas, r)

## UD 95 - MPA coverage ## -----------------------------------------------------
# cells of ud in mpa
ud95mpa <- mask(ud95, mpas_rast)

pixArea <- res(ud95)[1]

ncells_ud95  <- sum(!is.na(raster::getValues(ud95[[1]])))
area_95  <- (ncells_ud95 * pixArea^2) / (1000^2)  # area of full range (sq.km)
ncells_ud95mpa  <- sum(!is.na(raster::getValues(ud95mpa[[1]])))
area_ud95mpa5  <- (ncells_ud95mpa * pixArea^2) / (1000^2)  # area of full range (sq.km)

## % MPA coverage of 95% UD area
hr95_over <- (area_ud95mpa5/area_95) * 100

## UD 50 - MPA coverage ## -----------------------------------------------------
# cells of ud in mpa
ud50mpa <- mask(ud50, mpas_rast)

pixArea <- res(ud50)[1]

ncells_ud50  <- sum(!is.na(raster::getValues(ud50[[1]])))
area_50  <- (ncells_ud50 * pixArea^2) / (1000^2)  # area of full range (sq.km)
ncells_ud50mpa  <- sum(!is.na(raster::getValues(ud50mpa[[1]])))
area_ud50mpa5  <- (ncells_ud50mpa * pixArea^2) / (1000^2)  # area of full range (sq.km)

## % MPA coverage of 95% UD area
hr50_over <- (area_ud50mpa5/area_50) * 100

hr_over_gen <- data.frame(
      period = period,
          ud = c("50", "95"),
  perc_cover = c(hr50_over, hr95_over)
)

## Split foraging regions ## ---------------------------------------------------

# 1. split CDF raster into three, by destination
# 2. calc. overlap one at a time

## region-specific bounding boxes

xtnt <- sf::st_as_sf(mpas) %>% st_transform(crs=4326) %>% st_bbox
if(one == "Bijagos"){
  xtnt[1] <- -16.81
  xtnt[2] <- 10.375
  xtnt[3] <- -14.9
  xtnt[4] <- 11.875

  
  
  depth_crp <- crop(depth, extent(xtnt[1]-.22, xtnt[3]+.1, xtnt[2],xtnt[4]))
  # axis labels

} else if(one=="Senegal-Gambia"){
  xtnt[1] <- -17.3
  xtnt[2] <- 13.18
  xtnt[3] <- -16.35
  xtnt[4] <- 14.4
  w <- 7
  h <- 10
  depth_crp <- crop(depth, extent(xtnt[1]-.23, xtnt[3], xtnt[2]-0.1, xtnt[4] + 0.1))
  # axis labels
  ytix <- seq(13.2, 14.4, by=0.4)
  xtix <- seq(-17.1, -16.5, by=0.2)
} else if(one=="Mauritania"){
  xtnt[1] <- -17.25
  xtnt[2] <- 19.3
  xtnt[3] <- -16.15
  xtnt[4] <- 20.71
  w <- 7
  h <- 10
  depth_crp <- crop(depth, extent(xtnt[1]-.23, xtnt[3], xtnt[2]-0.1, xtnt[4]+0.1))
  # axis labels
  ytix <- seq(19.4, 20.6, by=0.4)
  xtix <- seq(-17.1, -16.3, by=0.2)
}

