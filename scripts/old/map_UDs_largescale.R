## Create maps of tracking locations and KDEs zoomed in on each region ##
pacman::p_load(raster, mapview, sf, track2KBA, dplyr, ggplot2, stringr, rasterVis)

period <- "foraging"

## Load data ## 
## UDs 
# 1 is arithmetic 2 is weighted
percUDs <- readRDS(
  list.files(
    "data/analysis/UDs/CDF", full.names=T)[str_detect(list.files("data/analysis/UDs/CDF"), pattern=period)][1]) 
# percUDs <- readRDS("data/analysis/UDs/w_group_percUDs_h2.08_foraging.rds") # weighted mean
# percUDs <- readRDS("data/analysis/UDs/a_group_percUDs_h2.08_foraging.rds") # arithmetic mean

## Load bathymetry data ~~~~~~~~~~~
depth <- raster("data/geodata/ETOPO1.tiff")

## Load land data ~~~~~~~~~~~~~~~~~
# land <- rnaturalearth::ne_download(
#   scale=10, category = "cultural", returnclass = "sf")
land <- raster::shapefile("data/geodata/WA_adminboundaries/wca_admbnda_adm0_ocha_18022021.shp")
land <- st_as_sf(land)

### load MPA polygons ##----------------
folder <- "data/geodata/WA_conservation area_shp/"
files <- str_subset(list.files("data/geodata/WA_conservation area_shp/", full.names = T),  pattern = fixed(".shp"))

## load MPA shapefiles ## 
shp_list <- lapply(files, function(x) shapefile(x) ) # list of shps

## bind all areas of conservation interest from bijagos to mauritania ##
cons_areas <- bind(shp_list) %>% st_as_sf()
# mapview(cons_areas)

# separate biosphere reserve #
mpas <- subset(cons_areas, NAME != "Bijagos Archipelago Biosphere Reserve")
babr <- subset(cons_areas, NAME == "Bijagos Archipelago Biosphere Reserve")

## Tracking data ~~~~~~~~~~~~~~~~~~
folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods

tracks <- readRDS(paste0(folder, period, ".rds"))
tracks <- formatFields(tracks,
                       fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date")

# convert %UDs to polygons ~~~~~~~~~~~~~~~~~~~~~~~~~
UD95 <- percUDs$`95`
UD50 <- percUDs$`50`

UD95p <- raster::rasterToPolygons(UD95) %>% st_as_sf() %>% st_union() %>% st_transform(crs=4326)
UD50p <- raster::rasterToPolygons(UD50) %>% st_as_sf() %>% st_union() %>% st_transform(crs=4326)

UDs <- c(UD95p, UD50p) %>% st_sf() %>% mutate(UD=factor(c("95", "50"), levels = c("95", "50")))
# mapview::mapview(UD95p) + mapview::mapview(UD50p, col.regions="red")

## location of Poilão
poilao <- data.frame(label="Poilão", "Longitude" = -15.726667, "Latitude" = 10.864722)
poilao <- st_as_sf(poilao, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

## Plot one region at a time ## ~~~~~~~~~~~~~~~~~~~~
regions <- unique(tracks$destination[tracks$destination!="unknown"])

for(i in 1:length(regions)){
  print(i)
  one <- regions[i]
  
  regtrcks <- tracks %>% filter(destination == one) %>% st_as_sf(coords = c("Longitude", "Latitude"), 
                                                                 crs = 4326, agr = "constant")
  xtnt <- st_bbox(regtrcks)
  if(one == "Bijagos"){
    # xtnt[1] <- xtnt[1]-.15 ## original
    xtnt[1] <- xtnt[1]-.25
    xtnt[3] <- xtnt[3]+.1
    # c(xtnt[1]-.23, xtnt[3]+.1), ylim = c(xtnt[2], xtnt[4])
    w <- 9
    h <- 8
    depth_crp <- crop(depth, extent(xtnt[1]-.22, xtnt[3]+.1, xtnt[2],xtnt[4]))
  } else if(one=="Senegal-Gambia"){
    xtnt[1] <- xtnt[1]-.23
    xtnt[3] <- xtnt[3]+.1
    # xtnt[2] <- xtnt[2] - .15
    xtnt[4] <- xtnt[4] + .15
    w <- 7
    h <- 10
    depth_crp <- crop(depth, extent(xtnt[1]-.23, xtnt[3], xtnt[2],xtnt[4]))
  } else if(one=="Mauritania"){
    xtnt[1] <- xtnt[1]-.23
    xtnt[3] <- xtnt[3]+.1
    xtnt[4] <- xtnt[4] + .1
    w <- 7
    h <- 10
    depth_crp <- crop(depth, extent(xtnt[1]-.23, xtnt[3], xtnt[2],xtnt[4]))
  }
  
  # regUD95 <- st_crop(UD95p, xtnt)
  # regUD50 <- st_crop(UD50p, xtnt)
  
  cntr_lvls <- c(-25) ## depth contours to display
  cntrs <- st_as_sf(rasterToContour(depth_crp, levels=cntr_lvls))
  
  values(depth_crp) <- ifelse(values(depth_crp)>0, 0, values(depth_crp))
  # values(depth_crp) <- abs(values(depth_crp))
  
  # breaks <- c(0,-25,-50,-100,-200,-500)
  # breaks <- c(0,-25,-50,-100,-200,-500)
  # breaks <- c(0,25,50,100,200,500)
  # breaks <- c(500,200,100,50,25,0)
  breaks <- c(-1500,-500,-100,-25,-15,0)
  
  # map <- ggplot() +    # no depth 
  map <- gplot(depth_crp) +  # w/ depth
    geom_tile(aes(fill = value)) +
    scale_fill_fermenter(breaks=breaks, direction=-1, limits = c(-2000,0), name = "Depth (m)") +
    ggnewscale::new_scale_fill() + 
    geom_sf(data=UDs, aes(fill=UD), inherit.aes = FALSE, color=NA) +
    geom_sf(data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
    geom_sf(
      data=mpas, inherit.aes = FALSE, aes(linetype="solid"), color="black", size=1, fill=NA, show.legend = "line") +
    {if(one == "Bijagos")
      geom_sf(
        data=babr, inherit.aes = FALSE, aes(linetype="dashed"), color="black",  size=1, fill=NA, show.legend = "line")} +
    geom_sf(data = regtrcks, color="black", alpha=0.05, size=0.25, inherit.aes = FALSE) +
    {if(one == "Bijagos")
      geom_sf(
        data = poilao, inherit.aes = FALSE, fill="gold1", color="black", size = 5, stroke=1.5, shape=23)
    } +
    geom_sf_label(data = poilao, aes(label = label), inherit.aes = FALSE, nudge_x = .15) +
    # geom_sf(data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
    coord_sf(
      xlim = c(xtnt[1], xtnt[3]), ylim = c(xtnt[2], xtnt[4]), expand = F) +
    # {if(one == "Mauritania")
    #   coord_sf(
    #     xlim = c(xtnt[1]-.23, xtnt[3]+.1), ylim = c(xtnt[2], xtnt[4]), expand = F)} +
    # {if( one %in% c("Bijagos", "Senegal-Gambia") )
    #   coord_sf(
    #     xlim = c(xtnt[1]-.23, xtnt[3]+.1), ylim = c(xtnt[2], xtnt[4]), expand = F)} +
    scale_linetype_identity(
      name = NULL,
      breaks = c("solid", "dashed"),
      labels = c("MPA", "BABR"),
      guide = "legend") +
    theme(
      legend.position = c(1.125, .5),
      legend.background = element_rect(fill=NA, colour=NA),
      panel.background=element_rect(fill="white", colour="black"),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.text=element_text(size=12, color="black"),
      axis.title=element_text(size=16)) +
    {if(one != "Bijagos")
      theme(legend.position = c("none"))} +
    {if(one == "Bijagos")
      guides(
        linetype = guide_legend(keywidth = unit(2,"line")),
        fill = guide_legend(override.aes = list(linetype = 0))
      )} +
    ylab("") +
    xlab("") +
    ggspatial::annotation_scale(style="ticks", height=unit(0.5, "cm"), text_cex = 1)
  # map
  
  ## Save ##
  ggsave(paste0("figures/largescale_UDs_", one, "X.png"), plot=map, width = w, height=h )
  
  if(one == "Bijagos"){ p1 <- map } else if(one == "Senegal-Gambia"){p2 <- map} else {p3 <- map}
}


## patch together maps ## 
library(patchwork)

layout <- 
"
AAAAAA
AAAAAA
AAAAAA
AAAAAA
BBBCCC
BBBCCC
BBBCCC
BBBCCC
"

maps <- p1 + p2 + p3 +
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 28), plot.tag.position = c(.21, .94))
# & theme(plot.tag = element_text(size = 30), plot.tag.position = c(.925, .90))
ggsave( "C:/Users/Martim Bill/Desktop/test/1111111111111.png", plot=maps, width = 11, height=13 )


