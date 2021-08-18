## MAPS ## =================================================================
#===========================================================================

pacman::p_load(stringr, ggplot2, raster, sf)

folder <- "data/analysis/UDs"

UD50 <- readRDS(list.files(folder, pattern="UD50", full.names = T)[1])
UD95 <- readRDS(list.files(folder, pattern="UD95", full.names = T)[1])

UD50 <- UD50$UDPolygons
UD95 <- UD95$UDPolygons
## quick maps
# mapKDE(UD50) 
# mapKDE(UD95)

land <- rnaturalearth::ne_download(
  scale=10, category = "cultural", returnclass = "sf")

## Metadata ##
tracksum <- read.csv("data/tracks/tracking_summary.csv", stringsAsFactors = F) # deployment summary 

### load MPA polygons ##----------------
folder <- "data/geodata/shp/"
files <- str_subset(list.files("data/geodata/shp/", full.names = T),  pattern = fixed(".shp"))

## load MPA shapefiles ## 
shp_list <- lapply(files, function(x) shapefile(x) ) # list of shps

## bind all areas of conservation interest from bijagos to mauritania ##
cons_areas <- bind(shp_list) %>% st_as_sf()
# mapview(cons_areas)

# filter out biosphere reserve for overlaying points
mpas <- subset(cons_areas, NAME != "Bijagos Archipelago Biosphere Reserve")
babr <- subset(cons_areas, NAME == "Bijagos Archipelago Biosphere Reserve")

# custom map of individual ranges # 
coordsets <- st_bbox(UD95)
# coordsets[1] <- coordsets[1] - .85
# coordsets[3] <- coordsets[1] + 1

## union together individual %UD areas 
lvl50 <- st_union(UD50)
lvl95 <- st_union(UD95)

## location of Poilão
poilao <- data.frame("Longitude" = -15.726667, "Latitude" = 10.864722)
poilao <- st_as_sf(poilao, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")


UDPLOT <- ggplot() +
  # geom_sf(data=lvl95, color="orangered", size=1, fill=NA) +
  geom_sf(data=lvl95, color="orangered", fill="orangered") +
  geom_sf(data=lvl50, color="orangered4", fill="orangered4") +
  # geom_sf(data=lvl50, color="orangered3", size=1, fill="orangered3") +
  # geom_sf(data=lvl95, color="orangered2", size=1, fill="orangered") +
  geom_sf(data=mpas, color="black", size=1, fill=NA) +
  geom_sf(data = poilao, inherit.aes = FALSE, color="black", fill="gold1", size = 5.5, stroke=1.5, shape=23) +
  coord_sf(xlim = c(coordsets$xmin - 1.75, coordsets$xmax + .1), ylim = c(coordsets$ymin, coordsets$ymax), expand = TRUE) +
  borders("world",fill="grey65", colour="grey40") +
  theme(
    panel.background=element_rect(fill="white", colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text=element_text(size=12, color="black"),
    axis.title=element_text(size=16),
    legend.position = "none") +
  ylab("") +
  xlab("")
UDPLOT

ggsave("figures/foraging_UDs_MPAs.png", UDPLOT, width = 6, height=12 )


## separate by 'destination' ## ----------------------------------------------

# add destination info to UD datasets
tracksum$PTT <- as.character(tracksum$PTT)
UD50 <- UD50 %>% left_join(tracksum[,c("PTT", "destination")], by=c("id"="PTT"))
UD95 <- UD95 %>% left_join(tracksum[,c("PTT", "destination")], by=c("id"="PTT"))

split(UD50, UD50$destination)

lvl50_g <- UD50 %>% group_by(destination) %>% st_union(UD50)
lvl95_g <- UD95 %>% group_by(destination) %>% st_union(UD95)


