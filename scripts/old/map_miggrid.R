### map migration route density ### 

pacman::p_load(dplyr, sf, ggplot2, sp, stringr, rasterVis)

## Migration route density grid ##
migrid <- readRDS("data/analysis/mig_grid/mig_grid_10x10_n22.rds")
prj <- proj4string(migrid)

## location of Poilão
poilao <- data.frame("Longitude" = -15.726667, "Latitude" = 10.864722) 
poilao <- poilao %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>% 
  st_transform(crs=prj)

## land polyons ## 
land <- raster::shapefile("data/geodata/WA_adminboundaries/wca_admbnda_adm0_ocha_18022021.shp")
land <- st_as_sf(land) %>% st_transform(crs=prj)

## load MPA shapefiles ## 
folder <- "data/geodata/shp/"
files <- str_subset(list.files("data/geodata/WA_conservation area_shp/", full.names = T),  pattern = fixed(".shp"))

shp_list <- lapply(files, function(x) shapefile(x) ) # list of shps

## bind all areas of conservation interest from bijagos to mauritania ##
cons_areas <- bind(shp_list) %>% st_as_sf() %>% st_transform(crs = prj)

# separate thebiosphere reserve
mpas <- subset(cons_areas, NAME != "Bijagos Archipelago Biosphere Reserve")
babr <- subset(cons_areas, NAME == "Bijagos Archipelago Biosphere Reserve")

# xtnt <- extent(migrid)
xtnt <- extent(raster::projectExtent(raster("data/geodata/ETOPO1.tiff"), crs=prj))


## Tracking Data input for marking migration end points ## ~~~~~~~~~~~~~~~~~~
folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods

tracks <- readRDS(paste0(folder, "migration", ".rds"))
tracks <- track2KBA::formatFields(tracks,
                                  fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date")

## filter to only turtles leaving Bijagos ##
tracks <- tracks %>% filter(destination != "Bijagos")

## filter out points w/in Bijos ## 
tracks <- tracks %>% filter(Latitude > 11.6)

lastpnts <- tracks %>% group_by(ID) %>% summarise(
  Longitude=last(Longitude), Latitude=last(Latitude))

n_distinct(lastpnts$ID) # one individuals drops out (182454)

# project tracks to data-centered projection
# TD <- track2KBA::projectTracks(tracks, projType = "azim", custom=T) # equal area azimuthal proj
TD <- track2KBA::projectTracks(lastpnts, projType = "azim", custom=T)  %>% 
  st_as_sf() %>% mutate(end_point = "End point")


## map it ## -----------------------------------------------------------------
map <- gplot(migrid) + 
  geom_tile(aes(fill = value), na.rm = T) +
  scale_fill_viridis_c(option="B", direction=-1, na.value = NA) +
  geom_sf(data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
  geom_sf(
    data=babr, inherit.aes = FALSE, aes(linetype="dashed"), color="black",  size=1, fill=NA, show.legend = "line") +
  geom_sf(
    data=mpas, inherit.aes = FALSE, aes(linetype="solid"), color="black", size=1, fill=NA, show.legend = "line") +
  geom_sf(
    data = poilao, inherit.aes = FALSE, fill="gold1", aes(color="Poilão"), size = 5.5, stroke=1.5, shape=23) +
  # ggnewscale::new_scale_fill() + 
  geom_sf(
    # data = TD, inherit.aes = FALSE, fill="red", color="red", size = 3.5, alpha=0.75, stroke=1.5, shape=3) +
    data = TD, inherit.aes = FALSE, color="red", size = 3.5, alpha=0.75, stroke=1.5, shape=3) +
  coord_sf(
    xlim = c(xtnt@xmin, xtnt@xmax), ylim = c(xtnt@ymin, xtnt@ymax), expand = F) +
  scale_x_continuous(breaks = seq(-19, -15, by = 1)) + 
  scale_color_manual(
    name = NULL,
    values = c("Poilão"="black")) +
  scale_linetype_identity(
    name = NULL,
    breaks = c("solid", "dashed"),
    labels = c("MPA", "BABR"),
    guide = "legend") +
  theme(
    legend.title = element_text(size=11),
    legend.position = c(0.85, .5),
    legend.background = element_rect(fill=NA, colour=NA),
    legend.key = element_rect(fill=NA, colour=NA),
    # legend.key.width = unit(2.5,"line"),
    panel.background=element_rect(fill="steelblue1", colour="black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text=element_text(size=12, color="black"),
    axis.title=element_text(size=16)) +
  ylab("") +
  xlab("") + 
  labs(fill = "N IDs") + guides(
    linetype = guide_legend(keywidth = unit(3,"lines"))) +
  ggspatial::annotation_scale(style="ticks", height=unit(0.5, "cm"), text_cex = 1)

## save ##
ggsave("figures/migrid_10x10_n22XxXxX.png", plot=map, width = 5, height=10 )

