### map migration route density ### 

pacman::p_load(dplyr, sf, ggplot2, sp, stringr, rasterVis)

datatypes <- c("raw_filtered","interpolated")
# datatypes <- c("raw_filtered","interpolated", "sensitivity_analysis")

migrid <- readRDS("data/analysis/mig_grid/rawdata_dbbmm_mig_grid_10x10_UD99_18_perc.rds")

prj <- proj4string(migrid)

## location of Poilão
poilao <- data.frame("Longitude" = -15.726667, "Latitude" = 10.864722) 
poilao <- poilao %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>% 
  st_transform(crs=prj)

## land polyons ## 
land <- raster::shapefile("data/geodata/WA_adminboundaries/wca_admbnda_adm0_ocha_18022021.shp")
land <- st_as_sf(land) %>% st_transform(crs=prj)

### load MPA polygons ##----------------
mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp")
mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021.shp")
babr <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/BABR_polygon.shp")

mpas <- st_as_sf(mpas)
babr <- st_as_sf(babr)

# xtnt <- extent(migrid)
xtnt <- extent(raster::projectExtent(raster("data/geodata/ETOPO1.tiff"), crs=prj))

## Tracking Data input for marking migration end points ## ~~~~~~~~~~~~~~~~~~
# if(datatype == "interpolated"){
#   folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods
# } else if(datatype == "raw_filtered"){
#   folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods
# } else if(datatype == "sensitivity_analysis"){
#   # folder <- "data/sensitivity_analysis/raw_filtered_gpsids/" # data for sensitivity analysis (only GPS data for 2019/2020 data)
# }
# tracks <- readRDS(paste0(folder, "migration", "_", satfilt, ".rds"))
# 
# if("id" %in% colnames(tracks)){
#   tracks <- formatFields(tracks,
#                          fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date") 
# }
# 
# ## filter to only turtles leaving Bijagos ##
# tracks <- tracks %>% filter(destination != "Bijagos")
# 
# ## filter out points w/in Bijos ## 
# tracks <- tracks %>% filter(Latitude > 11.6)
# 
# lastpnts <- tracks %>% group_by(ID) %>% summarise(
#   Longitude=last(Longitude), Latitude=last(Latitude))
# 
# n_distinct(lastpnts$ID) # one individuals drops out (182454)
# 
# # project tracks to data-centered projection
# # TD <- track2KBA::projectTracks(tracks, projType = "azim", custom=T) # equal area azimuthal proj
# TD <- track2KBA::projectTracks(lastpnts, projType = "azim", custom=T)  %>% 
#   st_as_sf() %>% mutate(end_point = "End point")

migrid_fac <- migrid 

migrid_fac <- reclassify(migrid, c(0, 24.999, 1, 25, 65, 2))

## map it ## -----------------------------------------------------------------
map <- gplot(migrid) +
  geom_tile(aes(fill = value), na.rm = T) +
  scale_fill_viridis_c(option="B", direction=-1, na.value = NA) +
  # gplot(migrid_fac) +
  # geom_tile(aes(fill = factor(value)), na.rm = T) +
  # scale_fill_manual(values = viridis(n=2, option="B", begin = .2, end = 1, direction=-1 ))
  geom_sf(data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
  geom_sf(
    data=babr, inherit.aes = FALSE, aes(linetype="dashed"), color="black",  size=1, fill=NA, show.legend = "line") +
  geom_sf(
    data=mpas, inherit.aes = FALSE, aes(linetype="solid"), color="black", size=1, fill=NA, show.legend = "line") +
  geom_sf(
    data = poilao, inherit.aes = FALSE, fill="gold1", aes(color="Poilão"), size = 5.5, stroke=1.5, shape=23) +
  # ggnewscale::new_scale_fill() + 
  # geom_sf(
  #   data = TD, inherit.aes = FALSE, color="red", size = 3.5, alpha=0.75, stroke=1.5, shape=3) +
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
  labs(fill = "% routes") + guides(
    linetype = guide_legend(keywidth = unit(3,"lines"))) +
  ggspatial::annotation_scale(style="ticks", height=unit(0.5, "cm"), text_cex = 1)

## save ##
# ggsave(paste0("figures/migrid_10x10_n22XxXxX_", datatype, ".png"), 
#        plot=map, width = 5, height=10 )
ggsave(paste0("figures/migrid_10x10_n18_dbbmm_perc99.png"),
       plot=map, width = 5, height=10 )


### binary (hi / low)
migrid_fac <- migrid
values(migrid_fac) <- ifelse(values(migrid) > 25, 1, NA)

## map it ## -----------------------------------------------------------------
map <- gplot(migrid_fac) +
  geom_tile(aes(fill = factor(value)), na.rm = T) +
  scale_fill_manual(
    values = "red",
    na.value = NA) +
  geom_sf(data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
  geom_sf(
    data=babr, inherit.aes = FALSE, aes(linetype="dashed"), color="black",  size=1, fill=NA, show.legend = "line") +
  geom_sf(
    data=mpas, inherit.aes = FALSE, aes(linetype="solid"), color="black", size=1, fill=NA, show.legend = "line") +
  geom_sf(
    data = poilao, inherit.aes = FALSE, fill="gold1", aes(color="Poilão"), size = 5.5, stroke=1.5, shape=23) +
  # ggnewscale::new_scale_fill() + 
  # geom_sf(
  #   data = TD, inherit.aes = FALSE, color="red", size = 3.5, alpha=0.75, stroke=1.5, shape=3) +
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
  labs(fill = "% routes") + guides(
    linetype = guide_legend(keywidth = unit(3,"lines"))) +
  ggspatial::annotation_scale(style="ticks", height=unit(0.5, "cm"), text_cex = 1)

