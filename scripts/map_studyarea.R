## make map of study area, with bathymetry and MPA borders ## 

pacman::p_load(ggplot2, sf, stringr, dplyr, raster, rasterVis)

### load MPA polygons ##----------------
# mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp")
mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021.shp")
babr <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/BABR_polygon.shp")

mpas <- st_as_sf(mpas)
babr <- st_as_sf(babr)

## load bathymetry data ##
# depth <- raster("data/geodata/GEBCO_2020_W Africa/gebco_2020_n20.9_s10.1_w-18.9_e-14.8.tif")
# depth <- raster("data/geodata/ETOPO1.tiff")
depth     <- raster("data/geodata/ETOPO1_patchedw_bathyrastPNBA.tif")

# nodepthzone <- raster::shapefile("data/geodata/Bathy_BA/no_depth_info_zone.shp")
# nodepthzone_sf  <- st_as_sf(nodepthzone)

# plot(depth)
# plot(nodepthzone, add=T)

## political borders/coast
land <- raster::shapefile("data/geodata/WA_adminboundaries/wca_admbnda_adm0_ocha_18022021.shp")
land <- st_as_sf(land)

## location of Poilão
poilao <- data.frame("Longitude" = -15.726667, "Latitude" = 10.864722)
poilao <- st_as_sf(poilao, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

## create depth raster where cells falling within unsampled (untrustworthy) zone are changed to NA ##
# depth2 <- mask(depth, nodepthzone, inverse=T)

xtnt <- extent(depth)

# plot(depth2)
# cntr_lvls <- c(-24)
cntr_lvls <- c(-13)
# cntr_lvls <- c(-50, -100, -500, -1000, -2000, -3000 -4000, -5000)
cntrs <- st_as_sf(rasterToContour(depth, levels=cntr_lvls))

gplot(depth) + 
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(direction=-1) +
  # geom_sf(data=nodepthzone_sf, inherit.aes = FALSE, fill="grey90", color=NA) +
  geom_sf(data=cntrs, inherit.aes = FALSE, color="grey20", size=.35, fill=NA) +
  # geom_sf(data=babr, inherit.aes = FALSE, aes(color="grey1"), linetype="twodash",  size=.75, fill=NA, show.legend = T) +
  # geom_sf(data=mpas, inherit.aes = FALSE, aes(color="black"), size=1.25, fill=NA) +
  geom_sf(
    data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
  geom_sf(
    data=babr, inherit.aes = FALSE, aes(linetype="dashed"), color="black",  size=1, fill=NA, show.legend = "line") +
  geom_sf(
    data=mpas, inherit.aes = FALSE, aes(linetype="solid"), color="black", size=1, fill=NA, show.legend = "line") +
  geom_sf(
    data = poilao, inherit.aes = FALSE, fill="gold1", color="black", 
    size = 5.5, stroke=1.5, shape=23,
    show.legend = FALSE) +
  coord_sf(
    xlim = c(xtnt@xmin, xtnt@xmax), ylim = c(xtnt@ymin, xtnt@ymax), expand = F) +
  scale_x_continuous(breaks = seq(-19, -15, by = 1)) +
  # scale_color_manual(
  #   name = NULL,
  #   values = c("Poilão"="black")) +
  scale_linetype_identity(
    name = NULL,
    breaks = c("solid", "dashed"),
    labels = c("MPA", "BABR"),
    guide = "legend") +
  theme(
    legend.title = element_text(size=16),
    legend.text  = element_text(size=16),
    legend.position = c(0.82, .5),
    legend.background = element_rect(fill=NA, colour=NA),
    legend.key = element_rect(fill=NA, colour=NA),
    # legend.key.width = unit(2.5,"line"),
    panel.background=element_rect(fill=NA, colour="black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text=element_text(size=16, color="black")) +
  ylab("") +
  xlab("") + 
  labs(fill = "Depth (m)") + guides(linetype = guide_legend(keywidth = unit(2.5,"line"))) +
  ggspatial::annotation_scale(style="ticks", height=unit(0.5, "cm"), text_cex = 1)

## SAVE ## 
ggsave("figures/study_area/smallscale_WAfrica_13m_2xx.png", width = 5, height=12 )

