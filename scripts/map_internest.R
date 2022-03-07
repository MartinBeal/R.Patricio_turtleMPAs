## Create maps of tracking locations and KDEs zoomed in on each region ##
pacman::p_load(raster, mapview, sf, track2KBA, dplyr, ggplot2, stringr, rasterVis)

datatypes <- c("raw_filtered","interpolated")
period <- "internest"
satfilt <- "satfilt6"

## Double loop, first through datatypes, then throug regions ##
for(y in seq_along(datatypes)){
  datatype <- datatypes[y]
  
  ## load tracking and UD (CDF) datasets ## 
  if(datatype == "interpolated"){
    folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods
    tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
    tracks <- formatFields(tracks,
                           fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date") 
    folder <- "data/analysis/UDs/interpolated/CDF"
    percUDs <- readRDS(
      list.files(
        folder, full.names=T)[str_detect(
          list.files(folder), pattern=fixed(period))][1])   
  } else if(datatype == "raw_filtered"){
    folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods
    tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
    folder <- "data/analysis/UDs/raw_filtered/CDF"
    percUDs <- readRDS(
      list.files(
        folder, full.names=T)[str_detect(
          list.files(folder), pattern=fixed(period))][1]) 
  }
  
  
  ## Load bathymetry data ~~~~~~~~~~~
  # depth <- raster("data/geodata/ETOPO1.tiff")
  depth <- raster("data/geodata/ETOPO1_patchedw_bathyrastPNBA.tif")
  
  ## Load land data ~~~~~~~~~~~~~~~~~
  # land <- rnaturalearth::ne_download(
  #   scale=10, category = "cultural", returnclass = "sf")
  land <- raster::shapefile("data/geodata/WA_adminboundaries/wca_admbnda_adm0_ocha_18022021.shp")
  land <- st_as_sf(land)
  
  ### load MPA polygons ##----------------
  # mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021_dissolve.shp")
  mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_Nov2021_noBoba.shp")
  babr <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/BABR_polygon.shp")
  
  mpas <- st_as_sf(mpas)
  babr <- st_as_sf(babr)
  
  # convert %UDs to polygons ~~~~~~~~~~~~~~~~~~~~~~~~~
  UD95 <- percUDs$`95`
  UD50 <- percUDs$`50`
  
  UD95p <- raster::rasterToPolygons(UD95) %>% st_as_sf() %>% st_union() %>% st_transform(crs=4326)
  UD50p <- raster::rasterToPolygons(UD50) %>% st_as_sf() %>% st_union() %>% st_transform(crs=4326)
  
  UDs <- c(UD95p, UD50p) %>% st_sf() %>% mutate(UD=factor(c("95", "50"), levels = c("95", "50")))
  # mapview::mapview(UD95p) + mapview::mapview(UD50p, col.regions="red")
  
  ## location of Poilão
  poilao_meio <- data.frame(label=c("Poilão", "Meio"), 
                            "Longitude" = c(-15.726667, -15.666024), 
                            "Latitude" = c(10.864722, 10.976397)
  )
  poilao_meio <- st_as_sf(poilao_meio, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")
  
  ## Plot one region at a time ## ~~~~~~~~~~~~~~~~~~~~
  regions <- unique(tracks$destination[tracks$destination!="unknown"])
  
  # xtnt <- st_bbox(tracks)
  xtnt <- st_bbox(c(xmin = -16.51626, xmax = -15.02416, ymax = 11.96899, ymin = 10.5), crs = st_crs(4326)) #w/out inset
  xtnt2 <- st_bbox(
    c(xmin = -15.93, xmax = -15.49, ymax = 11.11, ymin = 10.75), 
    crs = st_crs(4326)) # if w/ inset
  
  regtrcks <- tracks %>% st_as_sf(coords = c("Longitude", "Latitude"), 
                                  crs = 4326, agr = "constant") # spatialize TD
  
  # if(one == "Bijagos"){ 
  depth_crp <- crop(depth, extent(xtnt[1]-.25, xtnt[3]+.1, xtnt[2],xtnt[4]))
  # } else if(one=="Senegal-Gambia"){
  # xtnt[2] <- xtnt[2] - .15
  #   xtnt[4] <- xtnt[4] + .1
  #   w <- 7
  #   h <- 10
  #   depth_crp <- crop(depth, extent(xtnt[1]-.22, xtnt[3], xtnt[2],xtnt[4]))
  # } else if(one=="Mauritania"){
  #   w <- 7
  #   h <- 10
  #   depth_crp <- crop(depth, extent(xtnt[1]-.22, xtnt[3], xtnt[2],xtnt[4]))
  # }
  
  # regUD95 <- st_crop(UD95p, xtnt)
  # regUD50 <- st_crop(UD50p, xtnt)
  
  cntr_lvls <- c(-25) ## depth contours to display
  cntrs <- st_as_sf(rasterToContour(depth_crp, levels=cntr_lvls))
  
  values(depth_crp) <- ifelse(values(depth_crp)>0, 0, values(depth_crp))
  # values(depth_crp) <- abs(values(depth_crp))
  
  breaks <- c(-1500,-500,-100,-25,-15,0)
  
  ## bbox polygon for plotting (if with inset)
  crnrs <- rbind(c(xtnt2$xmin,xtnt2$ymin), c(xtnt2$xmax,xtnt2$ymin), 
                 c(xtnt2$xmax,xtnt2$ymax), c(xtnt2$xmin,xtnt2$ymax), c(xtnt2$xmin,xtnt2$ymin))
  bbox_poly <- st_polygon(list(crnrs)) %>% sf::st_sfc( crs=4326) 
  
  # map <- ggplot() +    # no depth 
  map <- gplot(depth_crp) +  # w/ depth
    geom_tile(aes(fill = value)) +
    scale_fill_fermenter(breaks=breaks, direction=-1, limits = c(-2000,0), name = "Depth (m)") +
    ggnewscale::new_scale_fill() + 
    geom_sf(data=UDs, aes(fill=UD), inherit.aes = FALSE, color=NA) +
    geom_sf(
      data = regtrcks, color="black", alpha=0.025, size=0.25, inherit.aes = FALSE) +
    geom_sf(
      data=bbox_poly, inherit.aes = FALSE, color="red3", size=1, fill=NA) + ## inset bbox polygon
    geom_sf(
      data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
    geom_sf(
      data=mpas, inherit.aes = FALSE, aes(linetype="solid"), color="black", size=1, fill=NA, show.legend = "line") +
    geom_sf(
      data=babr, inherit.aes = FALSE, aes(linetype="dashed"), color="black",  size=1, fill=NA, show.legend = "line") +
    geom_sf(
      data = poilao_meio, inherit.aes = FALSE, fill="gold1", color="black", size = 2.5, stroke=1.5, shape=23) +
    # geom_sf_label(
    #   data = poilao_meio, aes(label = label), inherit.aes = FALSE, nudge_x = .15) +
    coord_sf(
      xlim = c(xtnt[1]-.25, xtnt[3]+.1), ylim = c(xtnt[2], xtnt[4]), expand = F) + 
    scale_y_continuous(breaks = seq(10.6, 11.8, by = .4)) +
    scale_linetype_identity(
      name = NULL,
      breaks = c("solid", "dashed"),
      labels = c("MPA", "BABR"),
      guide = "legend") +
    theme(
      # legend.position = c(1.125, .5),
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.background = element_rect(fill=NA, colour=NA),
      panel.background=element_rect(fill="white", colour="black"),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      legend.text = element_text(size=16, color="black"),
      legend.title = element_text(size=16, color="black"),
      axis.text = element_text(size=16, color="black"),
      axis.title=element_blank()
    ) +
    guides(
      linetype = guide_legend(keywidth = unit(2,"line")),
      fill = guide_legend(override.aes = list(linetype = 0))
    ) +
    ggspatial::annotation_scale(
      style="ticks", height=unit(0.5, "cm"), text_cex = 1.35)
  map
  
  ## Save ##
  # filename <- paste0("figures/internest/internest_", datatype, "X.png")
  # ggsave(filename, plot=map, width = 9, height=7 )
  
  ### Create inset zoomed in on Poilao ###
  
  inmap <- gplot(depth_crp) +  # w/ depth
    geom_tile(aes(fill = value)) +
    scale_fill_fermenter(breaks=breaks, direction=-1, limits = c(-2000,0), name = "Depth (m)") +
    ggnewscale::new_scale_fill() + 
    geom_sf(data=UDs, aes(fill=UD), inherit.aes = FALSE, color=NA) +
    geom_sf(
      data=mpas, inherit.aes = FALSE, aes(linetype="solid"), color="black", size=1, fill=NA, show.legend = "line") +
    geom_sf(
      data=babr, inherit.aes = FALSE, aes(linetype="dashed"), color="black",  size=1, fill=NA, show.legend = "line") +
    geom_sf(data = regtrcks, color="black", alpha=0.025, size=0.25, inherit.aes = FALSE) +
    geom_sf(data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
    geom_sf(
      data = poilao_meio, inherit.aes = FALSE, fill="gold1", color="black", size = 2.5, stroke=1.5, shape=23) +
    geom_sf_label(data = poilao_meio, aes(label = label), inherit.aes = FALSE, nudge_x = .05, alpha=.65) +
    coord_sf(
      xlim = c(xtnt2[1], xtnt2[3]), ylim = c(xtnt2[2], xtnt2[4]), expand = F) +
    scale_linetype_identity(
      name = NULL,
      breaks = c("solid", "dashed"),
      labels = c("MPA", "BABR"),
      guide = "legend") +
    theme(
      legend.position = c("none"),
      plot.margin = unit(c(0,0,0,0), "cm"),
      # legend.background = element_rect(fill=NA, colour=NA),
      panel.background=element_rect(fill=NA, colour=NA),
      plot.background=element_rect(fill=NA, colour=NA),
      panel.border = element_rect(colour = "red3", fill=NA, size=2),
      # axis.text=element_text(size=16, color="black"),
      axis.text=element_blank(),
      axis.title=element_blank(),
      axis.ticks=element_blank()
    ) +
    # guides(
    #   linetype = guide_legend(keywidth = unit(2,"line")),
    #   fill = guide_legend(override.aes = list(linetype = 0))
    # ) +
    ggspatial::annotation_scale(
      style="ticks", 
      height=unit(0.5, "cm"), 
      text_cex = 1.35)
  
  # inmap
  
  ## Combine maps ##
  library(patchwork)
  
  # comb <- map + inset_element(inmap,  0.7, 0, 1, .5, align_to = 'panel')
  # comb <- map + inset_element(inmap,  0, 0, .3, .5, align_to = 'panel')
  # comb <- map + inset_element(inmap,   0.6, 0.6, 1, 1, align_to = 'plot')
  # comb$patches$layout$widths  <- 1
  # comb$patches$layout$heights <- 1
  # comb
  
  comb <- map + inmap + 
    plot_layout(guides="collect") + 
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 28))
  
  # save #
  filename <- paste0("figures/internest/UDs_internest_", datatype, "_sideby_XX.png")
  ggsave(plot=comb, filename, width = 12, height=5)  
}