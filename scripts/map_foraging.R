## Create maps of tracking locations and KDEs zoomed in on each region ##
pacman::p_load(raster, mapview, sf, track2KBA, dplyr, ggplot2, stringr, rasterVis)

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
mpas <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/WDPA_MPAs_Wafrica_May2021.shp")
babr <- raster::shapefile("data/geodata/WDPA_MPAs_Wafrica_May2021/BABR_polygon.shp")

mpas <- st_as_sf(mpas) %>% filter(!NAME %in% c("Banc dÃ¢â‚¬â„¢Arguin", "Banc d'Arguin"))
babr <- st_as_sf(babr)

### set datatype loop ## --------------------------------------------------

datatypes <- c("raw_filtered","interpolated")
period <- "foraging"
satfilt <- "satfilt6"

## Double loop, first through datatypes, then through regions ##
for(y in seq_along(datatypes)){
  datatype <- datatypes[y]
  weighttype <- c("a_group", "w_group") # arithmetic or weighted
  
  for(z in 1:2){ # first unweighted then weighted
    print(paste("z =", z))
    ## load tracking and UD (CDF) datasets ## 
    if(datatype == "interpolated"){
      folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods
      tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
      tracks <- formatFields(tracks,
                             fieldID = "id", fieldLat="y", fieldLon="x", fieldDateTime="date") 
      folder <- "data/analysis/UDs/interpolated/CDF"
      files  <- list.files(
        folder, full.names=T, pattern=period)
      percUDs <- readRDS(
        files[str_detect(
          files, pattern=fixed(weighttype[z])
        )]
      ) # 1 = arithmetic, 2 = weighted
    } else if(datatype == "raw_filtered"){
      folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods
      tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
      folder <- "data/analysis/UDs/raw_filtered/CDF"
      files  <- list.files(
        folder, full.names=T, pattern=period)
      percUDs <- readRDS(
        files[str_detect(
          files, pattern=fixed(weighttype[z])
          )]
        ) # 1 = arithmetic, 2 = weighted
    }
    
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
    
    for(i in 1:length(regions)) {
      print(i)
      one <- regions[i]
      
      regtrcks <- tracks %>% filter(destination == one) %>% st_as_sf(coords = c("Longitude", "Latitude"), 
                                                                     crs = 4326, agr = "constant")
      xtnt <- st_bbox(regtrcks)
      if(one == "Bijagos"){
        # xtnt[1] <- xtnt[1]-.15 ## original
        # xtnt[1] <- xtnt[1]-.25
        xtnt[1] <- -16.81
        xtnt[2] <- 10.375
        # xtnt[3] <- xtnt[3]+.1
        xtnt[3] <- -14.9
        xtnt[4] <- 11.875
        # c(xtnt[1]-.23, xtnt[3]+.1), ylim = c(xtnt[2], xtnt[4])
        w <- 9
        h <- 8
        depth_crp <- crop(depth, extent(xtnt[1]-.22, xtnt[3]+.1, xtnt[2],xtnt[4]))
        # axis labels
        ytix <- seq(10.5, 11.7, by=0.4)
        xtix <- seq(-16.5, -15.0, by=0.5)
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
      
      # regUD95 <- st_crop(UD95p, xtnt)
      # regUD50 <- st_crop(UD50p, xtnt)
      
      cntr_lvls <- c(-12) ## depth contours to display
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
        geom_sf(data = UDs, aes(fill=UD), inherit.aes = FALSE, color=NA) +
        geom_sf(data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
        geom_sf(
          data=mpas, inherit.aes = FALSE, aes(linetype="solid"), 
          color="black", size=1, fill=NA, show.legend = "line") +
        {if(one == "Bijagos")
          geom_sf(
            data=babr, inherit.aes = FALSE, aes(linetype="dashed"), 
            color="black",  size=1, fill=NA, show.legend = "line")
            } +
        geom_sf(data = regtrcks, color="black", alpha=0.05, size=0.25, inherit.aes = FALSE) +
        {if(one == "Bijagos")
          geom_sf(
            data = poilao, inherit.aes = FALSE, 
            fill="gold1", color="black", size = 5, stroke=1.5, shape=23)
        } +
        scale_y_continuous(breaks = ytix) + # custom tick labels
        scale_x_continuous(breaks = xtix) + # custom tick labels
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
          legend.title = element_text(size=16, color="black"),
          legend.text = element_text(size=16, color="black"),
          legend.position = c(1.125, .5),
          legend.background = element_rect(fill=NA, colour=NA),
          panel.background=element_rect(fill="white", colour="black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text = element_text(size=16, color="black"),
          axis.title = element_blank()
          ) +
        {if(one != "Bijagos")
          theme(legend.position = c("none"))} +
        {if(one == "Bijagos")
          guides(
            linetype = guide_legend(keywidth = unit(2,"line")),
            fill = guide_legend(override.aes = list(linetype = 0))
          )} +
        ylab("") +
        xlab("") +
        ggspatial::annotation_scale(
          style="ticks", height=unit(0.5, "cm"), 
          line_width = 1.5, text_cex = 1.35)
      map
      
      ## Save ##
      if(z == 1){
        filename <- paste0("figures/foraging/largescale_UDs_", one, "_", datatype, "_arithmeticX.png")
        ggsave(filename, plot=map, width = w, height=h ) 
      } else if(z==2) {
        filename <- paste0("figures/foraging/largescale_UDs_", one, "_", datatype, "_weightedX.png")
        ggsave(filename, plot=map, width = w, height=h )
      }
      
      if(one == "Bijagos") {p1 <- map} else if(one == "Senegal-Gambia") {p2 <- map} else {p3 <- map}
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
      plot_annotation(tag_levels = 'A') & theme(
        plot.tag = element_text(size = 28), 
        plot.tag.position = c(.22, .94),
        legend.box.margin=margin(-10,-10,-10,-30),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length=unit(.15, "cm")
      )
  # & theme(plot.tag = element_text(size = 30), plot.tag.position = c(.925, .90))
    # ggsave("C:/Users/Martim Bill/Desktop/test/AA11.png", plot=maps, width = 11, height=13 )
    
    if(z == 1){
    filename <- paste0("figures/foraging/combined_forage_", datatype, "_arithmeticX.png")
    ggsave(filename, plot=maps, width = 11, height=13 )
    } else if(z == 2){
      filename <- paste0("figures/foraging/combined_forage_", datatype, "_weightedX.png")
      ggsave(filename, plot=maps, width = 11, height=13 )
    }
  }
 
}
