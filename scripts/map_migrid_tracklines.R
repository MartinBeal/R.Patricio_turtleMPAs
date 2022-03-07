## make map of study area, with Tracking Data ## 

pacman::p_load(ggplot2, sf, stringr, dplyr, raster, rasterVis)

datatypes <- c("raw_filtered","interpolated")

satfilt <- "satfilt6"

datatype <- "raw_filtered"
folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods


## dbbmm grid ## 
migrid <- readRDS("data/analysis/mig_grid/rawdata_dbbmm_mig_grid_10x10_UD95_18_perc.rds")

prj <- proj4string(migrid)

## migration tracks ##
period <- "migration"
mig_tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
mig_tracks <- mig_tracks %>% filter(ID != "205279" & UTC_Time != "06:35:52")


## map ## ---------------------------------------------------------------------

tracks <- mig_tracks

if(datatype == "interpolated"){
  tracks <- tracks %>% rename(ID = id, Longitude=x, Latitude=y)
}

tracks <- tracks %>% st_as_sf(
  coords = c("Longitude", "Latitude"),
  crs = 4326, agr = "constant"
) %>% st_transform(crs=prj)
## create linestrings from points
track_lines <- tracks %>% 
  group_by(ID) %>% summarize(do_union=FALSE) %>% st_cast("LINESTRING")

## political borders/coast
land <- raster::shapefile("data/geodata/WA_adminboundaries/wca_admbnda_adm0_ocha_18022021.shp")
land <- st_as_sf(land) %>% 
  st_transform(crs=prj)

## location of Poilão and Meio
poilao_meio <- data.frame(label=c("Poilão", "Meio"), 
                          "Longitude" = c(-15.726667, -15.666024), 
                          "Latitude" = c(10.864722, 10.976397)
)
poilao_meio <- poilao_meio %>% st_as_sf(
  coords = c("Longitude", "Latitude"),
  crs = 4326, agr = "constant") %>% 
  st_transform(crs=prj)

xtnt <- extent(raster::projectExtent(depth, crs=prj))


gplot(migrid) +
  geom_tile(aes(fill = value), na.rm = T) +
  scale_fill_viridis_c(option="B", direction=-1, na.value = NA) +
  geom_sf(
    data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
  geom_sf(data = track_lines, 
          color="black", 
          show.legend = FALSE,
          size=0.5, inherit.aes = FALSE) +
  geom_sf(
    data = poilao_meio, inherit.aes = FALSE, fill="gold1", color="black", 
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
    panel.background=element_rect(fill="steelblue1", colour="black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text=element_text(size=16, color="black")) +
  ylab("") +
  xlab("") + 
  labs(fill = "% routes") + guides(linetype = guide_legend(keywidth = unit(2.5,"line"))) +
  ggspatial::annotation_scale(style="ticks", height=unit(0.5, "cm"), text_cex = 1)

## save
ggsave("figures/study_area/study_area_wMigrid_wTracks_rawfiltered.png", width = 5, height=12 )
