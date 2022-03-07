## make map of study area, with Tracking Data ## 

pacman::p_load(ggplot2, sf, stringr, dplyr, raster, rasterVis)

datatypes <- c("raw_filtered","interpolated")

satfilt <- "satfilt6"

datatype <- "raw_filtered"
folder <- "data/analysis/raw_filtered/" # repository w/ datasets split into periods

# datatype <- "interpolated"
# folder <- "data/analysis/interpolated/" # repository w/ datasets split into periods

## internesting
period <- "internest"
int_tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
int_tracks <- int_tracks %>% filter(ID != "182455" & UTC_Time != "18:00:03")
## migration
period <- "migration"
mig_tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))
mig_tracks <- mig_tracks %>% filter(ID != "205279" & UTC_Time != "06:35:52")

## foraging
period <- "foraging"
for_tracks <- readRDS(paste0(folder, period, "_", satfilt, ".rds"))

## map ## ---------------------------------------------------------------------

tracks <- mig_tracks
tracks <- rbind(int_tracks, mig_tracks, for_tracks) %>% arrange(ID, DateTime)

if(datatype == "interpolated"){
  tracks <- tracks %>% rename(ID = id, Longitude=x, Latitude=y)
}

tracks <- tracks %>% st_as_sf(
  coords = c("Longitude", "Latitude"),
  crs = 4326, agr = "constant"
)
## create linestrings from points
track_lines <- tracks %>% 
  group_by(ID) %>% summarize(do_union=FALSE) %>% st_cast("LINESTRING")

## load bathymetry data ##
depth     <- raster("data/geodata/ETOPO1_patchedw_bathyrastPNBA.tif")

## political borders/coast
land <- raster::shapefile("data/geodata/WA_adminboundaries/wca_admbnda_adm0_ocha_18022021.shp")
land <- st_as_sf(land)

## location of Poilão
poilao <- data.frame("Longitude" = -15.726667, "Latitude" = 10.864722)
poilao <- st_as_sf(poilao, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

xtnt <- extent(depth)

# plot(depth2)
# cntr_lvls <- c(-24)
# cntr_lvls <- c(-13)
cntr_lvls <- c(-15)
# cntr_lvls <- c(-50, -100, -500, -1000, -2000, -3000 -4000, -5000)
cntrs <- st_as_sf(rasterToContour(depth, levels=cntr_lvls))

gplot(depth) + 
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(direction=-1) +
  geom_sf(data=cntrs, inherit.aes = FALSE, color="grey20", size=.35, fill=NA) +
  geom_sf(
    data = land, inherit.aes = FALSE, fill="grey65", colour="grey40") +
  geom_sf(data = track_lines, 
          color="black", 
          show.legend = FALSE,
          size=0.5, inherit.aes = FALSE) +
  # {if(datatype == "raw_filtered")
  #   geom_sf(data = tracks, color="black", size=0.5, inherit.aes = FALSE)
  # } +
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

## save
ggsave("figures/study_area/study_area_wTracks_rawfiltered_alldata.png", width = 5, height=12 )
ggsave("figures/study_area/study_area_wTracks_rawfiltered_alldata_wPoints.png", width = 5, height=12 )
  