### quick map different datasets along the analysis chain ###

library(mapview)

da_sf = st_as_sf(one, coords = c("Longitude", "Latitude"), 
                 crs = 4326, agr = "constant")
# mapview::mapview(da_sf)

in_sf = st_as_sf(internest, coords = c("Longitude", "Latitude"), 
                 crs = 4326, agr = "constant")
# mapview::mapview(da_sf)

mi_sf = st_as_sf(migration, coords = c("Longitude", "Latitude"), 
                 crs = 4326, agr = "constant")
# mapview::mapview(da_sf)

fg_sf = st_as_sf(foraging, coords = c("Longitude", "Latitude"), 
                 crs = 4326, agr = "constant")
# mapview::mapview(da_sf)

mapview(in_sf) + mapview(mi_sf, col.regions="red") + mapview(fg_sf, col.regions="orange")


da_sf = st_as_sf(tracks_i, coords = c("Longitude", "Latitude"), 
                 crs = 4326, agr = "constant")
mapview(da_sf)

da_sf = st_as_sf(tracks_re, coords = c("x", "y"), 
                 crs = 4326, agr = "constant")

mapview(da_sf)
