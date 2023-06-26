### Bowen Island Lidar ###

library(lidR)
library(rlas)
library(sf)
library(terra)
library(tidyverse)
library(mapview)
library(raster)
library(RCSF)
library(ggpubr)
library(rayshader)
library(osmdata)


## Load catalog and then clip

# Note that "keep first" keeps only the first return.
ctg_pc <- readLAScatalog("bowen/dat/point_cloud", filter = "drop_z_below 0 -drop_class 7 ")
las_check(ctg_pc)
ctg_pc
plot(ctg_pc, mapview = TRUE)

# Refuge boundary
refuge <- st_read("bowen/dat/shp/Park Site - Outline.shp")
refuge
refuge
mapview(refuge)

# Clip Area
aoi_extent <- st_bbox(refuge)
buff <- 50
aoi_extent[1] <- aoi_extent[1] - buff
aoi_extent[2] <- aoi_extent[2] - buff
aoi_extent[3] <- aoi_extent[3] + buff
aoi_extent[4] <- aoi_extent[4] + buff
aoi_extent <- aoi_extent %>% 
  st_as_sfc()
mapview(aoi_extent) + mapview(refuge)


# Clip
pc_aoi <- clip_roi(ctg_pc, aoi_extent)
plot(pc_aoi, bg = "white", size = 5)

# Filter out Class 7 noise - not needed now.
pc_aoi <- filter_poi(pc_aoi, Classification != 7L)
table(pc_aoi$Classification)

hist(pc_aoi$Z)
plot(pc_aoi, bg = "white")

writeLAS(pc_aoi, "bowen/dat/point_cloud/processed/pc_refuge.laz")
pc_aoi <- readLAS("bowen/dat/point_cloud/processed/pc_refuge.laz")

plot(pc_aoi, bg = "white")

# Look at a cross-section

## Functions
plot_crossection <- function(las,
                             p1 = c(min(las@data$X), mean(las@data$Y)),
                             p2 = c(max(las@data$X), mean(las@data$Y)),
                             width = 4, colour_by = NULL) {
  colour_by <- rlang::enquo(colour_by)
  data_clip <- clip_transect(las, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + 
    geom_point(size = 0.5) + 
    coord_equal() + 
    theme_minimal() +
    xlab("UTM E") + 
    ylab("Elevation (m)")
  
  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")
  
  return(p)
}

transect_length <- 400
p1 <- c(469250, 5464928)

pc_cs-plot_crossection(pc_aoi, p1 , c(p1[1]+transect_length, p1[2]), colour_by = factor(Classification))
pc_cs

first <- filter_first(pc_aoi)
plot(first, size = 3, bg = "white", color = "Intensity", legend = TRUE)


# Create DTM --------------------------------------------------------------

# This is what ground looks like. 
gnd <- filter_ground(pc_aoi)
plot(gnd, size = 3, bg = "white")

# This is a test to see if ground classification will improve class values.

# pmf - classifies some canopy areas as grand 
pc_gnd0 <- classify_ground(pc_aoi, algorithm = pmf(ws = 5, th = 3))
filter_gnd <- filter_ground(pc_gnd0)
plot(filter_gnd, bg = "white")
p_gnd0 <- plot_crossection(pc_gnd0, p1 , p2, colour_by = factor(Classification))
p_gnd0

# Alternate pmf - similar to no classification
ws <- seq(3, 12, 3)
th <- seq(0.1, 1.5, length.out = length(ws))
pc_gnd1 <- classify_ground(pc_aoi, algorithm = pmf(ws = ws, th = th))
filter_gnd <- filter_ground(pc_gnd1)
plot(filter_gnd, bg = "white")
p_gnd1 <- plot_crossection(pc_gnd1, p1 , p2, colour_by = factor(Classification))
p_gnd1

# Cloth - missing ground spots
pc_gnd2 <- classify_ground(pc_aoi, algorithm = csf())
filter_gnd <- filter_ground(pc_gnd2)
plot(filter_gnd, bg = "white")
p_gnd2 <- plot_crossection(pc_gnd2, p1 , p2, colour_by = factor(Classification))
p_gnd2

# Compare
ggarrange(pc_cs, p_gnd0, p_gnd1,  p_gnd2, ncol = 1)


# Original point cloud data is best.
# Some errors. Other methods may be useful
dtm <- rasterize_terrain(pc_aoi, 0.5, knnidw())

plot(dtm)
plot(refuge, add = TRUE, colour = "")

writeRaster(dtm, "bowen/dat/dtm.tif", overwrite=TRUE)
dtm <- rast("bowen/dat/dtm.tif")


# Make a hillshade
dtm_prod <- terrain(dtm, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col =gray(0:30/30), legend = FALSE)
plot(refuge, add = TRUE, colour = "")
#dtm_hillshade <- terra::crop(dtm_hillshade, clip_area)

# Rayshader

dtm <- raster(dtm)
dtm_m <- raster_to_matrix(dtm)

dtm_m %>%
  height_shade() %>% 
  plot_3d(heightmap = dtm_m, zscale = 0.5)

dtm_m %>%
  height_shade() %>% 
  add_water(detect_water(dtm_m < 0)) %>% 
  plot_3d(heightmap = dtm_m, zscale = 0.5)

dtm_m %>%
  height_shade() %>% 
  add_overlay(sphere_shade(dtm_m), alphalayer=0.5) %>% 
  add_water(detect_water(dtm_m < 0)) %>% 
  plot_map(heightmap = dtm_m, zscale = 0.5)

dtm_m %>%
  height_shade() %>% 
  add_overlay(sphere_shade(dtm_m), alphalayer=0.5) %>% 
  add_water(detect_water(dtm_m < 0)) %>% 
  plot_3d(solid = TRUE, shadow = TRUE)

dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
  add_shadow(ambient_shade(dtm_m)) %>%
  add_water(detect_water(dtm_m < 1)) %>%
  plot_3d(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 0.5)

# Add Lambertian shading
dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
  add_shadow(lamb_shade(dtm_m, zscale = 0.5),0) %>%
  add_water(detect_water(dtm_m < 1)) %>%
  plot_3d(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 0.5)

dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
  add_shadow(lamb_shade(dtm_m, zscale = 0.5),0) %>%
  add_water(detect_water(dtm_m < 1)) %>%
  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 0.5)

dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
  add_shadow(lamb_shade(dtm_m, zscale = 0.5),0) %>%
  add_shadow(texture_shade(dtm_m,detail=8/10,contrast=9,brightness = 11)) %>%
  add_water(detect_water(dtm_m < 1)) %>%
  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 0.5)

# Add ambient shade
dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
  # add_shadow(lamb_shade(dtm_m, zscale = 0.5),0) %>%
  # add_shadow(texture_shade(dtm_m,detail=8/10,contrast=9,brightness = 11)) %>%
  # add_shadow(ambient_shade(dtm_m), 0) %>%
  # add_water(detect_water(dtm_m < 1)) %>%
  generate_polygon_overlay(refuge, extent = refuge, 
                                    heightmap = dtm_m) %>% 
  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 0.5)





# Add boundary from OSM data

available_features(bb_ll)

bbox <- st_bbox(aoi_extent)

bb_ll = bbox %>%
  st_as_sfc() %>%
  st_transform(crs = 4326) %>%
  st_bbox()
bb_ll

osm_bbox <- osmdata(bb_ll)
osm_bbox


protect = opq(osm_bbox) %>% 
  add_osm_feature("parking") %>% 
  osmdata_sf() 
protect


ref_line <- st_read("bowen/dat/shp/ref_line.shp")
mapview(ref_line)
ref_line <- st_as_sf(ref_line)
# ref_line <- st_cast(st_geometry(refuge), "linestring")

dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
#  add_shadow(lamb_shade(dtm_m, zscale = 0.5),0) %>%
#  add_shadow(texture_shade(dtm_m,detail=8/10,contrast=9,brightness = 11)) %>%
#  add_shadow(ambient_shade(dtm_m), 0) %>%
  add_water(detect_water(dtm_m < 1)) %>%
  add_overlay(generate_line_overlay(ref_line, extent = refuge, 
                                    heightmap = dtm_m), 1) %>% 
  plot_3d(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 0.5)

rgl.snapshot('bowen/out/3dplot3.png', fmt = 'png')




# Vegetation Surface ------------------------------------------------------
# Using dtm
veg_n <- pc_aoi - dtm
plot(veg_n, size = 3, )

# using point cloud and interpolation. This is better. 
veg_n <- normalize_height(pc_aoi, knnidw())
plot(veg_n, size = 3, bg = "white")

chm <- rasterize_canopy(veg_n, res = 1, algorithm = p2r())
plot(chm, bg = "white")
writeRaster(chm, "bowen/chm.tif", overwrite = TRUE)

ttops_all <- locate_trees(veg_n, lmf(ws = 5))
# ttops using normalized creates more ttops! That makes sense.

# Filter out trees less than 10m 
ttops <- filter(ttops_all, Z > 10)

# Now add the absolute heights in for 3d visualization
ttops_elev <- extract(dtm, ttops)
ttops$Z_norm <- ttops$Z
ttops$Z <- ttops$Z_norm + ttops_elev$Z 

plot(pc_aoi, bg = "white") |> add_treetops3d(ttops, z = "Z")  
  
plot(chm)
plot(sf::st_geometry(ttops), add = TRUE, pch = 1)
plot(sf::st_geometry(refuge), add = TRUE, pch = 1)
st_write(ttops, "bowen/ttops.gpkg", append = FALSE)

ttops <- st_as_sf(ttops)
st_write(ttops, dsn = "bowen/ttops.shp", layer = "ttops")




# Terrain -----------------------------------------------------------------








