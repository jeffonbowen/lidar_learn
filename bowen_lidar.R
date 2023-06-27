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
library(bcmaps)


## Load catalog and then clip

# Note that "keep first" keeps only the first return.
ctg <- readLAScatalog("bowen/las", filter = "drop_z_below 0 -drop_class 7 ")
las_check(ctg)
ctg
plot(ctg, mapview = TRUE)

# Bowen boundary
# Comes in bc albers.Lidar data is in UTM Xone 10.
munis <- bcmaps::municipalities()
bowen <- filter(munis, ADMIN_AREA_ABBREVIATION == "Bowen Island") %>% 
  st_transform(crs = 3157)
mapview(bowen)

# Clip
pc <- clip_roi(ctg, bowen)
plot(pc_aoi, bg = "white", size = 5)

# Save for later
writeLAS(pc, "bowen/las_processed/pc.laz")
pc <- readLAS("bowen/las_processed/pc.laz")

plot(pc, bg = "white")


# Create DTM --------------------------------------------------------------

dtm <- rasterize_terrain(ctg, 0.5, knnidw())

plot(dtm)
plot(refuge, add = TRUE, colour = "")

writeRaster(dtm, "bowen/dat/dtm.tif", overwrite=TRUE)
dtm <- rast("bowen/dat/dtm.tif")



# Load existing dtm -------------------------------------------------------

ls <- list.files("bowen/tif", "tif$", full.names=TRUE) # the $ excludes aux-files
ic <- sprc(lapply(ls, rast))
dtm <- mosaic(ic, fun = "mean") %>% 
  crop(bowen) %>% 
  mask(bowen)  

# writeRaster(dtm, "bowen/tif/bowen_dtm.tif")
dtm <- rast("bowen/tif/bowen_dtm.tif")

# Get rid of the band of ocean
dtm <- ifel(dtm < 0, NA, dtm)

# Make a hillshade
dtm_prod <- terrain(dtm, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col =gray(0:30/30), legend = FALSE)


# Rayshader

dtm <- raster(dtm)
dtm_m <- raster_to_matrix(dtm)
# Resize for prototyping
dtm_m <- resize_matrix(dtm_m, 0.25)

# Simple sphereshade
dtm_m %>%
  height_shade() %>% 
  save_png("bowen/out/bowen_sphereshade.png")

# Add water
# dtm_m %>%
#   height_shade() %>% 
#   add_overlay(sphere_shade(dtm_m), alphalayer=0.5) %>% 
#   add_water(detect_water(dtm_m < 0)) %>% 
# #  plot_map(heightmap = dtm_m, zscale = 0.5) %>% 
#   save_png("bowen/out/bowen_sphereshade_water.png")

# Add rayshade and Lambertian shading
dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
  add_shadow(lamb_shade(dtm_m, zscale = 0.5),0) %>%
#  add_water(detect_water(dtm_m < 1)) %>%
#  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 0.5) %>% 
  save_png("bowen/out/bowen_lamshade.png")

# Add texture shade
# final version
dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
  add_shadow(lamb_shade(dtm_m, zscale = 1),0) %>%
  add_shadow(texture_shade(dtm_m,detail=8/10,contrast=9,brightness = 11)) %>%
  add_water(detect_water(dtm_m, min_area = 100), color = "imhof1") %>%
#  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 1) %>%
  save_png("bowen/out/bowen_terrain_shade.png")

# Add ambient shade
# I have not found this to be any better.
# dtm_m %>%
#   height_shade() %>%
#   add_shadow(ray_shade(dtm_m)) %>%
#   add_shadow(lamb_shade(dtm_m, zscale = 0.5),0) %>%
#   add_shadow(texture_shade(dtm_m,detail=8/10,contrast=9,brightness = 11)) %>%
#   add_shadow(ambient_shade(dtm_m), 0) %>%
#   add_water(detect_water(dtm_m, min_area = 200)) %>%
#   save_png("bowen/out/bowen_ambientshade_water.png")

#3d
dtm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dtm_m)) %>%
  add_shadow(lamb_shade(dtm_m, zscale = 0.5),0) %>%
  add_shadow(texture_shade(dtm_m,detail=8/10,contrast=9,brightness = 11)) %>%
  add_water(detect_water(dtm_m, min_area = 100), color = "imhof1") %>%
  plot_3d(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 10)

render_snapshot()
rgl.snapshot('bowen/out/3dplot1.png', fmt = 'png')


# Vegetation Surface ------------------------------------------------------

pc <- readLAS("bowen/las_processed/pc.laz")

dsm <- rasterize_canopy(pc)

dsm <- raster(dsm)
dsm_m <- raster_to_matrix(dsm)
# Resize for prototyping
dsm_m <- resize_matrix(dsm_m, 0.5)

dsm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dsm_m)) %>%
  add_shadow(lamb_shade(dsm_m, zscale = 1),0) %>%
  add_shadow(texture_shade(dsm_m,detail=8/10,contrast=9,brightness = 11)) %>%
  add_water(detect_water(dsm_m, min_area = 100), color = "imhof1") %>%
  #  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 1) %>%
  save_png("bowen/out/bowen_surface_shade.png")


# CHM
chm <- dsm - dtm
plot(chm)

# using point cloud and interpolation. This is better. 
chm2 <- normalize_height(pc, knnidw())
plot(veg_n, size = 3, bg = "white")

chm <- rasterize_canopy(veg_n, res = 1, algorithm = p2r())
plot(chm, bg = "white")
writeRaster(chm, "bowen/chm.tif", overwrite = TRUE)

