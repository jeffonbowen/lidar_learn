### Bowen Island Lidar ###

{
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
library(bcdata)
library(future)
}

## Load catalog and then clip

# Note that "keep first" keeps only the first return.
ctg <- readLAScatalog("bowen/las", 
                      select = "xyzc",
                      filter = "drop_z_below 0 -drop_class 7")
las_check(ctg)
ctg
plot(ctg, mapview = TRUE)

# Bowen boundary
# Comes in bc albers. Lidar data is in UTM Zone 10.
munis <- bcmaps::municipalities()
bowen <- filter(munis, ADMIN_AREA_ABBREVIATION == "Bowen Island") %>% 
  st_transform(crs = 3157)
mapview(bowen)

# Clip a region of interest and use parallel processing
# Not run - very long and no need yet. 
# plan(multisession)
# opt_chunk_size(ctg) <- 400
# opt_chunk_buffer(ctg) <- 40
# # opt_output_files(ctg) <- paste0(tempdir(), "/{*}_classified")
# ctg2 <- clip_roi(ctg, bowen)
# plot(pc_aoi, bg = "white", size = 5)

# Save for later
# writeLAS(pc, "bowen/las_processed/pc.laz")
# pc <- readLAS("bowen/las_processed/pc.laz")


# Create DTM From Lidar ---------------------------------------------------

dtm_lidar <- rasterize_terrain(ctg, 0.5, knnidw())

plot(dtm_lidar)
plot(refuge, add = TRUE, colour = "")

writeRaster(dtm_lidar, "bowen/dat/dtm_lidar.tif", overwrite=TRUE)


# Load existing dtm -------------------------------------------------------

# There is already a dtm available at 1 m. We can just use that rather than 
# create a new one. 
ls <- list.files("bowen/tif", "tif$", full.names=TRUE) # the $ excludes aux-files
ic <- terra::sprc(lapply(ls, terra::rast))
dtm <- mosaic(ic, fun = "mean") %>% 
  crop(bowen) %>% 
  mask(bowen)  

# Get rid of the band of ocean
dtm <- ifel(dtm < 0, NA, dtm)

writeRaster(dtm, "bowen/tif/bowen_dtm.tif", overwrite = TRUE)
dtm <- rast("bowen/tif/bowen_dtm.tif")
plot(dtm)

# Make a hillshade
dtm_prod <- terrain(dtm, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col =gray(0:30/30), legend = FALSE)


# Rayshader ---------------------------------------------------------------

dtm <- raster(dtm)
dtm_m <- raster_to_matrix(dtm)
# Resize for prototyping
dtm_m <- resize_matrix(dtm_m, scale = 0.25)

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
  add_water(detect_water(dtm_m < 0.1)) %>%
  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 0.5) 
#  save_png("bowen/out/bowen_lamshade.png")

# Add texture shade
# final version
dtm_m %>%
  height_shade(texture = topo.colors(256)) %>%
  add_shadow(ray_shade(dtm_m)) %>%
  add_shadow(lamb_shade(dtm_m, zscale = 1),0) %>%
  add_shadow(texture_shade(dtm_m,detail=10/10,contrast=9,brightness = 11)) %>%
  add_water(detect_water(dtm_m, min_area = 250), color = "imhof1") %>%
  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 1) 
#  save_png("bowen/out/bowen_terrain_shade7.png")

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
  add_water(detect_water(dtm_m, min_area = 200), color = "imhof1") %>%
  plot_3d(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 10)

render_snapshot()
rgl.snapshot('bowen/out/3dplot1.png', fmt = 'png')


# Vegetation Surface ------------------------------------------------------

dsm <- rasterize_canopy(ctg, res = 1, p2r())
dsm <- mask(dsm, bowen)
# Get rid of the band of ocean
dsm <- ifel(dsm < 0, NA, dsm)

writeRaster(dsm, "bowen/tif/bowen_dsm.tif", overwrite = TRUE)
dsm <- rast("bowen/tif/bowen_dsm.tif")
plot(dsm)

dsm <- raster(dsm)
dsm_m <- raster_to_matrix(dsm)
# Resize for prototyping
dsm_m <- resize_matrix(dsm_m, 0.1)

dsm_m %>%
  height_shade() %>%
  add_shadow(ray_shade(dsm_m)) %>%
  add_shadow(lamb_shade(dsm_m, zscale = 1),0) %>%
  add_shadow(texture_shade(dsm_m,detail=8/10,contrast=9,brightness = 11)) %>%
  add_water(detect_water(dsm_m, min_area = 100), color = "imhof1") %>%
  #  plot_map(heightmap = dtm_m, solid = TRUE, shadow = TRUE, zscale = 1) %>%
  save_png("bowen/out/bowen_surface_shade.png")


# CHM --------------------------------------------------------------------

# Test 1 tile
las <- readLAS("bowen/las/bc_092g034_3_1_4_xyes_8_utm10_2019.laz", 
               select = "xyzc",
               filter = "drop_z_below 0 -drop_class 7")
las_check(las)
las
plot(las)

dtm <- rast("bowen/tif/bowen_dtm.tif")

# using hybrid method
nlas <- normalize_height(las, dtm)
chm_las <- rasterize_canopy(nlas, res = 1, algorithm = p2r())
chm_las_0 <- ifel(chm_las < 0, 0, chm_las)
col <- height.colors(25)
plot(chm_las_0, col = col)

# Seems that simple is ok. 
# Simple method for whole island.
chm <- dsm-dtm
writeRaster(chm, "bowen/tif/chm.tif", overwrite = TRUE)

# Load
chm <- rast("bowen/tif/chm.tif")
col <- height.colors(50)
plot(chm, col = col)

# Hybrid method. Has never completed.
# Load the dtm if not already in memory
dtm <- rast("bowen/tif/bowen_dtm.tif")

# To use normalize_height on ctg, must specify temp directory.
opt_output_files(ctg) <-  paste0(tempdir(), "/{*}_norm")
ctg_norm <- normalize_height(ctg, dtm)

chm <- rasterize_canopy(ctg_norm, res = 1, algorithm = p2r())
writeRaster(chm, "bowen/tif/chm.tif", overwrite = TRUE)


# Plot
col <- height.colors(25)
plot(chm, col = col)

# Some clean-up options. Not much difference. 
fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)
filled <- terra::focal(chm, w, fun = fill.na)
smoothed <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
plot(filled, col = col)
plot(smoothed, col = col)

# Rayshader on canopy?
chm <- raster(chm)
chm_m <- raster_to_matrix(chm)
chm_m <- resize_matrix(chm_m, scale = 0.5)

chm_m %>%
  height_shade(texture = height.colors(256)) %>%
  add_shadow(ray_shade(chm_m)) %>%
  add_shadow(lamb_shade(chm_m, zscale = 1),0) %>%
  add_shadow(texture_shade(chm_m,detail=10/10,contrast=9,brightness = 11)) %>%
#  plot_map(heightmap = chm_m, solid = TRUE, shadow = TRUE, zscale = 1) 
  save_png("bowen/out/bowen_veg_height4.png")



