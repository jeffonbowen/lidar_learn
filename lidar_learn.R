### Playing Arpuind with LiDAR data ###

library(lidR)
library(rlas)
library(sf)
library(terra)
library(tidyverse)
library(mapview)
library(raster)
library(RCSF)
library(ggpubr)

# Using the Squamish Landfill area as test. 
# 092g075_3_4_3 is the north tile
# 092g075_3_4_1 is the south tile

# Point cloud
pc <- readLAS("dat/bc_092g075_3_4_3_xyes_12_bcalb_2018.laz")
las_check(las)
plot(dsm, colour = "RGB", bg = "white", size = 4, axis = TRUE, legend = TRUE)

ttops <- locate_trees(pc, lmf(ws = 5))
plot(ttops)
st_write(ttops, "dat/ttops.gpkg")


# DEM
dem <- rast("dat/bc_092g075_xli1m_utm10_2018.tif")
plot(dem)


plot(pc, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)



# Load the dsm
dsm <- readLAS("dat/bc_092g075_3_4_1_xyes_8_utm10_20170601_dsm.laz")
plot(dsm)

plot(dsm, color = "RGB", bg = "white", axis = TRUE, legend = TRUE)


dsmr <- rasterize_canopy(dsm, res = 1, algorithm = p2r())
plot(dsmr)

p2 <- plot(dsm)
add_treetops3d(p2, ttops)
add_treetops3d(ttops)

plot(dsm) |> add_treetops3d(ttops)



## Load catalog and then clip

# Note that keep first keeps only the first return.
ctg_pc <- readLAScatalog("dat/point_cloud/", filter = "-keep_first")
las_check(ctg_pc)
ctg_pc
plot(ctg_pc)

# Clip
pc_aoi <- clip_rectangle(ctg_pc, 
                         xleft = 1205240, 
                         ybottom = 533850, 
                         xright = 1206010,
                         ytop = 534840)
plot(pc_aoi, bg = "white")

hist(pc_aoi$Z)
summary(pc_aoi$Z)
table(pc_aoi$Classification)

# Get rid of class 7
pc_aoi <- filter_poi(pc_aoi, Classification != 7L)
table(test$Classification)
plot(pc_aoi, bg = "white")

writeLAS(pc_aoi, "dat/point_cloud/pc_aoi.laz")
pc_aoi <- readLAS("dat/point_cloud/pc_aoi.laz")

first <- filter_first(pc_aoi)
plot(first, size = 3, bg = "white", color = "Classification")

# Look at a cross-section

p1 <- c(1205240, 534375)
p2 <- c(1206006, 534000)
plot_crossection(pc_aoi, p1 , p2, colour_by = factor(Classification))


## Create DTM
# This is what ground looks like
gnd <- filter_ground(pc_aoi)
plot(gnd, size = 3, bg = "white")

dtm_aoi <- rasterize_terrain(pc_aoi, 1, knnidw())
plot(dtm_aoi)
writeRaster(dtm_aoi, "dat/point_cloud/dtm_aoi.tif", overwrite=TRUE)

# This is a test to see if ground classification will improve class values.
# It does not. It makes it worse. 
pc_aoi_gnd <- classify_ground(pc_aoi, algorithm = pmf(ws = 5, th = 3))
table(pc_aoi_gnd$Classification)
plot(pc_aoi_gnd, bg = "white")
plot_crossection(pc_aoi_gnd, p1 , p2, colour_by = factor(Classification))


## Normalize vegetation height
# Using dtm
veg_n <- pc_aoi - dtm_aoi
plot(veg_n1, size = 3, bg = "white")

# using point cloud and interpolation. This is better. 
veg_n <- normalize_height(pc_aoi, knnidw())
plot(veg_n1, size = 3, bg = "white")

chm <- rasterize_canopy(veg_n, res = 1, algorithm = p2r())
plot(chm)
writeRaster(chm, "dat/chm.tif", overwrite = TRUE)

chm2 <- rasterize_canopy(veg_n, res = 0.5, p2r(0.2, na.fill = tin()))
writeRaster(chm2, "dat/chm2.tif", overwrite = TRUE)

ttops <- locate_trees(veg_n, lmf(ws = 5))
ttops <- filter(ttops, Z >= 10)
plot(chm2, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
st_write(ttops, "dat/ttops.gpkg", append = FALSE)

x <- plot(pc_aoi, bg = "white", size = 4)
add_treetops3d(x, ttops)


## Functions
plot_crossection <- function(las,
                             p1 = c(min(las@data$X), mean(las@data$Y)),
                             p2 = c(max(las@data$X), mean(las@data$Y)),
                             width = 4, colour_by = NULL)
{
  colour_by <- rlang::enquo(colour_by)
  data_clip <- clip_transect(las, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 0.5) + coord_equal() + theme_minimal()
  
  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")
  
  return(p)
}







