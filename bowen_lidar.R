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


## Load catalog and then clip

# Note that "keep first" keeps only the first return.
ctg_pc <- readLAScatalog("bowen", filter = "-keep_first")
las_check(ctg_pc)
ctg_pc
plot(ctg_pc)

# Refuge boundary
refuge <- st_read("bowen/Park Site - Outline.shp")
refuge
plot

# Clip
pc_aoi <- clip_roi(ctg_pc, 
                         refuge)
plot(pc_aoi, bg = "white", size = 5)

# Filter out Class 7 noise
pc_aoi <- filter_poi(pc_aoi, Classification != 7L)
table(pc_aoi$Classification)

hist(pc_aoi$Z)
plot(pc_aoi, bg = "white")

writeLAS(pc_aoi, "bowen/pc_refuge.laz")
pc_aoi <- readLAS("bowen/pc_refuge.laz")

plot(pc_aoi, bg = "white", color = "Intensity")

# Look at a cross-section
p1 <- c(469400, 5464794)
p2 <- c(469400, 5465245)
# Eastwest
p1 <- c(469250, 5464928)
p2 <- c(469700, 5464928)

pc_cs <- plot_crossection(pc_aoi, p1 , p2, colour_by = factor(Classification))
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
dtm <- rasterize_terrain(pc_aoi, 1, knnidw())
plot(dtm)
writeRaster(dtm, "bowen/dtm.tif", overwrite=TRUE)


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

## Functions
plot_crossection <- function(las,
                             p1 = c(min(las@data$X), mean(las@data$Y)),
                             p2 = c(max(las@data$X), mean(las@data$Y)),
                             width = 4, colour_by = NULL)
{
  colour_by <- rlang::enquo(colour_by)
  data_clip <- clip_transect(las, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 0.5) + 
    coord_equal() + 
    theme_minimal()
  
  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")
  
  return(p)
}








