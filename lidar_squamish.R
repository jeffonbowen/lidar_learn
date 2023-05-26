### Squamish Landfill Expansion LiDAR Data ###
{
library(lidR)
library(rlas)
library(lidRviewer)
library(sf)
library(rgdal)
library(terra)
library(rayshader)
library(tidyverse)
library(mapview)
library(raster)
library(RCSF)
library(ggpubr)
}

# Load footprint

fp <- readOGR(dsn = "squamish/shp/footprint.shp") %>% st_as_sfc()
plot(fp)
mapview(fp)

# Make a rectangular aoi
aoi_extent <- st_bbox(fp)
buff <- 250
aoi_extent[1] <- aoi_extent[1] - buff
aoi_extent[2] <- aoi_extent[2] - buff
aoi_extent[3] <- aoi_extent[3] + buff
aoi_extent[4] <- aoi_extent[4] + buff

aoi <- st_as_sfc(aoi_extent)
mapview(aoi) + mapview(fp)

# Import and Tile ---------------------------------------------------------

# Squamish Landfill Expansion area: 
# 092g075_3_4_3 is the north tile
# 092g075_3_4_1 is the south tile

# Raw data for this area is point cloud, DSM (.laz), and DEM (raster).
# The DSM is not complete for the Squamish landfill area. Will need to create
# DSM from the point cloud.
# Note that point cloud is in BC Albers. Best to work in UTM as soon as possible. 

## Will work with point cloud and DSM only.

# Have a look at the DSM
dsm <- readLAS("squamish/dsm/raw/bc_092g075_3_4_1_xyes_8_utm10_20170601_dsm.laz")
plot(dsm, color = "Classification", bg = "white")

# Read one of the tiles and have a look at the data. 

pc <- readLAS("squamish/point_cloud/raw/bc_092g075_3_4_3_xyes_12_bcalb_2018.laz")
pc <- readLAS("squamish/point_cloud/raw/bc_092g075_3_4_3_xyes_12_bcalb_2018.laz",
              filter = "-keep_first -drop_z_below 0 -drop_class 7 ")
las_check(pc)
names(pc)
table(pc$Classification)
table(pc$ReturnNumber)
plot(pc, color = "Intensity", bg = "white", size = 4)
plot(pc, colour = "Intensity", bg = "white", size = 4, backend = "lidRviewer")

# Load all tiles as a catalog
ctg_pc <- readLAScatalog("squamish/point_cloud/raw", 
                         filter = "-keep_first -drop_z_below 0 -drop_class 7")

# Clip pc now so that we are not transforming the entire catalog. 
aoi_albers <- st_transform(aoi, crs = 3005)
pc <- clip_roi(ctg_pc, aoi_albers)

# Now transform in to UTM Zone 10N
pc <- st_transform(pc, crs = 26910)

# Not a lot of duplicates, but helps
pc <- filter_duplicates(pc)

plot(pc, bg = "white")

# Save
writeLAS(pc, "squamish/point_cloud/processed/pc_aoi.laz")
pc <- readLAS("squamish/point_cloud/processed/pc_aoi.laz")

las_check(pc)
plot(pc, color = "Classification", bg = "white", legend = TRUE)

# Can get rid of transmission line by removing unclassified.
# "vg" is veg / ground.
# However, much of teh unclassified is vegetation and it means losing
# a lof points that are veg. For this project area, best to leave with 
# unclassified in. 
pc_v_g <- filter_poi(pc, Classification != 1L)
plot(pc_v_g, bg = "white")

# Look at a cross-section
transect_length <- 400
p1 <- c(489407, 5515288)
plot_crossection(pc, p1 , c(p1[1]+transect_length, p1[2]), colour_by = factor(Classification))

p2 <- c(489591, 5515263)
plot_crossection(pc, p2 , c(p2[1]+transect_length, p2[2]), colour_by = factor(Classification))

## DEM / DTM
# I compared the dem provided in the lidar catalogue with making one from
# the point cloud. They seem identical. Either way. 

dem <- rast("squamish/dem/bc_092g075_xli1m_utm10_2018.tif")
dem <- terra::crop(dem, aoi)
plot(dem)

# Make our own. Can use either pc or pc_vg. Doesn't matter
dtm <- rasterize_terrain(pc, res = 1, algorithm = tin())
plot(dtm)

x <- plot(pc, bg = "white")
add_dtm3d(x, dtm)

# Make a hillshade
dtm_prod <- terrain(dtm, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col =gray(0:30/30), legend = FALSE)
plot(fp, add = TRUE)

## Canopy
# Using dtm
veg_n <- pc_v_g - dtm
plot(veg_n, size = 3, bg = "white")

# using point cloud and interpolation. This is better. 
veg_n <- normalize_height(pc, knnidw())
y <- plot(veg_n, size = 3, bg = "white")
add_dtm3d(y, dtm)

# CHM
col = height.colors(25)
chm <- rasterize_canopy(veg_n, res = 0.5, 
                        algorithm = p2r(subcircle = 0.15, na.fill = tin()))
plot(chm, legend = TRUE, col = col)
plot(fp, add = TRUE)

# Tree tops 
# ttops using normalized creates more ttops! That makes sense.
ttops <- locate_trees(veg_n, lmf(ws = 5))
ttops <- filter(ttops, Z >= 10)

# Now add the absolute heights in for 3d visualization
ttops_elev <- extract(dtm, ttops)
ttops$Z_norm <- ttops$Z
ttops$Z <- ttops$Z_norm + ttops_elev$Z 

# Not interested in trees outside of footprint
ttops <- st_intersection(ttops, fp)

# Summarize data for footprint
ggplot(ttops, aes(Z_norm)) + 
  geom_histogram(binwidth = 10, color = "white", just = 1) +
  theme_bw() +
  xlab("Tree Height (m)") +
  ylab("Number of Trees") +
  stat_bin(binwidth=10, geom='text', color='black', size=3,
           aes(label=..count..), position=position_dodge(), 
           vjust=-1.1,
           hjust=2.5)

chm_fp <- terra::crop(chm, fp)
plot(chm_fp, col = height.colors(25))
plot(sf::st_geometry(ttops), add = TRUE, pch = 1)

# to Zoom in
e <- drawExtent()
# or
class      : Extent 
xmin       : 489630 
xmax       : 489730 
ymin       : 5515200 
ymax       : 5515331 

zoom(chm_fp, e, col = height.colors(25))
plot(sf::st_geometry(ttops), add = TRUE, pch = 1)

plot(e)

e <- drawExtent() ## interactive!
plot(e, asp = 1, xlab = "", ylab = "", axes = FALSE)


zoom

st_write(ttops, "squamish/out/ttops.gpkg", append = FALSE)

x <- plot(pc, bg = "white", size = 4)
add_treetops3d(x, ttops)


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







