### Playing Arpuind with LiDAR data ###

library(lidR)
library(rlas)


# 092G075 is teh main site. 

las <- readLAS("dat/bc_092g075_3_4_1_xyes_12_bcalb_2018.laz")

# or//
las <- readLAS("dat/bc_092g075_3_4_1_xyes_12_bcalb_2018.laz", select = "xyz")


las_check(las)

p1 <- plot(las, bg = "white", size = 3)


ttops <- locate_trees(las, lmf(ws = 5))

add_treetops3d(p1, ttops)

plot(las, color = "ScanAngleRank", bg = "white", axis = TRUE, legend = TRUE)




# Load the dsm
dsm <- readLAS("dat/bc_092g075_3_4_1_xyes_8_utm10_20170601_dsm.laz")
p2 <- plot(dsm)
add_treetops3d(p2, ttops)



# Test scripts from package

LASfile <- system.file("extdata", "Topography.laz", package="lidR")
las <- readLAS(LASfile)

dtm <- rasterize_terrain(las, algorithm = tin())
ttops <- locate_trees(las, lmf(ws = 5))

plot_dtm3d(dtm)

x <- plot(las)
add_dtm3d(x, dtm)
add_treetops3d(x, ttops)

plot(las) |> add_dtm3d(dtm) |> add_treetops3d(ttops)





