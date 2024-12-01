library(terra)

elev_1_2_300 <- rast("outputs/elevation_buffers_1-300.tif")
elev_400 <- rast("outputs/elevation_buffers_400.tif")
elev_500 <- rast("outputs/elevation_buffer_500.tif")
elev_600 <- rast("outputs/elevation_buffer_600.tif")
elev_700 <- rast("outputs/elevation_buffer_700.tif")
elev_800 <- rast("outputs/elevation_buffer_800.tif")
elev_900 <- rast("outputs/elevation_buffer_900.tif")
elev_1000 <- rast("outputs/elevation_buffer_1000.tif")

names(elev_1_2_300) <- c("elev_0", "elev_100", "elev_200", "elev_300")

elev_1_2_300$elev_400 <- elev_400$Elev_400
elev_1_2_300$elev_500 <- elev_500$Elev_500
elev_1_2_300$elev_600 <- elev_600$Elev_600
elev_1_2_300$elev_700 <- elev_700$Elev_700
elev_1_2_300$elev_800 <- elev_800$Elev_800
elev_1_2_300$elev_900 <- elev_900$Elev_900
elev_1_2_300$elev_1000 <- elev_1000$Elev_1000

elev_buffers <- elev_1_2_300

writeRaster(elev_buffers, "outputs/elevation_buffers.tif",
            overwrite = T)