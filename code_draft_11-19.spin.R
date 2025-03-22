## Olivia Guswiler
## Final Project Script
### 2024-12-17

library(tidyverse)
library(terra)
library(amt)
library(raster)
library(MuMIn)
library(glmmTMB)
library(doParallel)
library(ggplot2)
library(paletteer)
library(cowplot)

## Tidy GPS Data ####

# reference data with study site info
ref <- read_csv("data_raw/prugh_movebank_ref.csv") %>% 
  dplyr::select("animal-id", "study-site") %>% 
  rename(ID = "animal-id", study_site =`study-site`)

# all meso telemetry data
all_data <- readRDS("data_raw/Prugh_provided_files/WPPP_mesopredator_GPS_locations.rds") %>%
  rename(date_time = Acquisition.Time,
         Lat = GPS.Latitude, Long = GPS.Longitude, ID = AnimalID) %>% 
  left_join(ref, "ID")

# filter by study site, remove disperal data, and tidy
data_ok <- all_data %>%
  filter(study_site == "Okanogan") %>% 
  group_by(ID) %>%
  arrange(date_time) %>%
  ungroup() %>%
  mutate(Sp = ifelse(Species == "BOB", "Bobcat", "Coyote")) %>%
  # excluding MVBOB71M after dispersal date
  filter(ID != "MVBOB71M" | date_time < as.POSIXct("2019-09-24 00:00:00"))

# save as new RDS
write_rds(data_ok, "outputs/Okanogan_GPS_locations.rds")

## Tidy Covariate Data ####

## Hiking Trails 
trail <- rast("rasters/trail_raster.tif") %>% 
  # reproject to appropriate crs
  project("EPSG:32610", method = "bilinear")

# make data binary
trail[!is.na(trails)] <- 1  # change any values to 1

# save raster
writeRaster(trail, "rasters/trail.tif",
            overwrite = T)

# give data meaningful name
names(trail) <- "trail_dist"

# Calculate Euclidean distance to trails
trail_dist <- distance(trails,
                            filename = "rasters/trail_dist.tif",
                            filetype = "GTiff", gdal = c("COMPRESS=DEFLATE", "BIGTIFF=YES"),
                            overwrite = T)


## Human Footprint Index (2019)
hfi <- rast("data_raw/Prugh_provided_files/human_footprint_2019.tif") %>% 
  project("EPSG:32610", method = "bilinear") %>% 
  # crop to trails extent
  crop(trails)

names(hfi) <- "hfi_0"

# save as new file
writeRaster(hfi, "rasters/hfi.tif",
            overwrite = T)


## Digital Elevation Model (2024)
# read in DEM sections
D1 <- rast("rasters/USGS_13_n48w120_20240617.tif")
D2 <- rast("rasters/USGS_13_n48w121_20240617.tif")
D3 <- rast("rasters/USGS_13_n49w120_20240617.tif")
D4 <- rast("rasters/USGS_13_n49w121_20240617.tif")

# merge sections into single raster layer
D12 <- merge(D1, D2)
D34 <- merge(D3, D4)
DEM <- merge(D12, D34) %>% 
  # reproject to appropriate crs
  project("EPSG:32610", method = "bilinear") %>% 
  crop(trails)

# rename layer for clarity
names(DEM) <- "elev_0"

writeRaster(DEM, "rasters/elevation.tif",
            overwrite = T)


## Percent Tree Canopy Cover (2019)
canopy <- rast("rasters/NLCD_CanopyCover.tiff") %>%
  project("EPSG:32610", method = "bilinear") %>% 
  crop(trails)
# change values over 87 to NA, these are an error from reprojection
canopy[canopy>87] <- NA

names(canopy) <- "canopy_0"

writeRaster(canopy, "rasters/canopy.tif",
            overwrite = T)

## Calculate Available Step Lengths and Turn Angles ####
rm(list=ls())

# load previously created file
meso <- readRDS("outputs/Okanogan_GPS_locations.rds")

# select variables, nest by individual animal, and make 'amt' track
meso_track <- meso %>%
  data.frame() %>%
  dplyr::select(ID, Species, Sex, date_time, E, N) %>%
  nest(.by = ID) %>%
  mutate(trk = lapply(data, function(d){
    make_track(d, .x = E, .y = N, .t = date_time, crs = "EPSG:32610")}))

# summarise sampling rate by individual
track_summary <- meso_track %>% 
  mutate(sr = lapply(trk, summarize_sampling_rate, time_unit = "hour")) %>%
  dplyr::select(ID, sr) %>% unnest(cols = c(sr)) %>% 
  left_join(distinct(data.frame(meso)[,c("ID", "Sex", "Species")])) %>%
  arrange(Species, median)
print(track_summary)


# separate by species, unnest
bobcat <- meso_track %>% filter(grepl("BOB", ID))
coyote <- meso_track %>% filter(grepl("COY", ID))


# calculate available points from sl and ta empirical distributions

# BOBCAT
bobcat1 <- bobcat[bobcat$ID != "MVBOB87M",] %>% # exclude due to few bursts
  # drop individuals with fewer than 50 points  
  filter(map_int(trk, nrow) > 50) %>% 
  mutate(stp = map(trk, ~ .x %>%
                     track_resample(rate = hours(4), tolerance = minutes(10)) %>%
                     filter_min_n_burst(min_n = 3) %>%
                     # calculate step lengths
                     steps_by_burst() %>% 
                     # draw 10 random points based on gamma and vonmises dists
                     random_steps(sl_distr = fit_distr(.$sl_, "gamma"),
                                  ta_distr = fit_distr(.$ta_, "vonmises"),
                                  n_control = 10) %>% 
                     # calculate log of step length and cos of turn angle
                     mutate(cos_ta_ = cos(ta_),
                            log_sl_ = log(sl_)))) %>% 
  dplyr::select(-data, -trk) %>%
  unnest(cols = stp) %>%
  left_join(distinct(data.frame(meso)[,c("ID", "Sex", "Species")])) %>%
  mutate(case_binary = ifelse(case_ == TRUE, 1, 0))



# COYOTE
coyote1 <- coyote %>%  
  filter(map_int(trk, nrow) > 50) %>% 
  mutate(stp = map(trk, ~ .x %>%
                     track_resample(rate = hours(4), tolerance = minutes(10)) %>%
                     filter_min_n_burst(min_n = 3) %>%
                     steps_by_burst() %>% 
                     random_steps(sl_distr = fit_distr(.$sl_, "gamma"),
                                  ta_distr = fit_distr(.$ta_, "vonmises"),
                                  n_control = 10) %>% 
                     mutate(cos_ta_ = cos(ta_),
                            log_sl_ = log(sl_)))) %>% 
  dplyr::select(-data, -trk) %>%
  unnest(cols = stp) %>%
  left_join(distinct(data.frame(meso)[,c("ID", "Sex", "Species")])) %>%
  mutate(case_binary = ifelse(case_ == TRUE, 1, 0))

# save data as RDS
write_rds(bobcat1, "outputs/bobcat_steps.rds")
write_rds(coyote1, "outputs/coyote_steps.rds")


coy1 <- coyote_extract %>% 
  filter(ID == "MVCOY68F")

sl_gamma <- fit_distr(coy1$sl_, "gamma")

hist(coy1$sl_, freq = FALSE,
     main = "Observed Steps",
     xlab = "Step Length (m)")

curve(dgamma(x, 
             shape = sl_gamma$params$shape,
             scale = sl_gamma$params$scale),
      add = T,
      lwd = 2,
      col = "blue")

hist(rgamma(651320, 
            shape = sl_gamma$params$shape,
            scale = sl_gamma$params$scale),
     freq = FALSE,
     main = "Random Steps",
     xlab = "Step Length (m)")


ta_VM <- fit_distr(coy1$ta_, "vonmises")
ta_VM

hist(coy1$ta_, freq = FALSE,
     main = "Observed Turn Angles",
     xlab = "Turn Angle (radians)")

TA_dist <- data.frame(x = seq(-3.15, 3.15, length.out = 651320))
TA_dist$y <- circular::dvonmises(
  x = TA_dist$x, 
  mu = ta_VM$params$mu,
  kappa = ta_VM$params$kappa)

lines(y~x, TA_dist, add=T, lwd = 2, col = "red")

# 500X550


## Create Raster Buffers ####
rm(list=ls())

trail <- rast("rasters/distance_to_trail.tif")
elev <- rast("rasters/elevation.tif")
canopy <- rast("rasters/canopy.tif")
hfi <- rast("rasters/hfi.tif")


# set range of buffers
buffer_sizes = seq(100, 1000, by = 100)


  ## Elevation
for (i in 1:length(buffer_sizes)) {
  buff_i <- buffer_sizes[i]
  cat("Starting buffer size =",buff_i,"\n")
  site_buffers_i <- focalMat(x = elev,           
                             d = buff_i,          
                             type="circle")
  elev_buffer_i <- focal(elev, w = site_buffers_i, 
                         fun = "mean")
  elev$TMP <- elev_buffer_i
  names(elev)[which(names(elev)=="TMP")] <- paste0("Elev_",buff_i)
}

# save
writeRaster(elev, "rasters/elevation_buffers.tif",
            overwrite = T)


  ## Human Footprint Index
for (i in 1:length(buffer_sizes)) {
  buff_i <- buffer_sizes[i]
  cat("Starting buffer size =",buff_i,"\n")
  site_buffers_i <- focalMat(x = hfi,           
                             d = buff_i,          
                             type="circle")
  hfi_buffer_i <- focal(hfi, w = site_buffers_i, 
                         fun = "mean")
  hfi$TMP <- hfi_buffer_i
  names(hfi)[which(names(hfi)=="TMP")] <- paste0("hfi_",buff_i)
}

writeRaster(hfi, "rasters/hfi_buffers.tif",
            overwrite = T)

  ## Distance to Trails
for (i in 1:length(buffer_sizes)) {
  buff_i <- buffer_sizes[i]
  cat("Starting buffer size =",buff_i,"\n")
  site_buffers_i <- focalMat(x = trail,           
                             d = buff_i,          
                             type="circle")
  trail_buffer_i <- focal(trail, w = site_buffers_i, 
                        fun = "mean")
  trail$TMP <- trail_buffer_i
  names(trail)[which(names(trail)=="TMP")] <- paste0("trail_",buff_i)
}

writeRaster(trail, "rasters/trail_buffers.tif",
            overwrite = T)


  ## Percent Canopy Cover
for (i in 1:length(buffer_sizes)) {
  buff_i <- buffer_sizes[i]
  cat("Starting buffer size =",buff_i,"\n")
  site_buffers_i <- focalMat(x = canopy,           
                             d = buff_i,          
                             type="circle")
  canopy_buffer_i <- focal(canopy, w = site_buffers_i, 
                          fun = "mean")
  canopy$TMP <- canopy_buffer_i
  names(canopy)[which(names(canopy)=="TMP")] <- paste0("canopy_",buff_i)
}

writeRaster(canopy, "rasters/canopy_buffers.tif",
            overwrite = T)



## Extract Covariates ####
rm(list=ls())

# load previously created files
trail <- rast("rasters/trail_dist.tif")
elev <- rast("rasters/elevation_buffers.tif")
canopy <- rast("rasters/canopy_buffers.tif")
hfi <- rast("rasters/hfi_buffers.tif")
bobcat_extract <- readRDS("outputs/bobcat_steps.rds")
coyote_extract <- readRDS("outputs/coyote_steps.rds")

# define season and set weights
bobcat_extract <- bobcat_extract %>% 
  mutate(season = ifelse(month(t2_) %in% 4:11, "summer", "winter"),
         w = ifelse(bobcat_extract$case_binary==1,1,10000))
coyote_extract <- coyote_extract %>% 
  mutate(season = ifelse(month(t2_) %in% 4:11, "summer", "winter"),
         w = ifelse(coyote_extract$case_binary==1,1,10000))


# extract covariates
bobcat_canopy <- raster::extract(canopy, bobcat_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  cbind(bobcat_extract)
coyote_canopy <- raster::extract(canopy, coyote_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  cbind(coyote_extract)

bobcat_elev <- raster::extract(elev, bobcat_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  cbind(bobcat_extract)
coyote_elev <- raster::extract(elev, coyote_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  cbind(coyote_extract)

bobcat_hfi <- raster::extract(hfi, bobcat_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  cbind(bobcat_extract)
coyote_hfi <- raster::extract(hfi, coyote_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  cbind(coyote_extract)


  ### Determine scales of effect ----
    #### Bobcat Non-Winter ----
## Human Footprint Index
bobcat_hfi_sum <- bobcat_hfi %>% # filter data by season  
  filter(season == "summer")

# create data frame to hold AIC values for comparison
bob_hfi_scales_sum <- data.frame(covariate = "hfi",
                           scale = seq(0,1000,by = 100),
                           AIC = NA,
                           delta_AIC = NA,
                           w = NA)
# list to hold models
hfi_models <- list()

# calculate AICs and enter in df
for(i in 1:nrow(bob_hfi_scales_sum)){
  cov_i <- paste0(bob_hfi_scales_sum$covariate[i],"_",bob_hfi_scales_sum$scale[i])
  data_i <- bobcat_hfi_sum[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_hfi_sum$w)
  hfi_models[[i]] <- model_i
  bob_hfi_scales_sum$AIC[i] <- AIC(model_i)
}

# calculate delta AIC values
for (i in 1:length(bob_hfi_scales_sum)) {
  min_i <- min(bob_hfi_scales_sum$AIC)
  bob_hfi_scales_sum$delta_AIC <- bob_hfi_scales_sum$AIC - min_i
}

# calculate weights and add to data frame
hfi_weights <- as.vector(Weights(AIC(hfi_models[[1]],hfi_models[[2]],
                                         hfi_models[[3]],hfi_models[[4]],
                                         hfi_models[[5]],hfi_models[[6]],
                                         hfi_models[[7]],hfi_models[[8]],
                                         hfi_models[[9]],hfi_models[[10]],
                                         hfi_models[[11]])))
bob_hfi_scales_sum <- bob_hfi_scales_sum %>% 
  mutate(w = hfi_weights)

# lowest AIC value
bob_hfi_scales_sum[which(bob_hfi_scales_sum$AIC==min(bob_hfi_scales_sum$AIC)),]
#   covariate scale    AIC delta_AIC         w
# 4       hfi   300 200435         0 0.1001909



## % Canopy Cover
bobcat_canopy_sum <- bobcat_canopy %>%
  filter(season == "summer") %>%         ###fix this
  mutate(canopy_0 = Canopy_Cover) %>% 
  dplyr::select(-Canopy_Cover)

bob_canopy_scales_sum <- data.frame(covariate = "canopy",
                                 scale = seq(0,1000,by = 100),
                                 AIC = NA,
                                 delta_AIC = NA,
                                 w = NA)

canopy_models <- list()

for(i in 1:nrow(bob_canopy_scales_sum)){
  cov_i <- paste0(bob_canopy_scales_sum$covariate[i],"_",bob_canopy_scales_sum$scale[i])
  data_i <- bobcat_canopy_sum[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_canopy_sum$w)
  canopy_models[[i]] <- model_i
  bob_canopy_scales_sum$AIC[i] <- AIC(model_i)
}

for (i in 1:length(bob_canopy_scales_sum)) {
  min_i <- min(bob_canopy_scales_sum$AIC)
  bob_canopy_scales_sum$delta_AIC <- bob_canopy_scales_sum$AIC - min_i
}

canopy_weights <- as.vector(Weights(AIC(canopy_models[[1]],canopy_models[[2]],
                                        canopy_models[[3]],canopy_models[[4]],
                                        canopy_models[[5]],canopy_models[[6]],
                                        canopy_models[[7]],canopy_models[[8]],
                                        canopy_models[[9]],canopy_models[[10]],
                                        canopy_models[[11]])))
bob_canopy_scales_sum <- bob_canopy_scales_sum %>% 
  mutate(w = canopy_weights)


bob_canopy_scales_sum[which(bob_canopy_scales_sum$AIC==min(bob_canopy_scales_sum$AIC)),]
#   covariate scale      AIC delta_AIC         w
# 1    canopy     0 200245.5         0 0.9999961



## Elevation
bobcat_elev_sum <- bobcat_elev %>%
  filter(season == "summer")

bob_elev_scales_sum <- data.frame(covariate = "elev",
                                  scale = seq(0,1000,by = 100),
                                  AIC = NA,
                                  delta_AIC = NA,
                                  w = NA)

elev_models <- list()

for(i in 1:nrow(bob_elev_scales_sum)){
  cov_i <- paste0(bob_elev_scales_sum$covariate[i],"_",bob_elev_scales_sum$scale[i])
  data_i <- bobcat_elev_sum[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_elev_sum$w)
  elev_models[[i]] <- model_i
  bob_elev_scales_sum$AIC[i] <- AIC(model_i)
}

for (i in 1:length(bob_elev_scales_sum)) {
  min_i <- min(bob_elev_scales_sum$AIC)
  bob_elev_scales_sum$delta_AIC <- bob_elev_scales_sum$AIC - min_i
}

elev_weights <- as.vector(Weights(AIC(elev_models[[1]],elev_models[[2]],
                                      elev_models[[3]],elev_models[[4]],
                                      elev_models[[5]],elev_models[[6]],
                                      elev_models[[7]],elev_models[[8]],
                                      elev_models[[9]],elev_models[[10]],
                                      elev_models[[11]])))
bob_elev_scales_sum <- bob_elev_scales_sum %>% 
  mutate(w = elev_weights)

bob_elev_scales_sum[which(bob_elev_scales_sum$AIC==min(bob_elev_scales_sum$AIC)),]
#    covariate scale      AIC delta_AIC        w
# 11      elev  1000 200436.1         0 0.091893


    #### Bobcat Winter ----
## Human Footprint Index
bobcat_hfi_win <- bobcat_hfi %>%  
  filter(season == "winter")

bob_hfi_scales_win <- data.frame(covariate = "hfi",
                                 scale = seq(0,1000,by = 100),
                                 AIC = NA,
                                 delta_AIC = NA,
                                 w = NA)
hfi_models <- list()

for(i in 1:nrow(bob_hfi_scales_win)){
  cov_i <- paste0(bob_hfi_scales_win$covariate[i],"_",bob_hfi_scales_win$scale[i])
  data_i <- bobcat_hfi_win[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_hfi_win$w)
  hfi_models[[i]] <- model_i
  bob_hfi_scales_win$AIC[i] <- AIC(model_i)
}

for (i in 1:length(bob_hfi_scales_win)) {
  min_i <- min(bob_hfi_scales_win$AIC)
  bob_hfi_scales_win$delta_AIC <- bob_hfi_scales_win$AIC - min_i
}

hfi_weights <- as.vector(Weights(AIC(hfi_models[[1]],hfi_models[[2]],
                                     hfi_models[[3]],hfi_models[[4]],
                                     hfi_models[[5]],hfi_models[[6]],
                                     hfi_models[[7]],hfi_models[[8]],
                                     hfi_models[[9]],hfi_models[[10]],
                                     hfi_models[[11]])))
bob_hfi_scales_win <- bob_hfi_scales_win %>% 
  mutate(w = hfi_weights)

bob_hfi_scales_win[which(bob_hfi_scales_win$AIC==min(bob_hfi_scales_win$AIC)),]
#   covariate scale      AIC delta_AIC        w
# 7       hfi   600 74551.64         0 0.131403


## % Canopy Cover
bobcat_canopy_win <- bobcat_canopy %>%
  filter(season == "winter") %>% 
  mutate(canopy_0 = Canopy_Cover) %>% 
  dplyr::select(-Canopy_Cover)

bob_canopy_scales_win <- data.frame(covariate = "canopy",
                                    scale = seq(0,1000,by = 100),
                                    AIC = NA,
                                    delta_AIC = NA,
                                    w = NA)
canopy_models <- list()

for(i in 1:nrow(bob_canopy_scales_win)){
  cov_i <- paste0(bob_canopy_scales_win$covariate[i],"_",bob_canopy_scales_win$scale[i])
  data_i <- bobcat_canopy_win[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_canopy_win$w)
  canopy_models[[i]] <- model_i
  bob_canopy_scales_win$AIC[i] <- AIC(model_i)
}

for (i in 1:length(bob_canopy_scales_win)) {
  min_i <- min(bob_canopy_scales_win$AIC)
  bob_canopy_scales_win$delta_AIC <- bob_canopy_scales_win$AIC - min_i
}

canopy_weights <- as.vector(Weights(AIC(canopy_models[[1]],canopy_models[[2]],
                                        canopy_models[[3]],canopy_models[[4]],
                                        canopy_models[[5]],canopy_models[[6]],
                                        canopy_models[[7]],canopy_models[[8]],
                                        canopy_models[[9]],canopy_models[[10]],
                                        canopy_models[[11]])))
bob_canopy_scales_win <- bob_canopy_scales_win %>% 
  mutate(w = canopy_weights)

bob_canopy_scales_win[which(bob_canopy_scales_win$AIC==min(bob_canopy_scales_win$AIC)),]
#   covariate scale      AIC delta_AIC         w
# 1    canopy     0 74550.29         0 0.4679327


## Elevation
bobcat_elev_win <- bobcat_elev %>%
  filter(season == "winter")

bob_elev_scales_win <- data.frame(covariate = "elev",
                                  scale = seq(0,1000,by = 100),
                                  AIC = NA,
                                  delta_AIC = NA,
                                  w = NA)
elev_models <- list()

for(i in 1:nrow(bob_elev_scales_win)){
  cov_i <- paste0(bob_elev_scales_win$covariate[i],"_",bob_elev_scales_win$scale[i])
  data_i <- bobcat_elev_win[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_elev_win$w)
  elev_models[[i]] <- model_i
  bob_elev_scales_win$AIC[i] <- AIC(model_i)
}

for (i in 1:length(bob_elev_scales_win)) {
  min_i <- min(bob_elev_scales_win$AIC)
  bob_elev_scales_win$delta_AIC <- bob_elev_scales_win$AIC - min_i
}

elev_weights <- as.vector(Weights(AIC(elev_models[[1]],elev_models[[2]],
                                      elev_models[[3]],elev_models[[4]],
                                      elev_models[[5]],elev_models[[6]],
                                      elev_models[[7]],elev_models[[8]],
                                      elev_models[[9]],elev_models[[10]],
                                      elev_models[[11]])))
bob_elev_scales_win <- bob_elev_scales_win %>% 
  mutate(w = elev_weights)

bob_elev_scales_win[which(bob_elev_scales_win$AIC==min(bob_elev_scales_win$AIC)),]
#   covariate scale      AIC delta_AIC         w
# 2      elev   100 74537.02         0 0.1931823



      ##### bind scales of effect for all covariates and z-score standardize ----
hfi_b <- dplyr::select(bobcat_hfi, hfi_300, hfi_600, hfi_800) %>% 
  mutate(zhfi_300 = scale(hfi_300),
         zhfi_600 = scale(hfi_600),
         zhfi_800 = scale(hfi_800))
canopy_b <- dplyr::select(bobcat_canopy, canopy_0) %>%
  mutate(zcanopy_0 = scale(canopy_0))
elev_b <- dplyr::select(.data = bobcat_elev, elev_1000, elev_100, elev_200) %>%
  mutate(zelev_100 = scale(elev_1000),
         zelev_1000 = scale(elev_100),
         zelev_200 = scale(elev_200))
trail_b <- raster::extract(trail, bobcat_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  mutate(ztrail_dist = scale(trail_dist))

bobcat_extract1 <- cbind(hfi_b, canopy_b, elev_b, trail_b, 
                         bobcat_extract)


      ##### plot scales of effect ----
plot_bob_hfi_sum <- ggplot(bob_hfi_scales_sum, aes(scale, delta_AIC)) +
  geom_line(color = "#de2c2c", lwd = 0.6) +
  geom_point(fill = "#de2c2c", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "Human Footprint Index (summer)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

plot_bob_hfi_win <- ggplot(bob_hfi_scales_win, aes(scale, delta_AIC)) +
  geom_line(color = "#771818", lwd = 0.6) +
  geom_point(fill = "#771818", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "Human Footprint Index (winter)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

plot_bob_canopy_sum <- ggplot(bob_canopy_scales_sum, aes(scale, delta_AIC)) +
  geom_line(color = "#53c153", lwd = 0.6) +
  geom_point(fill = "#53c153", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "% Canopy Cover (summer)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

plot_bob_canopy_win <- ggplot(bob_canopy_scales_win, aes(scale, delta_AIC)) +
  geom_line(color = "#334e33", lwd = 0.6) +
  geom_point(fill = "#334e33", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "% Canopy Cover (winter)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

plot_bob_elev_sum <- ggplot(bob_elev_scales_sum, aes(scale, delta_AIC)) +
  geom_line(color = "cornflowerblue", lwd = 0.6) +
  geom_point(fill = "cornflowerblue", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "Elevation (summer)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

plot_bob_elev_win <- ggplot(bob_elev_scales_win, aes(scale, delta_AIC)) +
  geom_line(color = "darkblue", lwd = 0.6) +
  geom_point(fill = "darkblue", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "Elevation (winter)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

plot_grid(plot_bob_hfi_sum, plot_bob_canopy_sum, plot_bob_elev_sum,
          plot_bob_hfi_win, plot_bob_canopy_win, plot_bob_elev_win,
          ncol = 3, nrow = 2) # saved 520x780

    #### Coyote ----
    ## Human Footprint Index
coy_hfi_scales <- data.frame(covariate = "hfi",
                             scale = seq(0,1000,by = 100),
                             AIC = NA,
                             delta_AIC = NA,
                             w = NA)
hfi_models <- list()

for(i in 1:nrow(coy_hfi_scales)){
  cov_i <- paste0(coy_hfi_scales$covariate[i],"_",coy_hfi_scales$scale[i])
  data_i <- coyote_hfi[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = coyote_hfi$w)
  hfi_models[[i]] <- model_i
  coy_hfi_scales$AIC[i] <- AIC(model_i)
}

for (i in 1:length(coy_hfi_scales)) {
  min_i <- min(coy_hfi_scales$AIC)
  coy_hfi_scales$delta_AIC <- coy_hfi_scales$AIC - min_i
}

hfi_weights <- as.vector(Weights(AIC(hfi_models[[1]],hfi_models[[2]],
                                     hfi_models[[3]],hfi_models[[4]],
                                     hfi_models[[5]],hfi_models[[6]],
                                     hfi_models[[7]],hfi_models[[8]],
                                     hfi_models[[9]],hfi_models[[10]],
                                     hfi_models[[11]])))
coy_hfi_scales <- coy_hfi_scales %>% 
  mutate(w = hfi_weights)

coy_hfi_scales[which(coy_hfi_scales$AIC==min(coy_hfi_scales$AIC)),]
#   covariate scale      AIC delta_AIC         w
# 9       hfi   800 814309.5         0 0.1955595


    ## % Canopy Cover
coyote_canopy <- coyote_canopy %>%
  mutate(canopy_0 = Canopy_Cover) %>% 
  dplyr::select(-Canopy_Cover)

coy_canopy_scales <- data.frame(covariate = "canopy",
                                scale = seq(0,1000,by = 100),
                                AIC = NA,
                                delta_AIC = NA,
                                w = NA)
canopy_models <- list()

for(i in 1:nrow(coy_canopy_scales)){
  cov_i <- paste0(coy_canopy_scales$covariate[i],"_",coy_canopy_scales$scale[i])
  data_i <- coyote_canopy[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = coyote_canopy$w)
  canopy_models[[i]] <- model_i
  coy_canopy_scales$AIC[i] <- AIC(model_i)
}

for (i in 1:length(coy_canopy_scales)) {
  min_i <- min(coy_canopy_scales$AIC)
  coy_canopy_scales$delta_AIC <- coy_canopy_scales$AIC - min_i
}

canopy_weights <- as.vector(Weights(AIC(canopy_models[[1]],canopy_models[[2]],
                                        canopy_models[[3]],canopy_models[[4]],
                                        canopy_models[[5]],canopy_models[[6]],
                                        canopy_models[[7]],canopy_models[[8]],
                                        canopy_models[[9]],canopy_models[[10]],
                                        canopy_models[[11]])))
coy_canopy_scales <- coy_canopy_scales %>% 
  mutate(w = canopy_weights)

coy_canopy_scales[which(coy_canopy_scales$AIC==min(coy_canopy_scales$AIC)),]
#   covariate scale      AIC delta_AIC w
# 1    canopy     0 813760.6         0 1


    ## Elevation
coy_elev_scales <- data.frame(covariate = "elev",
                              scale = seq(0,1000,by = 100),
                              AIC = NA,
                              delta_AIC = NA,
                              w = NA)
elev_models <- list()

for(i in 1:nrow(coy_elev_scales)){
  cov_i <- paste0(coy_elev_scales$covariate[i],"_",coy_elev_scales$scale[i])
  data_i <- coyote_elev[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = coyote_elev$w)
  elev_models[[i]] <- model_i
  coy_elev_scales$AIC[i] <- AIC(model_i)
}

for (i in 1:length(coy_elev_scales)) {
  min_i <- min(coy_elev_scales$AIC)
  coy_elev_scales$delta_AIC <- coy_elev_scales$AIC - min_i
}

elev_weights <- as.vector(Weights(AIC(elev_models[[1]],elev_models[[2]],
                                      elev_models[[3]],elev_models[[4]],
                                      elev_models[[5]],elev_models[[6]],
                                      elev_models[[7]],elev_models[[8]],
                                      elev_models[[9]],elev_models[[10]],
                                      elev_models[[11]])))
coy_elev_scales <- coy_elev_scales %>% 
  mutate(w = elev_weights)

coy_elev_scales[which(coy_elev_scales$AIC==min(coy_elev_scales$AIC)),]
#   covariate scale      AIC delta_AIC        w
# 7      elev   600 814337.3         0 0.178463



      ##### bind scales of effect for all covariates and z-score standardize ----
hfi_c <- dplyr::select(coyote_hfi, hfi_800, hfi_0) %>% 
  mutate(zhfi_800 = scale(hfi_800),
         zhfi_0 = scale(hfi_0))
canopy_c <- dplyr::select(coyote_canopy, canopy_0) %>%
  mutate(zcanopy_0 = scale(canopy_0))
elev_c <- dplyr::select(.data = coyote_elev, elev_700, elev_400, elev_600) %>%
  mutate(zelev_700 = scale(elev_700),
         zelev_400 = scale(elev_400),
         zelev_600 = scale(elev_600))
trail_c <- raster::extract(trail, coyote_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  mutate(ztrail_dist = scale(trail_dist))

coyote_extract1 <- cbind(hfi_c, canopy_c, elev_c, trail_c, 
                         coyote_extract)


      ##### plot scales of effect ----
plot_coy_hfi_soe <- ggplot(coy_hfi_scales, aes(scale, delta_AIC)) +
  geom_line(color = "firebrick", lwd = 0.6) +
  geom_point(fill = "firebrick", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "Human Footprint Index") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

plot_coy_canopy_soe <- ggplot(coy_canopy_scales, aes(scale, delta_AIC)) +
  geom_line(color = "darkgreen", lwd = 0.6) +
  geom_point(fill = "darkgreen", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "% Canopy Cover") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

plot_coy_elev_soe <- ggplot(coy_elev_scales, aes(scale, delta_AIC)) +
  geom_line(color = "blue", lwd = 0.6) +
  geom_point(fill = "blue", cex = 3, shape = 21) +
  scale_x_continuous(breaks = seq(0,1000, by = 200)) +
  labs(x = "", y = "", title = "Elevation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()


plot_grid(plot_coy_hfi_soe, plot_coy_canopy_soe,
          plot_coy_elev_soe, ncol = 3, nrow = 1) # saved 520x780


  ### Final touches before SSF ----
bobcat_extract2 <- bobcat_extract1 %>% 
  # create unique step ID for each animal
  mutate(step_id_ = paste(ID, step_id_, sep = "_")) %>%
  group_by(ID) %>% 
  # calculate sample size for each animal (dividing by 11 because of available points)
  mutate(n = n()/11) %>%
  ungroup()

coyote_extract2 <- coyote_extract1 %>% 
  # create unique step ID for each animal
  mutate(step_id_ = paste(ID, step_id_, sep = "_")) %>%
  group_by(ID) %>% 
  mutate(n = n()/11) %>%
  ungroup()

# save fully processed data for SSF
saveRDS(bobcat_extract2, "outputs/bobcat_covariates.rds")
saveRDS(coyote_extract2, "outputs/coyote_covariates.rds")




## Fit Models ####
rm(list=ls())

  ## all data
bob_dat_all <- readRDS("outputs/bobcat_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>%
  # select individuals with more than 100 fixes
  filter(n >= 100)
coy_dat_all <- readRDS("outputs/coyote_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>%
  filter(n >= 100)

  ## seasonal data
bob_dat_sum <- readRDS("outputs/bobcat_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>% 
  filter(season == "summer") %>% 
  filter(n >= 100)
bob_dat_win <- readRDS("outputs/bobcat_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>% 
  filter(season == "winter") %>% 
  filter(n >= 100)

coy_dat_sum <- readRDS("outputs/coyote_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>% 
  filter(season == "summer") %>% 
  filter(n >= 100)
coy_dat_win <- readRDS("outputs/coyote_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>% 
  filter(season == "winter") %>% 
  filter(n >= 100)


# set number of cores
nt <- min(detectCores()) - 2

# fit seasonal mixed-effects SSF with random slopes and intercepts
bob_ssf_sum <- glmmTMB(case_ ~ -1 + zhfi_300 + (0 + zhfi_300 | ID) +
                         I(zhfi_300^2) + (0 + I(zhfi_300^2)) +
                         ztrail_dist + (0 + ztrail_dist | ID) +
                         ztrail_dist:zhfi_300 + (0 + ztrail_dist:zhfi_300 | ID) +
                         zcanopy_0 + (0 + zcanopy_0 | ID) +
                         I(zcanopy_0^2) + (0 + I(zcanopy_0^2) | ID) +
                         zelev_1000 + (0 + zelev_1000 | ID) +
                         log_sl_ + (0 + log_sl_ | ID) + cos_ta_ + (0 + cos_ta_ | ID) +
                         # stratum
                         (1 | step_id_), 
                       family=poisson,
                       weights = w,
                       na.action = na.exclude,
                       doFit=T,
                       data = bob_dat_sum,
                       # Tell glmmTMB not to change the last standard deviation, all other values are estimated freely
                       map = list(theta = factor(c(NA,1:8))),
                       # Set the value of the standard deviation of the strata (the last ranef) to large constant value
                       start = list(theta = c(rep(0, times = 8),log(1000))),
                       control = glmmTMBControl(parallel = nt))
bob_ssf_win <- glmmTMB(case_ ~ -1 + zhfi_600 + (0 + zhfi_600 | ID) +
                         I(zhfi_600^2) + (0 + I(zhfi_600^2)) +
                         ztrail_dist + (0 + ztrail_dist | ID) +
                         ztrail_dist:zhfi_600 + (0 + ztrail_dist:zhfi_600 | ID) +
                         zcanopy_0 + (0 + zcanopy_0 | ID) +
                         I(zcanopy_0^2) + (0 + I(zcanopy_0^2) | ID) +
                         zelev_100 + (0 + zelev_100 | ID) +
                         log_sl_ + (0 + log_sl_ | ID) + cos_ta_ + (0 + cos_ta_ | ID) +
                         (1 | step_id_), 
                       family=poisson,
                       weights = w,
                       na.action = na.exclude,
                       doFit=T,
                       data = bob_dat_win,
                       map = list(theta = factor(c(NA,1:8))),
                       start = list(theta = c(rep(0, times = 8),log(1000))),
                       control = glmmTMBControl(parallel = nt))

coy_ssf_sum <- glmmTMB(case_ ~ -1 + zhfi_800 + (0 + zhfi_800 | ID) +
                         I(zhfi_800^2) + (0 + I(zhfi_800^2)) +
                         ztrail_dist + (0 + ztrail_dist | ID) +
                         ztrail_dist:zhfi_800 + (0 + ztrail_dist:zhfi_800 | ID) +
                         zcanopy_0 + (0 + zcanopy_0 | ID) +
                         I(zcanopy_0^2) + (0 + I(zcanopy_0^2) | ID) +
                         zelev_700 + (0 + zelev_700 | ID) +
                         log_sl_ + (0 + log_sl_ | ID) + cos_ta_ + (0 + cos_ta_ | ID) +
                         (1 | step_id_), 
                       family=poisson,
                       weights = w,
                       na.action = na.exclude,
                       doFit=T,
                       data = coy_dat_sum,
                       map = list(theta = factor(c(NA,1:8))),
                       start = list(theta = c(rep(0, times = 8),log(1000))),
                       control = glmmTMBControl(parallel = nt))
coy_ssf_win <- glmmTMB(case_ ~ -1 + zhfi_0 + (0 + zhfi_0 | ID) +
                         I(zhfi_0^2) + (0 + I(zhfi_0^2)) +
                         ztrail_dist + (0 + ztrail_dist | ID) +
                         ztrail_dist:zhfi_0 + (0 + ztrail_dist:zhfi_0 | ID) +
                         zcanopy_0 + (0 + zcanopy_0 | ID) +
                         I(zcanopy_0^2) + (0 + I(zcanopy_0^2) | ID) +
                         zelev_400 + (0 + zelev_400 | ID) +
                         log_sl_ + (0 + log_sl_ | ID) + cos_ta_ + (0 + cos_ta_ | ID) +
                         (1 | step_id_), 
                       family=poisson,
                       weights = w,
                       na.action = na.exclude,
                       doFit=T,
                       data = coy_dat_win,
                       map = list(theta = factor(c(NA,1:8))),
                       start = list(theta = c(rep(0, times = 8),log(1000))),
                       control = glmmTMBControl(parallel = nt))


# fit full mixed-effects SSF with random slopes and intercepts
coy_ssf <- glmmTMB(case_ ~ -1 + zhfi_800 + (0 + zhfi_800 | ID) +
                     I(zhfi_800^2) + (0 + I(zhfi_800^2)) +
                     ztrail_dist + (0 + ztrail_dist | ID) +
                     I(ztrail_dist^2) + (0 + I(ztrail_dist^2) | ID) +
                     ztrail_dist:zhfi_800 + (0 + ztrail_dist:zhfi_800 | ID) +
                     zcanopy_0 + (0 + zcanopy_0 | ID) +
                     I(zcanopy_0^2) + (0 + I(zcanopy_0^2) | ID) +
                     zelev_600 + (0 + zelev_600 | ID) +
                     log_sl_ + (0 + log_sl_ | ID) + cos_ta_ + (0 + cos_ta_ | ID) +
                     (1 | step_id_), 
                   family=poisson,
                   weights = w,
                   na.action = na.exclude,
                   doFit=T,
                   data = coy_dat_all,
                   map = list(theta = factor(c(NA,1:9))),
                   start = list(theta = c(rep(0, times = 9),log(1000))),
                   control = glmmTMBControl(parallel = nt))


## Model Summaries ----
bob_ssf_sum <- readRDS("outputs/bobcat_ssf_summer.rds")
bob_ssf_win <- readRDS("outputs/bobcat_ssf_winter.rds")
coy_ssf_sum <- readRDS("outputs/coyote_ssf_summer.rds")
coy_ssf_win <- readRDS("outputs/coyote_ssf_winter.rds")
coy_ssf <- readRDS("outputs/coyote_ssf.rds")
  ### seasonal ----
summary(bob_ssf_sum)
    # Random effects:
    #   
    #   Conditional model:
    # Groups   Name                 Variance  Std.Dev.
    # ID       zhfi_300              1.000000 1.00000 
    # ID.1     ztrail_dist           4.964413 2.22810 
    # ID.2     ztrail_dist:zhfi_300  4.741590 2.17752 
    # ID.3     zcanopy_0             0.048832 0.22098 
    # ID.4     I(zcanopy_0^2)        0.002358 0.04856 
    # ID.5     zelev_1000            0.499570 0.70680 
    # ID.6     log_sl_               0.015345 0.12387 
    # ID.7     cos_ta_               0.001597 0.03996 
    # step_id_ (Intercept)          73.872790 8.59493 
    # Number of obs: 87989, groups:  ID, 16; step_id_, 7999
    # 
    # Conditional model:
    #                      Estimate Std. Error z value Pr(>|z|)    
    # zhfi_300              0.63263    0.26477   2.389 0.016876 *  
    # I(zhfi_300^2)        -0.12498    0.01572  -7.949 1.88e-15 ***
    # ztrail_dist           1.07332    0.56511   1.899 0.057522 .  
    # zcanopy_0             0.23435    0.06077   3.856 0.000115 ***
    # I(zcanopy_0^2)       -0.24870    0.02159 -11.519  < 2e-16 ***
    # zelev_1000           -0.11072    0.18601  -0.595 0.551702    
    # log_sl_              -0.34926    0.03407 -10.253  < 2e-16 ***
    # cos_ta_               0.01884    0.02139   0.880 0.378638    
    # zhfi_300:ztrail_dist  0.85977    0.55468   1.550 0.121138      
summary(bob_ssf_win)
    # Random effects:
    #   
    #   Conditional model:
    # Groups   Name                 Variance Std.Dev.
    # ID       zhfi_600             1.00000  1.0000  
    # ID.1     ztrail_dist          9.23071  3.0382  
    # ID.2     I(ztrail_dist^2)     8.57725  2.9287  
    # ID.3     ztrail_dist:zhfi_600 1.33361  1.1548  
    # ID.4     zcanopy_0            0.07807  0.2794  
    # ID.5     I(zcanopy_0^2)       0.03582  0.1893  
    # ID.6     zelev_100            2.14622  1.4650  
    # ID.7     log_sl_              0.11164  0.3341  
    # ID.8     cos_ta_              0.07875  0.2806  
    # step_id_ (Intercept)          5.27988  2.2978  
    # Number of obs: 32725, groups:  ID, 16; step_id_, 2975
    # 
    # Conditional model:
    #                      Estimate Std. Error z value Pr(>|z|)    
    # zhfi_600              0.98782    0.27076   3.648 0.000264 ***
    # I(zhfi_600^2)        -0.32578    0.02235 -14.577  < 2e-16 ***
    # ztrail_dist           1.04977    0.78020   1.346 0.178460    
    # I(ztrail_dist^2)     -3.06467    0.76562  -4.003 6.26e-05 ***
    # zcanopy_0            -0.12592    0.08057  -1.563 0.118076    
    # I(zcanopy_0^2)       -0.42856    0.05852  -7.323 2.42e-13 ***
    # zelev_100             1.90359    0.37759   5.041 4.62e-07 ***
    # log_sl_              -1.10676    0.08501 -13.019  < 2e-16 ***
    # cos_ta_              -0.41815    0.08014  -5.218 1.81e-07 ***
    # zhfi_600:ztrail_dist -0.11016    0.31676  -0.348 0.728011    

summary(coy_ssf_sum)
    # Random effects:
    #   
    #   Conditional model:
    # Groups   Name                 Variance  Std.Dev.
    # ID       zhfi_800               1.00000  1.0000 
    # ID.1     ztrail_dist          311.62596 17.6529 
    # ID.2     I(ztrail_dist^2)     420.62140 20.5091 
    # ID.3     ztrail_dist:zhfi_800  39.86232  6.3137 
    # ID.4     zcanopy_0              0.16896  0.4110 
    # ID.5     I(zcanopy_0^2)         0.28826  0.5369 
    # ID.6     zelev_700              2.27327  1.5077 
    # ID.7     log_sl_                0.06351  0.2520 
    # ID.8     cos_ta_                0.03888  0.1972 
    # step_id_ (Intercept)            2.45030  1.5653 
    # Number of obs: 238216, groups:  ID, 16; step_id_, 21656
    # 
    # Conditional model:
    #                      Estimate Std. Error z value Pr(>|z|)    
    # zhfi_800              1.14314    0.25710    4.45 8.74e-06 ***
    # I(zhfi_800^2)        -0.74245    0.01320  -56.26  < 2e-16 ***
    # ztrail_dist          20.31961    4.41608    4.60 4.20e-06 ***
    # I(ztrail_dist^2)      9.47714    5.13037    1.85  0.06471 .  
    # zcanopy_0            -0.33524    0.10416   -3.22  0.00129 ** 
    # I(zcanopy_0^2)       -0.59045    0.13523   -4.37 1.26e-05 ***
    # zelev_700            -0.59088    0.37800   -1.56  0.11801    
    # log_sl_              -0.29001    0.06311   -4.60 4.32e-06 ***
    # cos_ta_              -0.25844    0.05086   -5.08 3.75e-07 ***
    # zhfi_800:ztrail_dist  1.84949    1.58141    1.17  0.24220    
summary(coy_ssf_win)
    # Random effects:
    #   
    #   Conditional model:
    # Groups   Name               Variance  Std.Dev.
    # ID       zhfi_0               1.00000  1.0000 
    # ID.1     ztrail_dist        241.56768 15.5424 
    # ID.2     I(ztrail_dist^2)   203.80484 14.2760 
    # ID.3     ztrail_dist:zhfi_0   2.70783  1.6455 
    # ID.4     zcanopy_0            0.10084  0.3175 
    # ID.5     I(zcanopy_0^2)       0.26048  0.5104 
    # ID.6     zelev_400            2.94372  1.7157 
    # ID.7     log_sl_              0.04949  0.2225 
    # ID.8     cos_ta_              0.03047  0.1746 
    # step_id_ (Intercept)          3.95764  1.9894 
    # Number of obs: 119218, groups:  ID, 15; step_id_, 10838
    # 
    # Conditional model:
    #                    Estimate Std. Error z value Pr(>|z|)    
    # zhfi_0              0.77019    0.26551    2.90  0.00372 ** 
    # I(zhfi_0^2)        -0.64800    0.01693  -38.26  < 2e-16 ***
    # ztrail_dist        16.77794    4.01898    4.17 2.98e-05 ***
    # I(ztrail_dist^2)    9.50861    3.69366    2.57  0.01004 *  
    # zcanopy_0          -0.14949    0.08453   -1.77  0.07700 .  
    # I(zcanopy_0^2)     -0.60879    0.13350   -4.56 5.11e-06 ***
    # zelev_400          -0.05026    0.44437   -0.11  0.90995    
    # log_sl_            -0.44247    0.05767   -7.67 1.68e-14 ***
    # cos_ta_            -0.25133    0.04864   -5.17 2.38e-07 ***
    # zhfi_0:ztrail_dist -0.14882    0.43512   -0.34  0.73233    

exp(confint(bob_ssf_sum))
exp(summary(bob_ssf_sum)$coefficients$cond[,2])
    #                           2.5 %       97.5 %   Estimate Std. Error
    # zhfi_300              1.1204146    3.1631376  1.8825582  1.303127             
    # I(zhfi_300^2)         0.8557314    0.9101310  0.8825122  1.015847             
    # ztrail_dist           0.9663101    8.8543697  2.9250755  1.759637             
    # zcanopy_0             1.1221424    1.4239951  1.2640907  1.062658
    # I(zcanopy_0^2)        0.7475043    0.8135228  0.7798152  1.021825             
    # zelev_1000            0.6217041    1.2889900  0.8951930  1.204436             
    # log_sl_               0.6596603    0.7538982  0.7052068  1.034652             
    # cos_ta_               0.9771688    1.0626508  1.0190138  1.021624 
    # zhfi_300:ztrail_dist  0.7966072    7.0071464  2.3626137  1.741392
exp(confint(bob_ssf_win))
exp(summary(bob_ssf_win)$coefficients$cond[,2])
    #                           2.5 %      97.5 %    Estimate Std. Error
    # zhfi_600             1.57955626   4.5653481  2.68537224   1.310957             
    # I(zhfi_600^2)        0.69102722   0.7542932  0.72196754   1.022599             
    # ztrail_dist          0.61915502  13.1832555  2.85700522   2.181919             
    # I(ztrail_dist^2)     0.01040711   0.2092815  0.04666922   2.150336 
    # zcanopy_0            0.75288911   1.0325087  0.88168279   1.083905             
    # I(zcanopy_0^2)       0.58085083   0.7306173  0.65144431   1.060267             
    # zelev_100            3.20122581  14.0643488  6.70992969   1.458759             
    # log_sl_              0.27988676   0.3905724  0.33062975   1.088727 
    # cos_ta_              0.56258608   0.7702175  0.65826563   1.083435             
    # zhfi_600:ztrail_dist 0.48142792   1.6664209  0.89569053   1.372673 

exp(confint(coy_ssf_sum))
exp(summary(coy_ssf_sum)$coefficients$cond[,2])
    #                           2.5 %     97.5 %    Estimate Std. Error
    # zhfi_600             1.57955626  4.5653481  2.68537224   1.293172             
    # I(zhfi_600^2)        0.69102722  0.7542932  0.72196754   1.013284            
    # ztrail_dist          0.61915502 13.1832555  2.85700522   82.770948           
    # I(ztrail_dist^2)     0.01040711  0.2092815  0.04666922   169.080328 
    # zcanopy_0            0.75288911  1.0325087  0.88168279   1.109781             
    # I(zcanopy_0^2)       0.58085083  0.7306173  0.65144431   1.144803             
    # zelev_100            3.20122581 14.0643488  6.70992969   1.459362             
    # log_sl_              0.27988676  0.3905724  0.33062975   1.065142 
    # cos_ta_              0.56258608  0.7702175  0.65826563   1.052176             
    # zhfi_600:ztrail_dist 0.48142792  1.6664209  0.89569053   4.861830 
exp(confint(coy_ssf_win))
exp(summary(coy_ssf_win)$coefficients$cond[,2])
    #                           2.5 %       97.5 %     Estimate Std. Error
    # zhfi_0             1.283770e+00 3.634855e+00 2.160166e+00   1.304092           
    # I(zhfi_0^2)        5.060117e-01 5.407425e-01 5.230890e-01   1.017079          
    # ztrail_dist        7.338381e+03 5.099520e+10 1.934482e+07   55.644486
    # I(ztrail_dist^2)   9.671340e+00 1.877538e+07 1.347528e+04   40.191771           
    # zcanopy_0          7.296617e-01 1.016327e+00 8.611474e-01   1.088210 
    # I(zcanopy_0^2)     4.187678e-01 7.067069e-01 5.440093e-01   1.142818           
    # zelev_400          3.980402e-01 2.272051e+00 9.509825e-01   1.559502           
    # log_sl_            5.737881e-01 7.193208e-01 6.424467e-01   1.059361           
    # cos_ta_            7.070351e-01 8.555647e-01 7.777624e-01   1.049847           
    # zhfi_0:ztrail_dist 3.672770e-01 2.021810e+00 8.617218e-01   1.545146 



  ### year-round ----
summary(coy_ssf)
# Random effects:
# Conditional model:
# Groups   Name                 Variance  Std.Dev.
# ID       zhfi_800               1.00000  1.0000 
# ID.1     ztrail_dist          297.61577 17.2515 
# ID.2     I(ztrail_dist^2)     344.78079 18.5683 
# ID.3     ztrail_dist:zhfi_800  31.04218  5.5716 
# ID.4     zcanopy_0              0.12012  0.3466 
# ID.5     I(zcanopy_0^2)         0.26605  0.5158 
# ID.6     zelev_600              2.20698  1.4856 
# ID.7     log_sl_                0.05812  0.2411 
# ID.8     cos_ta_                0.03032  0.1741 
# step_id_ (Intercept)            2.90002  1.7029 
# Number of obs: 357434, groups:  ID, 16; step_id_, 32494
# 
# Conditional model:
#                      Estimate Std. Error z value Pr(>|z|)    
# zhfi_800              1.15559    0.25475    4.54 5.73e-06 ***
# I(zhfi_800^2)        -0.84857    0.01103  -76.92  < 2e-16 ***
# ztrail_dist          19.79579    4.31486    4.59 4.48e-06 ***
# I(ztrail_dist^2)     10.04452    4.64416    2.16  0.03055 *  
# zcanopy_0            -0.26101    0.08781   -2.97  0.00295 ** 
# I(zcanopy_0^2)       -0.53969    0.12970   -4.16 3.17e-05 ***
# zelev_600            -0.37556    0.37213   -1.01  0.31287    
# log_sl_              -0.33032    0.06035   -5.47 4.41e-08 ***
# cos_ta_              -0.24822    0.04476   -5.55 2.93e-08 ***
# zhfi_800:ztrail_dist  1.35320    1.39499    0.97  0.33202    

exp(confint(coy_ssf))
exp(summary(coy_ssf)$coefficients$cond[,2])
#                              2.5 %       97.5 %     Estimate
# zhfi_800              1.927605e+00 5.232552e+00 3.175893e+00  1.290145             
# I(zhfi_800^2)         4.188690e-01 4.373813e-01 4.280251e-01  1.011094            
# ztrail_dist           8.402105e+04 1.862149e+12 3.955499e+08  74.802948           
# I(ztrail_dist^2)      2.565442e+00 2.067282e+08 2.302931e+04  103.976314             
# zcanopy_0             6.484879e-01 9.149285e-01 7.702727e-01  1.091779 
# I(zcanopy_0^2)        4.520777e-01 7.516536e-01 5.829287e-01  1.138489             
# zelev_600             3.312386e-01 1.424475e+00 6.869068e-01  
# log_sl_               6.385270e-01 8.089319e-01 7.186966e-01  
# cos_ta_               7.146643e-01 8.517269e-01 7.801915e-01  
# zhfi_800:ztrail_dist  2.513448e-01 5.958088e+01 3.869799e+00  



# save files
saveRDS(bob_ssf, "outputs/bobcat_ssf.rds")
saveRDS(bob_ssf_sum, "outputs/bobcat_ssf_summer.rds")
saveRDS(bob_ssf_win, "outputs/bobcat_ssf_winter.rds")
saveRDS(coy_ssf, "outputs/coyote_ssf.rds")
saveRDS(coy_ssf_sum, "outputs/coyote_ssf_summer.rds")
saveRDS(coy_ssf_win, "outputs/coyote_ssf_winter.rds")


## Plots ----

# read in data if needed
coy_dat_all <- readRDS("outputs/coyote_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>%
  filter(n >= 100)
bob_dat_sum <- readRDS("outputs/bobcat_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>% 
  filter(season == "summer") %>% 
  filter(n >= 100)
bob_dat_win <- readRDS("outputs/bobcat_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   
         cos_ta_ = cos(ta_)) %>% 
  filter(season == "winter") %>% 
  filter(n >= 100)

coy_ssf <- readRDS("outputs/coyote_ssf.rds")
bob_ssf_sum <- readRDS("outputs/bobcat_ssf_summer.rds")
bob_ssf_win <- readRDS("outputs/bobcat_ssf_winter.rds")

  ### Hex plots ----
hex_plot_fun <- function(dat, zhfi, ylab, title_, legend) {
  ggplot(dat , aes(get(zhfi), ztrail_dist)) +
    geom_bin_2d(aes(fill = after_stat(log(density)))) +
    scale_fill_continuous(type = "viridis", limits=c(-14,-1)) +
    theme_bw() +
    xlim(NA,9) + 
    labs(x = "Human Footprint", y = ylab, title = title_)+
    theme(plot.title = element_text(hjust = 0.5), legend.position = legend) +
    stat_density_2d(colour = "red", breaks = c(0.05), n = 15, size = 1)+
    geom_vline(xintercept = 0, colour = "gray70", size = 1) +
    geom_hline(yintercept = 0, colour = "gray70", size = 1)
}

hex_bob_all <- hex_plot_fun(bob_dat_all, zhfi = "zhfi_800",
                            ylab = "Trail from Distance",
                            title_ = "Human Footprint and Trail Distance",
                            legend = "")
hex_coy_all <- hex_plot_fun(coy_dat_win, zhfi = "zhfi_800",
                            ylab = "Trail from Distance",
                            title_ = "Human Footprint and Trail Distance",
                            legend = "")

hex_bob_all
hex_coy_all

ggsave("plots/hex_bob_all.png", plot = hex_bob_all)
ggsave("plots/hex_coy_all.png", plot = hex_coy_all)


  ### RSS ----
    #### Bobcat Summer ----
b_hfi_pred_sum <- readRDS("outputs/b_hfi_pred_sum.rds")
b_trail_pred_sum <- readRDS("outputs/b_trail_pred_sum.rds")
b_can_pred_sum <- readRDS("outputs/b_can_pred_sum.rds")
b_elev_pred_sum <- readRDS("outputs/b_elev_pred_sum.rds")
## HFI
b_hfi_pred_sum <- data.frame(zhfi_300 = seq(min(bob_dat_sum$zhfi_300),
                                            max(bob_dat_sum$zhfi_300),
                                            length.out=40),
                             ztrail_dist = mean(bob_dat_sum$ztrail_dist),
                             zcanopy_0 = mean(bob_dat_sum$zcanopy_0),
                             zelev_1000 = mean(bob_dat_sum$zelev_1000),
                             log_sl_ = mean(bob_dat_sum$log_sl_),
                             cos_ta_ = mean(bob_dat_sum$cos_ta_),
                             ID = NA,
                             step_id_ = NA,
                             w = NA)

y <- predict(bob_ssf_sum, newdata = b_hfi_pred_sum,
             type = "link", se = TRUE)

b_hfi_pred_sum$y <- exp(as.vector(y$fit))
b_hfi_pred_sum$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
b_hfi_pred_sum$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

b_rss_hfi_sum <- ggplot(b_hfi_pred_sum, aes(zhfi_300)) +
  geom_line(aes(y=y), color="#de2c2c", lwd = 1.2) +
  geom_line(aes(y=ciu), color="#de2c2c", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#de2c2c", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Human Footprint", y = "RSS", title = "Bobcat Summer") +
  theme_bw()

## trail
b_trail_pred_sum <- data.frame(zhfi_300 = mean(bob_dat_sum$zhfi_300),
                               ztrail_dist = seq(min(bob_dat_sum$ztrail_dist),
                                                 max(bob_dat_sum$ztrail_dist),
                                                 length.out=40),
                               zcanopy_0 = mean(bob_dat_sum$zcanopy_0),
                               zelev_1000 = mean(bob_dat_sum$zelev_1000),
                               log_sl_ = mean(bob_dat_sum$log_sl_),
                               cos_ta_ = mean(bob_dat_sum$cos_ta_),
                               ID = NA,
                               step_id_ = NA,
                               w = NA)

y <- predict(bob_ssf_sum, newdata = b_trail_pred_sum,
             type = "link", se = TRUE)

b_trail_pred_sum$y <- exp(as.vector(y$fit))
b_trail_pred_sum$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
b_trail_pred_sum$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

b_rss_trail_sum <- ggplot(b_trail_pred_sum, aes(ztrail_dist)) +
  geom_line(aes(y=y), color="#99795b", lwd = 1.2) +
  geom_line(aes(y=ciu), color="#99795b", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#99795b", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,26) +
  labs(x = "Trail Proximity", y = "", title = "Bobcat Summer") +
  theme_bw()


## canopy
b_can_pred_sum <- data.frame(zhfi_300 = mean(bob_dat_sum$zhfi_300),
                             ztrail_dist = mean(bob_dat_sum$ztrail_dist),
                             zcanopy_0 = seq(min(bob_dat_sum$zcanopy_0),
                                             max(bob_dat_sum$zcanopy_0),
                                             length.out=40),
                             zelev_1000 = mean(bob_dat_sum$zelev_1000),
                             log_sl_ = mean(bob_dat_sum$log_sl_),
                             cos_ta_ = mean(bob_dat_sum$cos_ta_),
                             ID = NA,
                             step_id_ = NA,
                             w = NA)

y <- predict(bob_ssf_sum, newdata = b_can_pred_sum,
             type = "link", se = TRUE)

b_can_pred_sum$y <- exp(as.vector(y$fit))
b_can_pred_sum$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
b_can_pred_sum$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

b_rss_can_sum <- ggplot(b_can_pred_sum, aes(zcanopy_0)) +
  geom_line(aes(y=y), color="#53c153", lwd = 1.2) +
  geom_line(aes(y=ciu), color="#53c153", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#53c153", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "% Canopy Cover", y = "", title = "Bobcat Summer") +
  theme_bw()


## elevation
b_elev_pred_sum <- data.frame(zhfi_300 = mean(bob_dat_sum$zhfi_300),
                              ztrail_dist = mean(bob_dat_sum$ztrail_dist),
                              zcanopy_0 = mean(bob_dat_sum$zcanopy_0),
                              zelev_1000 = seq(min(bob_dat_sum$zelev_1000),
                                               max(bob_dat_sum$zelev_1000),
                                               length.out=40),
                              log_sl_ = mean(bob_dat_sum$log_sl_),
                              cos_ta_ = mean(bob_dat_sum$cos_ta_),
                              ID = NA,
                              step_id_ = NA,
                              w = NA)

y <- predict(bob_ssf_sum, newdata = b_elev_pred_sum,
             type = "link", se = TRUE)

b_elev_pred_sum$y <- exp(as.vector(y$fit))
b_elev_pred_sum$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
b_elev_pred_sum$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

b_rss_elev_sum <- ggplot(b_elev_pred_sum, aes(zelev_1000)) +
  geom_line(aes(y=y), color="cornflowerblue", lwd = 1.2) +
  geom_line(aes(y=ciu), color="cornflowerblue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="cornflowerblue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Elevation", y = "", title = "Bobcat Summer") +
  theme_bw()

# combine plots
plot_grid(b_rss_hfi_sum, b_rss_trail_sum, b_rss_can_sum, b_rss_elev_sum, ncol = 4)
# 1300x350

# save predicted data
saveRDS(b_hfi_pred_sum,
        "outputs/b_hfi_pred_sum.rds")
saveRDS(b_trail_pred_sum,
        "outputs/b_trail_pred_sum.rds")
saveRDS(b_can_pred_sum,
        "outputs/b_can_pred_sum.rds")
saveRDS(b_elev_pred_sum,
        "outputs/b_elev_pred_sum.rds")

    #### Bobcat Winter ----
b_hfi_pred_win <- readRDS("outputs/b_hfi_pred_win.rds")
b_trail_pred_win <- readRDS("outputs/b_trail_pred_win.rds")
b_can_pred_win <- readRDS("outputs/b_can_pred_win.rds")
b_elev_pred_win <- readRDS("outputs/b_elev_pred_win.rds")
## HFI
b_hfi_pred_win <- data.frame(zhfi_600 = seq(min(bob_dat_win$zhfi_600),
                                            max(bob_dat_win$zhfi_600),
                                            length.out=40),
                             ztrail_dist = mean(bob_dat_win$ztrail_dist),
                             zcanopy_0 = mean(bob_dat_win$zcanopy_0),
                             zelev_100 = mean(bob_dat_win$zelev_100),
                             log_sl_ = mean(bob_dat_win$log_sl_),
                             cos_ta_ = mean(bob_dat_win$cos_ta_),
                             ID = NA,
                             step_id_ = NA,
                             w = NA)

y <- predict(bob_ssf_win, newdata = b_hfi_pred_win,
             type = "link", se = TRUE)

b_hfi_pred_win$y <- exp(as.vector(y$fit))
b_hfi_pred_win$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
b_hfi_pred_win$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

b_rss_hfi_win <- ggplot(b_hfi_pred_win, aes(zhfi_600)) +
  geom_line(aes(y=y), color="#771818", lwd = 1.2) +
  geom_line(aes(y=ciu), color="#771818", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#771818", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Human Footprint", y = "RSS", title = "Bobcat Winter") +
  theme_bw()

## trail
b_trail_pred_win <- data.frame(zhfi_600 = mean(bob_dat_win$zhfi_600),
                               ztrail_dist = seq(min(bob_dat_win$ztrail_dist),
                                                 max(bob_dat_win$ztrail_dist),
                                                 length.out=40),
                               zcanopy_0 = mean(bob_dat_win$zcanopy_0),
                               zelev_100 = mean(bob_dat_win$zelev_100),
                               log_sl_ = mean(bob_dat_win$log_sl_),
                               cos_ta_ = mean(bob_dat_win$cos_ta_),
                               ID = NA,
                               step_id_ = NA,
                               w = NA)

y <- predict(bob_ssf_win, newdata = b_trail_pred_win,
             type = "link", se = TRUE)

b_trail_pred_win$y <- exp(as.vector(y$fit))
b_trail_pred_win$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
b_trail_pred_win$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

b_rss_trail_win <- ggplot(b_trail_pred_win, aes(ztrail_dist)) +
  geom_line(aes(y=y), color="#55371b", lwd = 1.2) +
  geom_line(aes(y=ciu), color="#55371b", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#55371b", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,.0012) +
  labs(x = "Trail Proximity", y = "", title = "Bobcat Winter") +
  theme_bw()


## canopy
b_can_pred_win <- data.frame(zhfi_600 = mean(bob_dat_win$zhfi_600),
                             ztrail_dist = mean(bob_dat_win$ztrail_dist),
                             zcanopy_0 = seq(min(bob_dat_win$zcanopy_0),
                                             max(bob_dat_win$zcanopy_0),
                                             length.out=40),
                             zelev_100 = mean(bob_dat_win$zelev_100),
                             log_sl_ = mean(bob_dat_win$log_sl_),
                             cos_ta_ = mean(bob_dat_win$cos_ta_),
                             ID = NA,
                             step_id_ = NA,
                             w = NA)

y <- predict(bob_ssf_win, newdata = b_can_pred_win,
             type = "link", se = TRUE)

b_can_pred_win$y <- exp(as.vector(y$fit))
b_can_pred_win$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
b_can_pred_win$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

b_rss_can_win <- ggplot(b_can_pred_win, aes(zcanopy_0)) +
  geom_line(aes(y=y), color="#334e33", lwd = 1.2) +
  geom_line(aes(y=ciu), color="#334e33", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#334e33", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,0.0012) +
  labs(x = "% Canopy Cover", y = "", title = "Bobcat Winter") +
  theme_bw()


## elevation
b_elev_pred_win <- data.frame(zhfi_600 = mean(bob_dat_win$zhfi_600),
                              ztrail_dist = mean(bob_dat_win$ztrail_dist),
                              zcanopy_0 = mean(bob_dat_win$zcanopy_0),
                              zelev_100 = seq(min(bob_dat_win$zelev_100),
                                              max(bob_dat_win$zelev_100),
                                              length.out=40),
                              log_sl_ = mean(bob_dat_win$log_sl_),
                              cos_ta_ = mean(bob_dat_win$cos_ta_),
                              ID = NA,
                              step_id_ = NA,
                              w = NA)

y <- predict(bob_ssf_win, newdata = b_elev_pred_win,
             type = "link", se = TRUE)

b_elev_pred_win$y <- exp(as.vector(y$fit))
b_elev_pred_win$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
b_elev_pred_win$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

b_rss_elev_win <- ggplot(b_elev_pred_win, aes(zelev_100)) +
  geom_line(aes(y=y), color="darkblue", lwd = 1.2) +
  geom_line(aes(y=ciu), color="darkblue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="darkblue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,0.5) +
  labs(x = "Elevation", y = "", title = "Bobcat Winter") +
  theme_bw()

plot_grid(b_rss_hfi_win, b_rss_trail_win, b_rss_can_win, b_rss_elev_win, ncol = 4)
# 1300x350

saveRDS(b_hfi_pred_win,
        "outputs/b_hfi_pred_win.rds")
saveRDS(b_trail_pred_win,
        "outputs/b_trail_pred_win.rds")
saveRDS(b_can_pred_win,
        "outputs/b_can_pred_win.rds")
saveRDS(b_elev_pred_win,
        "outputs/b_elev_pred_win.rds")
    #### Coyote ----
c_hfi_pred <- readRDS("outputs/c_hfi_pred.rds")
c_trail_pred <- readRDS("outputs/c_trail_pred.rds")
c_can_pred <- readRDS("outputs/c_can_pred.rds")
c_elev_pred <- readRDS("outputs/c_elev_pred.rds")
## HFI
c_hfi_pred <- data.frame(zhfi_800 = seq(min(coy_dat_all$zhfi_800),
                                        max(coy_dat_all$zhfi_800),
                                        length.out=40),
                         ztrail_dist = mean(coy_dat_all$ztrail_dist),
                         zcanopy_0 = mean(coy_dat_all$zcanopy_0),
                         zelev_600 = mean(coy_dat_all$zelev_600),
                         log_sl_ = mean(coy_dat_all$log_sl_),
                         cos_ta_ = mean(coy_dat_all$cos_ta_),
                         ID = NA,
                         step_id_ = NA,
                         w = NA)


y <- predict(coy_ssf, newdata = c_hfi_pred,
             type = "link", se = TRUE)

c_hfi_pred$y <- exp(as.vector(y$fit))
c_hfi_pred$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
c_hfi_pred$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

c_rss_hfi <- ggplot(c_hfi_pred, aes(zhfi_800)) +
  geom_line(aes(y=y), color="firebrick", lwd = 1.2) +
  geom_line(aes(y=ciu), color="firebrick", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="firebrick", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Human Footprint", y = "RSS", title = "Coyote") +
  theme_bw()


## trail
c_trail_pred <- data.frame(zhfi_800 = mean(coy_dat_all$zhfi_800),
                           ztrail_dist = seq(min(coy_dat_all$ztrail_dist),
                                             max(coy_dat_all$ztrail_dist),
                                             length.out=40),
                           zcanopy_0 = mean(coy_dat_all$zcanopy_0),
                           zelev_600 = mean(coy_dat_all$zelev_600),
                           log_sl_ = mean(coy_dat_all$log_sl_),
                           cos_ta_ = mean(coy_dat_all$cos_ta_),
                           ID = NA,
                           step_id_ = NA,
                           w = NA)

y <- predict(coy_ssf, newdata = c_trail_pred,
             type = "link", se = TRUE)

c_trail_pred$y <- exp(as.vector(y$fit))
c_trail_pred$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
c_trail_pred$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

c_rss_trail <- ggplot(c_trail_pred, aes(ztrail_dist)) +
  geom_line(aes(y=y), color="tan4", lwd = 1.2) +
  geom_line(aes(y=ciu), color="tan4", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="tan4", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,1e+30) +
  labs(x = "Trail Proximity", y = "", title = "Coyote") +
  theme_bw()


## canopy
c_can_pred <- data.frame(zhfi_800 = mean(coy_dat_all$zhfi_800),
                         ztrail_dist = mean(coy_dat_all$ztrail_dist),
                         zcanopy_0 = seq(min(coy_dat_all$zcanopy_0),
                                         max(coy_dat_all$zcanopy_0),
                                         length.out=40),
                         zelev_600 = mean(coy_dat_all$zelev_600),
                         log_sl_ = mean(coy_dat_all$log_sl_),
                         cos_ta_ = mean(coy_dat_all$cos_ta_),
                         ID = NA,
                         step_id_ = NA,
                         w = NA)

y <- predict(coy_ssf, newdata = c_can_pred,
             type = "link", se = TRUE)

c_can_pred$y <- exp(as.vector(y$fit))
c_can_pred$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
c_can_pred$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

c_rss_can <- ggplot(c_can_pred, aes(zcanopy_0)) +
  geom_line(aes(y=y), color="darkgreen", lwd = 1.2) +
  geom_line(aes(y=ciu), color="darkgreen", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="darkgreen", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "% Canopy Cover", y = "", title = "Coyote") +
  theme_bw()


## elevation
c_elev_pred <- data.frame(zhfi_800 = mean(coy_dat_all$zhfi_800),
                          ztrail_dist = mean(coy_dat_all$ztrail_dist),
                          zcanopy_0 = mean(coy_dat_all$zcanopy_0),
                          zelev_600 = seq(min(coy_dat_all$zelev_600),
                                          max(coy_dat_all$zelev_600),
                                          length.out=40),
                          log_sl_ = mean(coy_dat_all$log_sl_),
                          cos_ta_ = mean(coy_dat_all$cos_ta_),
                          ID = NA,
                          step_id_ = NA,
                          w = NA)

y <- predict(coy_ssf, newdata = c_elev_pred,
             type = "link", se = TRUE)

c_elev_pred$y <- exp(as.vector(y$fit))
c_elev_pred$ciu <- exp(as.vector(y$fit + 1.96 * y$se.fit))
c_elev_pred$cil <- exp(as.vector(y$fit - 1.96 * y$se.fit))

c_rss_elev <- ggplot(c_elev_pred, aes(zelev_600)) +
  geom_line(aes(y=y), color="blue", lwd = 1.2) +
  geom_line(aes(y=ciu), color="blue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="blue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,2) +
  labs(x = "Elevation", y = "", title = "Coyote") +
  theme_bw()

plot_grid(c_rss_hfi, c_rss_trail, c_rss_can, c_rss_elev, ncol = 4)
# 1300x350

saveRDS(c_hfi_pred,
        "outputs/c_hfi_pred.rds")
saveRDS(c_trail_pred,
        "outputs/c_trail_pred.rds")
saveRDS(c_can_pred,
        "outputs/c_can_pred.rds")
saveRDS(c_elev_pred,
        "outputs/c_elev_pred.rds")


  ### Individual slope estimates ----
    #### Bobcat Non-Winter ----
bob_ranef_sum <- as.data.frame(coef(bob_ssf_sum)$cond$ID)

bob_ranef_sum_long <- bob_ranef_sum %>% 
  dplyr::select(ztrail_dist,
                zhfi_300,"I(zhfi_300^2)",
                "zhfi_300:ztrail_dist",
                zcanopy_0,"I(zcanopy_0^2)",
                zelev_1000) %>% 
  pivot_longer(cols = c(ztrail_dist,
                        zhfi_300,"I(zhfi_300^2)",
                        "zhfi_300:ztrail_dist",
                        zcanopy_0,"I(zcanopy_0^2)",
                        zelev_1000),
               names_to = "Cov", values_to = "Slope_Estimate")

ggplot(bob_ranef_sum_long, aes(Cov, Slope_Estimate, fill = Cov))+
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.color = "red", outlier.shape = 1, outlier.size = 3) +
  geom_jitter(width = 0) +
  ylim(-3,8) +
  theme_bw() +
  labs(x = "Covariate", y = "Slope Estimate", title = "Variation in Bobcat Slopes (Summer)") +
  scale_x_discrete(labels=c(expression("Canopy"^2),expression("HFI"^2),
                            'Canopy','Elev','HFI','TrailDist*HFI','TrailDist')) +
  scale_fill_brewer(palette = "YlOrRd") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "")

ggsave("plots/slopes_bob_summer.png")


    #### Bobcat Winter ----
bob_ranef_win <- as.data.frame(coef(bob_ssf_win)$cond$ID)

bob_ranef_win_long <- bob_ranef_win %>% 
  dplyr::select(ztrail_dist,"I(ztrail_dist^2)",
                zhfi_600,"I(zhfi_600^2)",
                "zhfi_600:ztrail_dist",
                zcanopy_0,"I(zcanopy_0^2)",
                zelev_100) %>% 
  pivot_longer(cols = c(ztrail_dist,"I(ztrail_dist^2)",
                        zhfi_600,"I(zhfi_600^2)",
                        "zhfi_600:ztrail_dist",
                        zcanopy_0,"I(zcanopy_0^2)",
                        zelev_100),
               names_to = "Cov", values_to = "Slope_Estimate")

ggplot(bob_ranef_win_long, aes(Cov, Slope_Estimate, fill = Cov))+
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.color = "red", outlier.shape = 1, outlier.size = 3) +
  geom_jitter(width = 0) +
  ylim(-3,8) +
  theme_bw() +
  labs(x = "Covariate", y = "Slope Estimate", title = "Variation in Bobcat Slopes (Winter)") +
  scale_x_discrete(labels=c(expression("Canopy"^2),expression("HFI"^2),expression("TrailDist"^2),
                            'Canopy','Elev','HFI','TrailDist*HFI','TrailDist')) +
  scale_fill_brewer(palette = "Blues") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "")

ggsave("plots/slopes_bob_winter.png")


    #### Coyote ----
coy_ranef <- as.data.frame(coef(coy_ssf)$cond$ID)

coy_ranef_long <- coy_ranef %>% 
  dplyr::select(ztrail_dist,"I(ztrail_dist^2)",
                zhfi_800,"I(zhfi_800^2)",
                "zhfi_800:ztrail_dist",
                zcanopy_0,"I(zcanopy_0^2)",
                zelev_600) %>% 
  pivot_longer(cols = c(ztrail_dist,
                        zhfi_800,"I(zhfi_800^2)",
                        "zhfi_800:ztrail_dist",
                        zcanopy_0,"I(zcanopy_0^2)",
                        zelev_600),
               names_to = "Cov", values_to = "Slope_Estimate")

ggplot(coy_ranef_long, aes(Cov, Slope_Estimate, fill = Cov))+
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.color = "red", outlier.shape = 1, outlier.size = 3) +
  geom_jitter(width = 0) +
  ylim(-8,20) +
  theme_bw() +
  labs(x = "Covariate", y = "Slope Estimate", title = "Variation in Coyote Slopes") +
  scale_x_discrete(labels=c(expression("Canopy"^2),expression("HFI"^2),
                            'Canopy','Elev','HFI','TrailDist*HFI','TrailDist')) +
  scale_fill_brewer(palette = "Spectral") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "")

coy_ranef_long <- coy_ranef %>% 
  dplyr::select(ztrail_dist,"I(ztrail_dist^2)") %>% 
  pivot_longer(cols = c(ztrail_dist,"I(ztrail_dist^2)"),
               names_to = "Cov", values_to = "Slope_Estimate") %>% 
  mutate(Slope_Estimate = exp(Slope_Estimate))

ggplot(coy_ranef_long, aes(Cov, Slope_Estimate, fill = Cov))+
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.color = "red", outlier.shape = 1, outlier.size = 3) +
  geom_jitter(width = 0) +
  # ylim(-8,20) +
  theme_bw() +
  labs(x = "", y = "Slope Estimate", title = "Variation in Coyote Slopes") +
  scale_x_discrete(labels=c(expression("Trail Proximity"^2),'Trail Proximity')) +
  scale_fill_manual(values = c("tan4","tan3")) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "")
# 300x500


  ### Histogram of Use ----
    #### Bobcat Non-Winter ----
bob_dat_sum_used <- bob_dat_sum %>% 
  filter(case_ == "TRUE")

hist_trail_sum <- ggplot(bob_dat_sum_used, aes(trail_dist)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="#99795b", color = "lightgrey") + 
  geom_vline(xintercept = mean(bob_dat_sum_used$trail_dist), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Proximity to Trail (m)",
       y = "") +
  theme_bw()

hist_hfi_sum <- ggplot(bob_dat_sum_used, aes(hfi_300)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="#de2c2c", color = "lightgrey") + 
  geom_vline(xintercept = mean(bob_dat_sum_used$hfi_300), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Human Footprint Index",
       y = "") +
  theme_bw()

hist_elev_sum <- ggplot(bob_dat_sum_used, aes(elev_1000)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="cornflowerblue", color = "lightgrey") + 
  geom_vline(xintercept = mean(bob_dat_sum_used$elev_1000), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Elevation (m)",
       y = "") +
  theme_bw()

hist_can_sum <- ggplot(bob_dat_sum_used, aes(canopy_0)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="#53c153", color = "lightgrey") + 
  geom_vline(xintercept = mean(bob_dat_sum_used$canopy_0), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Canopy Cover (%)",
       y = "") +
  theme_bw()

plot_grid(hist_hfi_sum, hist_trail_sum, hist_can_sum, hist_elev_sum, ncol = 4)
# 1300x350

    #### Bobcat Winter ----
bob_dat_win_used <- bob_dat_win %>% 
  filter(case_ == "TRUE")

hist_trail_win <- ggplot(bob_dat_win_used, aes(trail_dist)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="#55371b", color = "lightgrey") + 
  geom_vline(xintercept = mean(bob_dat_win_used$trail_dist), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Proximity to Trail (m)",
       y = "") +
  theme_bw()

hist_hfi_win <- ggplot(bob_dat_win_used, aes(hfi_600)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="#771818", color = "lightgrey") + 
  geom_vline(xintercept = mean(bob_dat_win_used$hfi_600), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Human Footprint Index",
       y = "") +
  theme_bw()

hist_elev_win <- ggplot(bob_dat_win_used, aes(elev_100)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="darkblue", color = "lightgrey") + 
  geom_vline(xintercept = mean(bob_dat_win_used$elev_100), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Elevation (m)",
       y = "") +
  theme_bw()

hist_can_win <- ggplot(bob_dat_win_used, aes(canopy_0)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="#334e33", color = "lightgrey") + 
  geom_vline(xintercept = mean(bob_dat_win_used$canopy_0), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Canopy Cover (%)",
       y = "") +
  theme_bw()


plot_grid(hist_hfi_win, hist_trail_win, hist_can_win, hist_elev_win, ncol = 4)
# 1300x350

    #### Coyote ----
coy_dat_used <- coy_dat_all %>% 
  filter(case_ == "TRUE")

hist_trail <- ggplot(coy_dat_used, aes(trail_dist)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="tan4", color = "lightgrey") + 
  geom_vline(xintercept = mean(coy_dat_used$trail_dist), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Proximity to Trail (m)",
       y = "") +
  theme_bw()

hist_hfi <- ggplot(coy_dat_used, aes(hfi_800)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="firebrick", color = "lightgrey") + 
  geom_vline(xintercept = mean(coy_dat_used$hfi_800), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Human Footprint Index",
       y = "") +
  theme_bw()

hist_elev <- ggplot(coy_dat_used, aes(elev_600)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="blue", color = "lightgrey") + 
  geom_vline(xintercept = mean(coy_dat_used$elev_600), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Elevation (m)",
       y = "") +
  theme_bw()

hist_can <- ggplot(coy_dat_used, aes(canopy_0)) +
  geom_histogram(bins = 20, na.rm = T, 
                 fill="darkgreen", color = "lightgrey") + 
  geom_vline(xintercept = mean(coy_dat_used$canopy_0), 
             color = "goldenrod1", lwd = 1) + 
  labs(x = "Canopy Cover (%)",
       y = "") +
  theme_bw()

plot_grid(hist_hfi, hist_trail, hist_can, hist_elev, ncol = 4)
# 1300x350


  ### Raster Maps ----
plot(trail, col = rev(map.pal("sepia")))
plot(hfi$hfi_0, col = rev(paletteer_c("ggthemes::Classic Red-Black",60)))
plot(canopy$Canopy_Cover, col = rev(map.pal("grass")))
plot(elev$elev_0, col = paletteer_c("ggthemes::Classic Blue", 60))
#700x500

  ### Species GPS Maps ----
meso <- readRDS("outputs/Okanogan_GPS_locations.rds")

bob_vect <- meso %>% 
  filter(Species == "BOB") %>% 
  vect(geom=c("Long", "Lat"), crs="EPSG:4326")  # using a different crs than data was analyzed in to match crs of trail vector
coy_vect <- meso %>% 
  filter(Species == "COY") %>% 
  vect(geom=c("Long", "Lat"), crs="EPSG:4326")

trail_vect <- vect("data_raw/Trails/Trails.shp")

hfi <- rast("data_raw/Prugh_provided_files/human_footprint_2019.tif") %>% 
  crop(trail_vect)

# bobcat
plot(hfi, col = rev(paletteer_c("ggthemes::Classic Red-Black",60)),
     main = "Bobcat Human Footprint")
points(bob_vect, col = "darkolivegreen4", cex = 0.3, alpha = 0.7)
lines(trail_vect, lwd = 0.4, col = "orange")

plot.new()
legend("center", legend = c("Bobcat Step", "Trail"),
       col = c("darkolivegreen4","orange"), pch = c(20,NA), lty = c(NA,1),
       cex = 0.7, pt.cex = 1, bty = "n")

# coyote
plot(hfi, col = rev(paletteer_c("ggthemes::Classic Red-Black",60)),
     main = "Coyote HFI", legend = F)
points(coy_vect, col = "cadetblue3", cex = 0.3, alpha = 0.7)
lines(trail_vect, lwd = 0.4, col = "orange")
#500x550

plot.new()
legend("center", legend = c("Coyote Step", "Trail"),
       col = c("cadetblue3","orange"), pch = c(20,NA), lty = c(NA,1),
       cex = 0.7, pt.cex = 1, bty = "n")
