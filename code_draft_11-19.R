### 2024-11-19 Draft
## Final Project Code

library(tidyverse)
library(terra)

#### Tidy GPS Data ####

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

#### Tidy Covariate Data ####

## Hiking Trails 
trail <- rast("rasters/trail_raster.tif") %>% 
  # reproject to appropriate crs
  project("EPSG:32610", method = "bilinear")

# make data binary
trail[!is.na(trails)] <- 1  # change any values to 1

# save raster
writeRaster(trail, "outputs/trail.tif",
            overwrite = T)

# give data meaningful name
names(trail) <- "trail_dist"

# Calculate Euclidean distance to trails
trail_dist <- distance(trails,
                            filename = "outputs/trail_dist.tif",
                            filetype = "GTiff", gdal = c("COMPRESS=DEFLATE", "BIGTIFF=YES"),
                            overwrite = T)



## Human Footprint Index (2019)
hfi <- rast("data_raw/Prugh_provided_files/human_footprint_2019.tif") %>% 
  project("EPSG:32610", method = "bilinear") %>% 
  # crop to trails extent
  crop(trails)

names(hfi) <- "hfi_0"

# save as new file
writeRaster(hfi, "outputs/hfi.tif",
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

writeRaster(DEM, "outputs/elevation.tif",
            overwrite = T)



## Percent Tree Canopy Cover (2019)
canopy <- rast("rasters/NLCD_CanopyCover.tiff") %>%
  project("EPSG:32610", method = "bilinear") %>% 
  crop(trails)

names(canopy) <- "canopy_0"

writeRaster(canopy, "outputs/canopy.tif",
            overwrite = T)

#### Calculate Available Step Lengths and Turn Angles ####

rm(list=ls())

library(amt)


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

#### Create Raster Buffers ####
rm(list=ls())

trail <- rast("outputs/distance_to_trail.tif")
elev <- rast("outputs/elevation.tif")
canopy <- rast("outputs/canopy.tif")
hfi <- rast("outputs/hfi.tif")
bobcat_end <- readRDS("outputs/bobcat_steps.rds") %>%
  select(x2_, y2_, case_binary)


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
writeRaster(elev, "outputs/elevation_buffers.tif",
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

writeRaster(hfi, "outputs/hfi_buffers.tif",
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

writeRaster(trail, "outputs/trail_buffers.tif",
            overwrite = T)


  ## Percent Canopy Cover 
# change values over 87 to NA, these are an error from reprojection
canopy[canopy>87] <- NA

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

writeRaster(canopy, "outputs/canopy_buffers.tif",
            overwrite = T)



#### Extract Covariates ####
rm(list=ls())

library(raster)

# load previously created files
trail <- rast("outputs/trail_dist.tif")
elev <- rast("outputs/elevation_buffers.tif")
canopy <- rast("outputs/canopy_buffers.tif")
hfi <- rast("outputs/hfi_buffers.tif")
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


# Bobcat: determine scales of effect
  ## hfi
# create data frame to hold AIC values for comparison
hfi_scales_b <- data.frame(covariate = "hfi",
                            scale = seq(0,1000,by = 100),
                            AIC = NA,
                           delta_AIC = NA)

# for loop to calculate AICs and enter in df
for(i in 1:nrow(hfi_scales_b)){
  cov_i <- paste0(hfi_scales_b$covariate[i],"_",hfi_scales_b$scale[i])
  data_i <- bobcat_hfi[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_hfi$w)
  hfi_scales_b$AIC[i] <- AIC(model_i)
}
# for loop to calculate delta AIC values
for (i in 1:length(hfi_scales_b)) {
  min_i <- min(hfi_scales_b$AIC)
  hfi_scales_b$delta_AIC <- hfi_scales_b$AIC - min_i
}

# lowest AIC value
hfi_scales_b[which(hfi_scales_b$AIC==min(hfi_scales_b$AIC)),]
# hfi_800


  ## canopy
canopy_scales_b <- data.frame(covariate = "canopy",
                         scale = seq(0,1000,by = 100),
                         AIC = NA,
                         delta_AIC = NA)

for(i in 1:nrow(canopy_scales_b)){
  cov_i <- paste0(canopy_scales_b$covariate[i],"_",canopy_scales_b$scale[i])
  data_i <- bobcat_canopy[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_canopy$w)
  canopy_scales_b$AIC[i] <- AIC(model_i)
}
for (i in 1:length(canopy_scales_b)) {
  min_i <- min(canopy_scales_b$AIC)
  canopy_scales_b$delta_AIC <- canopy_scales_b$AIC - min_i
}

canopy_scales_b[which(canopy_scales_b$AIC==min(canopy_scales_b$AIC)),]
# canopy_0


    ## elevation
elev_scales_b <- data.frame(covariate = "elev",
                            scale = seq(0,1000,by = 100),
                            AIC = NA,
                            delta_AIC = NA)

for(i in 1:nrow(elev_scales_b)){
  cov_i <- paste0(elev_scales_b$covariate[i],"_",elev_scales_b$scale[i])
  data_i <- bobcat_elev[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = bobcat_elev$w)
  elev_scales_b$AIC[i] <- AIC(model_i)
}
for (i in 1:length(elev_scales_b)) {
  min_i <- min(elev_scales_b$AIC)
  elev_scales_b$delta_AIC <- elev_scales_b$AIC - min_i
}

elev_scales_b[which(elev_scales_b$AIC==min(elev_scales_b$AIC)),]
# elev_200


# bind scales of effect for all covariates and z-score standardize
hfi_b <- dplyr::select(bobcat_hfi, hfi_800) %>% 
  mutate(zhfi_800 = scale(hfi_800))
canopy_b <- dplyr::select(bobcat_canopy, canopy_0) %>%
  mutate(zcanopy_0 = scale(canopy_0))
elev_b <- dplyr::select(.data = bobcat_elev, elev_200) %>%
  mutate(zelev_200 = scale(elev_200))
trail_b <- raster::extract(trail, bobcat_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  mutate(ztrail_dist = scale(trail_dist))

bobcat_extract1 <- cbind(hfi_b, canopy_b, elev_b, trail_b, 
                         bobcat_extract)


# plot scales of effect
par(mfrow=c(1,3))
plot(delta_AIC ~ scale, hfi_scales_b,
     type = "b",
     xlab = "",
     ylab = expression(paste(Delta, " AIC")),
     main = "Human Footprint Index")
plot(delta_AIC ~ scale, canopy_scales_b,
     type = "b",
     xlab = "Scale (buffer radius in meters)",
     ylab = "",
     main = "Percent Canopy Cover")
plot(delta_AIC ~ scale, elev_scales_b,
     type = "b",
     xlab = "",
     ylab = "",
     main = "Elevation")

mtext(expression(bold("Bobcat Scales of Effect")),
      at=-1000, line=3,
      cex=1)
par(mfrow=c(1,1))

save.image("outputs/bobcat_sof.png")

# Coyote: determine scales of effect
  ## hfi
hfi_scales_c <- data.frame(covariate = "hfi",
                           scale = seq(0,1000,by = 100),
                           AIC = NA,
                           delta_AIC = NA)

for(i in 1:nrow(hfi_scales_c)){
  cov_i <- paste0(hfi_scales_c$covariate[i],"_",hfi_scales_c$scale[i])
  data_i <- coyote_hfi[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = coyote_hfi$w)
  hfi_scales_c$AIC[i] <- AIC(model_i)
}
for (i in 1:length(hfi_scales_c)) {
  min_i <- min(hfi_scales_c$AIC)
  hfi_scales_c$delta_AIC <- hfi_scales_c$AIC - min_i
}

hfi_scales_c[which(hfi_scales_c$AIC==min(hfi_scales_c$AIC)),]
# hfi_800


  ## canopy
canopy_scales_c <- data.frame(covariate = "canopy",
                              scale = seq(0,1000,by = 100),
                              AIC = NA,
                              delta_AIC = NA)

for(i in 1:nrow(canopy_scales_c)){
  cov_i <- paste0(canopy_scales_c$covariate[i],"_",canopy_scales_c$scale[i])
  data_i <- coyote_canopy[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = coyote_canopy$w)
  canopy_scales_c$AIC[i] <- AIC(model_i)
}
for (i in 1:length(canopy_scales_c)) {
  min_i <- min(canopy_scales_c$AIC)
  canopy_scales_c$delta_AIC <- canopy_scales_c$AIC - min_i
}

canopy_scales_c[which(canopy_scales_c$AIC==min(canopy_scales_c$AIC)),]
# canopy_0


  ## elevation
elev_scales_c <- data.frame(covariate = "elev",
                            scale = seq(0,1000,by = 100),
                            AIC = NA,
                            delta_AIC = NA)

for(i in 1:nrow(elev_scales_c)){
  cov_i <- paste0(elev_scales_c$covariate[i],"_",elev_scales_c$scale[i])
  data_i <- coyote_elev[,c("case_",cov_i)]
  model_i <- glm(data_i[,1] ~ data_i[,2], data = data_i, family = poisson, weights = coyote_elev$w)
  elev_scales_c$AIC[i] <- AIC(model_i)
}
for (i in 1:length(elev_scales_c)) {
  min_i <- min(elev_scales_c$AIC)
  elev_scales_c$delta_AIC <- elev_scales_c$AIC - min_i
}

elev_scales_c[which(elev_scales_c$AIC==min(elev_scales_c$AIC)),]
# elev_600




# bind scales of effect for all covariates and z-score standardize
hfi_c <- dplyr::select(coyote_hfi, hfi_800) %>% 
  mutate(zhfi_800 = scale(hfi_800))
canopy_c <- dplyr::select(coyote_canopy, canopy_0) %>%
  mutate(zcanopy_0 = scale(canopy_0))
elev_c <- dplyr::select(.data = coyote_elev, elev_600) %>%
  mutate(zelev_600 = scale(elev_600))
trail_c <- raster::extract(trail, coyote_extract[,c("x2_","y2_")]) %>%
  dplyr::select(-ID) %>% 
  mutate(ztrail_dist = scale(trail_dist))

coyote_extract1 <- cbind(hfi_c, canopy_c, elev_c, trail_c, 
                         coyote_extract)


# plot scales of effect
par(mfrow=c(1,3))
plot(delta_AIC ~ scale, hfi_scales_c,
     type = "b",
     xlab = "",
     ylab = expression(paste(Delta, " AIC")),
     main = "Human Footprint Index")
plot(delta_AIC ~ scale, canopy_scales_c,
     type = "b",
     xlab = "Scale (buffer radius in meters)",
     ylab = "",
     main = "Percent Canopy Cover")
plot(delta_AIC ~ scale, elev_scales_c,
     type = "b",
     xlab = "",
     ylab = "",
     main = "Elevation")

mtext(expression(bold("Coyote Scales of Effect")),
      at=-1000, line=3,
      cex=1)
par(mfrow=c(1,1))

save.image("outputs/coyote_sof.png")


# final touches before SSF
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
  # calculate sample size for each animal (dividing by 11 because of available points)
  mutate(n = n()/11) %>%
  ungroup()

# save fully processed data for SSF
saveRDS(bobcat_extract2, "outputs/bobcat_covariates.rds")
saveRDS(coyote_extract2, "outputs/coyote_covariates.rds")




#### Fit Models ####
rm(list=ls())

library(glmmTMB)
library(doParallel)

bob_ssf_dat <- readRDS("outputs/bobcat_covariates.rds")  %>%
  mutate(log_sl_ = log(sl_),                                   # fix above and remove eventually
         cos_ta_ = cos(ta_)) %>% 
  # select individuals with more than 100 fixes
  filter(n >= 100)

# set number of cores
nt <- min(detectCores()) - 2

# fit mixed-effects SSF
bob_ssf <-  glmmTMB(case_ ~ zhfi_800 + zcanopy_0 + ztrail_dist + zelev_200 +
                      (1|step_id_) + (0 + zhfi_800 | ID) + (0 + zcanopy_0 | ID) + (0 + ztrail_dist | ID) + (0 + zelev_200 | ID) +
                      log_sl_ + (0 + log_sl_ | ID) + cos_ta_ + (0 + cos_ta_ | ID),
                    family=poisson,
                    weights = w,
                    doFit=T,
                    data = bob_ssf_dat,
                    control = glmmTMBControl(parallel = nt))

summary(bob_ssf)

