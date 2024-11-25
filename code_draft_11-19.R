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
trails <- rast("rasters/trail_raster.tif") %>% 
  # reproject to appropriate crs
  project("EPSG:32610", method = "bilinear")

# make data binary
trails[!is.na(trails)] <- 1  # change any values to 1

# give data meaningful name
names(trails) <- "trails"

# save raster
writeRaster(trails, "outputs/trails.tif",
            overwrite = T)

# Calculate Euclidean distance to trails
trail_dist <- distance(trails,
                            filename = "outputs/distance_from_trail.tif",
                            filetype = "GTiff", gdal = c("COMPRESS=DEFLATE", "BIGTIFF=YES"),
                            overwrite = T)



## Human Footprint Index (2019)
hfi <- rast("data_raw/Prugh_provided_files/human_footprint_2019.tif") %>% 
  project("EPSG:32610", method = "bilinear") %>% 
  # crop to trails extent
  crop(trails)

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
names(DEM) <- "Elevation"

writeRaster(DEM, "outputs/DEM.tif",
            overwrite = T)



## Percent Tree Canopy Cover (2019)
canopy <- rast("rasters/NLCD_CanopyCover.tiff") %>%
  project("EPSG:32610", method = "bilinear") %>% 
  crop(trails)

names(canopy) <- "Canopy_Cover"

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
                     # calculate cos(ta)
                     mutate(cos_ta = cos(ta_)))) %>% 
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
                     mutate(cos_ta = cos(ta_)))) %>% 
  dplyr::select(-data, -trk) %>%
  unnest(cols = stp) %>%
  left_join(distinct(data.frame(meso)[,c("ID", "Sex", "Species")])) %>%
  mutate(case_binary = ifelse(case_ == TRUE, 1, 0))

# save data as RDS
write_rds(bobcat1, "outputs/bobcat_steps.rds")
write_rds(coyote1, "outputs/coyote_steps.rds")

#### Determine Scales of Effect ####
rm(list=ls())

library(doParallel)
library(foreach)

trail <- rast("outputs/distance_from_trail.tif")
elev <- rast("outputs/DEM.tif")
canopy <- rast("outputs/canopy.tif")
hfi <- rast("outputs/hfi.tif")
bobcat_end <- readRDS("outputs/bobcat_steps.rds") %>%
  select(x2_, y2_, case_binary)


    # # extract covariate values at point
    # point_elev <- terra::extract(elev,
    #                              bobcat_end[,c("x2_","y2_")])
    # 
    # head(point_elev)
    # 
    # # add new column to bobcat_end data frame
    # bobcat_end$Elev_0 <- point_elev$Elevation
    # 


######
# set buffer range
buffer_sizes <- seq(100,300,by=100)


# buffer raster maps ####
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

#SAVE
writeRaster(elev, "outputs/elevation_buffers.tif",
            overwrite = T)


## hfi
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

## trail
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


## canopy #### TRY CANOPY FIRST, it has the lowest resolution ####
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


### Need to z-score standardize here ###

rm(list=ls())

library(geosphere)
library(raster)

# load previously created files
trail <- rast("outputs/distance_from_trail.tif")
elev <- rast("outputs/DEM.tif")
canopy <- rast("outputs/canopy.tif")
hfi <- rast("outputs/hfi.tif")
bobcat_extract <- readRDS("outputs/bobcat_steps.rds")
coyote2 <- readRDS("outputs/coyote_steps.rds")


# Change step lengths from degrees to meters, define season
bobcat_extract <- bobcat_extract %>% 
  dplyr::select(-sl_) %>%
  mutate(# step lengths from amt are in degrees; overwrite to meters
    sl_ = distGeo(.[,c("x1_","y1_")],.[,c("x2_","y2_")])) %>% 
  mutate(log_sl_ = log(sl_),
         # Define season
         season = ifelse(month(t2_) %in% 4:11, "summer", "winter"))



##### Trying to extract covariates #####

bobcat_extracted1 <- data.frame()
for(j in unique(bobcat_extract$season)) { # for each season within each year
  print(j)
  dat <- bobcat_extract %>% filter(season == j) %>% 
    # % forest cover
    mutate(hfi = raster::extract(hfi, .[,c("x2_","y2_")]))
  bobcat_extracted1 <- rbind(bobcat_extracted1, dat)
}



bob_test <- raster::extract(hfi, bobcat_extract[,c("x2_","y2_")])



# Extract predictors
bobcat2 <- data.frame()
for(j in unique(bobcat1$season)) { # for each season
  print(j)
  dat <- bobcat1 %>% filter(season == j) %>% 
    # forest cover
    mutate(canopy = raster::extract(canopy, .[,c("x2_","y2_")]),
           human_footprint = raster::extract(hfi, .[,c("x2_","y2_")]),
           elevation = raster::extract(elev, .[,c("x2_","y2_")]))
  bobcat2 <- rbind(bobcat2, dat)
}

coy_stp1 <- data.frame()
for(j in unique(coy_stp$season)) { # for each season
  print(j)
  dat <- coy_stp %>% filter(season == j) %>% 
    # forest cover
    mutate(canopy = raster::extract(canopy, .[,c("x2_","y2_")]),
           human_footprint = raster::extract(hfi, .[,c("x2_","y2_")]),
           elevation = raster::extract(elev, .[,c("x2_","y2_")]))
  coy_stp1 <- rbind(coy_stp1, dat)
}


# read fully pre-processed data for SSF
bob_stp2 <- bob_stp1 %>% 
  # create unique step ID for each animal
  mutate(step_id_ = paste(ID, step_id_, sep = "_"),
         human_footprint = human_footprint,
         canopy = canopy,
         elevation = elevation) %>%
  group_by(ID) %>% 
  # calculate sample size for each animal (dividing by 11 because of available points)
  mutate(n = n()/11) %>%
  ungroup()

coy_stp2 <- coy_stp1 %>% 
  # create unique step ID for each animal
  mutate(step_id_ = paste(ID, step_id_, sep = "_"),
         human_footprint = human_footprint,
         canopy = canopy,
         elevation = elevation) %>%
  group_by(ID) %>% 
  # calculate sample size for each animal (dividing by 11 because of available points)
  mutate(n = n()/11) %>%
  ungroup()

# Save fully processed data for SSF
saveRDS(bob_stp2, "outputs/bob_stp2.rds")
saveRDS(coy_stp2, "outputs/coy_stp2.rds")

#### Fit Models ####

library(parallel)
library(glmmTMB)

ncore <- min(parallel::detectCores()) - 2

# read fully processed data for SSFs
coy_ssf_dat <- readRDS("outputs/coy_stp2.rds")  %>% 
  # select coyotes with more than 100 fixes
  filter(n >= 100)



# fit SSF using glmmTB
glmmTMB_coy <-  glmmTMB(case_ ~ human_footprint + I(human_footprint^2) + (0 + human_footprint | ID) + (0 + I(human_footprint^2) | ID) + 
                          canopy + I(canopy^2) + (0 + canopy | ID) + (0 + I(canopy^2) | ID) +
                          elevation + I(elevation^2) + (0 + elevation | ID) + (0 + I(elevation^2) | ID) +
                          # control for step length
                          log_sl + (0 + log_sl | ID) + 
                          # strata
                          (1|step_id_),
                        family=poisson, doFit=T,
                        data = coy_ssf_dat, 
                        # Tell glmmTMB not to change the last standard deviation, all other values are estimated freely
                        map = list(theta = factor(c(1:9, NA))),
                        # Set the value of the standard deviation of the strata (the last ranef) to large constant value
                        start = list(theta = c(rep(0, times = 9),log(1000))),
                        control = glmmTMBControl(parallel = ncore)) 
#saveRDS(glmmTMB_coy, "Analysis/glmmTMB_coy_20230215_minusToPlusScale.rds")
#glmmTMB_coy <- readRDS("Analysis/glmmTMB_coy_20230215_minusToPlusScale.rds")
summary(glmmTMB_coy)

