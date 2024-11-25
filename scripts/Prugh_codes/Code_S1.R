# Main script for analysis of lethal human shields 
# Code written by C. Cunningham 2023-02-24
# Accompanies the paper: LR Prugh, CX Cunningham, RM Windell, BN Kertson, 
# R Ganz, SL Walker, AJ Wirsing. The paradox of the lethal human shield outside protected areas. 

# The analysis is broken into the following sections:

# 1. Load predictors
# 2. Wolf utilization distributions and resource selection functions
# 3. Cougar utilization distributions and resource selection functions
# 4. Prepare mesopredator locations for step selection function
# 5. Fit coyote step selection functions
# 6. Fit bobcat step selection functions
# 7. Create multipanel figure of selection strength & avg effects plots

# See the accompanying script "Functions for human shield analysis v4.R" for several custom functions 
# for constructing utilization distributions, calculating selection strength and graphing results.

# download some packages not on CRAN
#remotes::install_version("SDMTools", "1.1-221.2") # not available on CRAN anymore
#devtools::install_github("pratikunterwegs/atlastools")
#using development version of 'amt': devtools::install_github("jmsigner/amt", build_vignettes = T, force = T)

# load packages
pacman::p_load(tidyverse, sp, leaflet, amt, readr, dplyr, tidyr, purrr, 
               atlastools, spdplyr, data.table, pals, raster, RStoolbox, sf, cowplot, ggmap,
               plotly, future.apply, geosphere, furrr, GGally, MuMIn, suncalc, survival, 
               raster, magick, rphylopic, parallel, mgcv, glmmTMB, gratia, scales)

# load custom functions
source("CustomFunctionsS2.R")

options(scipen = 999) # turn off scientific notation
options(digits = 15) # set digits to 15 to ensure GPS coordinates aren't truncated
sf_use_s2(FALSE) # workaround for one of the amt functions 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 1.0 Load predictors ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# proportion forest within 250 m radius moving window. This is a stack with layers for 2018, 2019 & 2020
forest_pc <- stack("Spatial layers/Proportional Landcover Rasters/forestmix2prop_latlong_stack.gri")
#raster::plot(forest_pc)

# human footprint index; stack with layers for 2018, 2019 & 2020
human_footprint_scaled <- stack("Spatial layers/Human footprint/human_footprint_2020_scaled")
hfi <- raster("Spatial layers/Human footprint/human_footprint_2019_scaled.tif")
#raster::plot(human_footprint_scaled)

# caclulate mean HFI averaged across both study areas
MV_hfi <- crop(human_footprint_scaled, extent(c(-120.654, -119.568, 47.909, 49.03))) %>%
  values(.) %>% data.frame()
NE_hfi <- crop(human_footprint_scaled, extent(c(-118.234, -117.214, 47.762, 48.617))) %>%
  values(.) %>% data.frame()
hfi_both <- bind_rows(MV_hfi, NE_hfi) %>%
  mutate(meanHFI = (human_footprint_2018_scaled + human_footprint_2019_scaled + human_footprint_2020_scaled)/3)
mean(hfi_both$meanHFI)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 2.0 Wolf Utilization Distributions ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Raw location data for wolves and cougars is sensitive and cannot be made publicly available.
# However, on approx line 106 you can read a saved raster stack of the wolf utilization distributions.

# First, run script "Cleaning GPS data for lethal shields analysis.R", 
# which produces wolf_gps_sp2 & wolf_gps_sp3 (which contain the same data in different object classes)
wolf_gps_sp2
wolf_gps_sp3

# Read table specifying what data to use when none available in a given yea see references in "Wolf pack GPS details.xlsx"
lookup_table_UD <- read.csv("Wolf telemetry/lookup_table_UD.csv")

# calculate UDs; see script containing custom functions for details of "wolf_UD_function"
raster_template <- make_trast(wolf_gps_sp2, res = 0.005); raster_template # make raster template so all UDs on same grid and calculate utlization distribution
wolf_UDs <- wolf_UD_function(wolfPackInfo = lookup_table_UD, wolf_gps_data = wolf_gps_sp2, min_locs = 50) 

# create yearly UD rasters, and turn into brick
UD_2018_summer <- sum(brick(wolf_UDs$UD[c(grep(c("2018_summer"), names(wolf_UDs$UD)))]))
UD_2018_winter <- sum(brick(wolf_UDs$UD[c(grep(c("2018_winter"), names(wolf_UDs$UD)))]))
UD_2019_summer <- sum(brick(wolf_UDs$UD[c(grep("2019_summer", names(wolf_UDs$UD)))]))
UD_2019_winter <- sum(brick(wolf_UDs$UD[c(grep("2019_winter", names(wolf_UDs$UD)))]))
UD_2020_summer <- sum(brick(wolf_UDs$UD[c(grep("2020_summer", names(wolf_UDs$UD)))]))
UD_2020_winter <- sum(brick(wolf_UDs$UD[c(grep("2020_winter", names(wolf_UDs$UD)))]))
UD_stack <- brick(UD_2018_summer, UD_2018_winter, UD_2019_summer, UD_2019_winter, UD_2020_summer, UD_2020_winter)
names(UD_stack) <- c("yr_2018_summer", "yr_2018_winter", "yr_2019_summer", "yr_2019_winter", "yr_2020_summer", "yr_2020_winter")

# scaling UD values according to min/max across all years 
minWolf <- min(cellStats(UD_stack, stat = "min"))
maxWolf <- max(quantile(UD_stack, probs = 0.9995)) # (using 0.9995th percentile to slightly down-weight anomalously high pixels)
UD_stack_scale <- ((UD_stack-minWolf)/(maxWolf-minWolf))
for(i in 1:nlayers(UD_stack_scale)) {
  r <- UD_stack_scale[[i]]
  r[r > 1] <- 1
  UD_stack_scale[[i]] <- r - 0.5
}

# Save & read raster stack of wolf UDs
#writeRaster(UD_stack_scale, "Analysis/wolfUDStack")
UD_stack_scale <- stack("Analysis/wolfUDStack")

# plots for sanity checks only, blocking out area with very low wolf UD for easier visualization
raster::plot(UD_stack_scale)
zero_wolf <- UD_stack_scale[[1]]
zero_wolf[zero_wolf < -0.49] <- NA
raster::plot(zero_wolf, colNA = "blue")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 2.1 Wolf Resource Selection Function  ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# using wolf_rsf_dat which was created in "Cleaning GPS data for lethal shields analysis.R" 

# defining study area polygons based on a 100% MCP around all wolf locations at each site
NE_poly <- wolf_rsf_dat %>% filter(studyArea == "NE") %>% make_track(.x = Longitude, .y = Latitude, .t = date_time, crs = 4326) %>% hr_mcp(levels = 1)
MV_poly <- wolf_rsf_dat %>% filter(studyArea == "MV") %>% make_track(.x = Longitude, .y = Latitude, .t = date_time, crs = 4326) %>% hr_mcp(levels = 1)

# for each animal, distribute 10 random points across relevant study site
wolf_rsf_dat1 <- list()
for(i in unique(wolf_rsf_dat$ID)){
  print(i)
  dat <- wolf_rsf_dat[wolf_rsf_dat$ID == i,] # select animal
  study_area <- unique(dat$studyArea) # select that animal's study site
  tracks <- dat %>% make_track(.x = Longitude, .y = Latitude, .t = date_time, crs = 4326) # make track for animal
  # distribute 10 random points for each location
  wolf_rsf_dat1[[i]] <- random_points(get(paste(study_area, "poly", sep = "_")), n = nrow(dat)*10, presence = tracks) %>%
    mutate(ID = i, studyArea = study_area)
}

# extract spatial covariates + some extra formatting
wolf_rsf_dat2 <- wolf_rsf_dat1 %>% bind_rows() %>%
  # extract human footprint
  mutate(human_footprint = raster::extract(hfi, .[,c("x_","y_")]) , 
         case_binary = ifelse(case_ == T, 1, 0), ID = factor(paste(ID))) %>%
  # add sample size to use as filter (dividing by 11 because of 1:10 ratio of used:available)
  group_by(ID) %>% mutate(n = n()/11) %>% ungroup() %>%
  # add large constant weight for available points, as recommended for RSFs by Muff et al (2020) JAE
  mutate(weight = ifelse(case_ == TRUE, 1, 1000)) %>% 
  # select animals with more than 100 locations
  filter(n > 100)

# fit model with glmmTMB, following the mixed effects approach of Muff et al (2020) JAE
glmmTMB_wolf <- glmmTMB(case_ ~ 
                          # fixed effects
                          human_footprint + I(human_footprint^2) + 
                          # random intercepts and slopes
                          (1 | ID) +
                          (0 + human_footprint | ID) + (0 + I(human_footprint^2) | ID),
                        family=binomial(), doFit=T,
                        data = wolf_rsf_dat2, weights = weight) 
#saveRDS(glmmTMB_wolf, "Analysis/saved models/glmmTMB_wolf_20230216.rds")
#glmmTMB_wolf <- readRDS("Analysis/saved models/glmmTMB_wolf_20230216.rds")
summary(glmmTMB_wolf)

# create new data to predict log-RSS
wolf_x1 <- data.frame(
  human_footprint = seq(min(wolf_rsf_dat2$human_footprint), max(wolf_rsf_dat2$human_footprint), by = 0.01), 
  ID = "001X") %>% # ID is needed to use predict function but it will have no effect because we choose population-level predictions
  mutate(weight = 1)
# create dataframe representing the point of comparison, i.e. what selection is relative to.
# here we make selection relative to lowest observed value of HFI, akin to wildland areas
wolf_x2 <- wolf_x1[wolf_x1$human_footprint == min(wolf_x1$human_footprint),]

# predict population-level relative selectio strength
# this is a custom function: see "Functions for human shield analysis.R"
wolf_logRSS <- calc_logRSS_glmmTMB(x1 = wolf_x1, x2 = wolf_x2, model = glmmTMB_wolf, ci_level = 0.95, model_type = "RSF")
#saveRDS(wolf_logRSS, "Analysis/saved models/wolf_logRSS_20230216.rds")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 3.0 Format cougar data ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Raw location data for wolves and cougars is sensitive and cannot be made publicly available.
# However, on approx line 244 you can read a saved raster stack of the wolf utilization distributions.

# load data & format
cougar <- read.csv("Data from Lauren/Cougar_gps_combined.csv") %>% 
  mutate(date_time = as.POSIXct(date_time)) %>%
  group_by(ID) %>%
  arrange(date_time) %>%
  ungroup() %>%
  SpatialPointsDataFrame(.[,c("Long", "Lat")], data = ., proj4string = CRS("EPSG:4326")) %>% 
  mutate(StudyArea = ifelse(grepl("NEC", ID), "NE", "MV"),
         analysis_year = ifelse(month(date_time) == 12, year(date_time) + 1, year(date_time))) %>%
  filter(!is.na(analysis_year)) %>% 
  # omit animals after they dispersed from study area
  filter(ID != "NEC142M" | date_time < as.POSIXct("2019-07-03 00:00:00"), 
         ID != "NEC106M" | date_time < as.POSIXct("2017-05-03 00:00:00"),
         ID != "NEC143F" | date_time < as.POSIXct("2019-04-15 00:00:00"),
         ID != "NEC146M" | date_time < as.POSIXct("2019-06-24 00:00:00"),
         ID != "NEC106M" | date_time < as.POSIXct("2017-05-11 00:00:00"))

# plot cougar locations
leaflet() %>%
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
  addCircles(data = cougar, 
             color = ~ colorFactor(cols25(), domain = factor(ID))(factor(ID)),
             label = paste(cougar$ID,
                           cougar$date_time)) 
# make 'amt' track
cougar1 <- cougar %>%
  make_track(Long,Lat, date_time, id = ID, crs = 4326, all_cols = T) %>%  
  mutate(date_time = t_, 
         season = ifelse(month(t_) %in% 4:11, "summer", "winter"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 3.1 Cougar Utilization Distribution ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# calculate utilization distributions for summer and winter; see functions script for details
raster_template1 <- make_trast(cougar1, res = 0.005) # define common raster template
cougar_UDs_summer <- cougar_UD_function(cougar_gps_data = cougar1 %>% filter(season == "summer"))
names(cougar_UDs_summer$UD) <- paste(names(cougar_UDs_summer$UD), "summer", sep = "_")
cougar_UDs_winter <- cougar_UD_function(cougar_gps_data = cougar1 %>% filter(season == "winter"))
names(cougar_UDs_winter$UD) <- paste(names(cougar_UDs_winter$UD), "winter", sep = "_")

# stack rasters
cougar_UD_stack <- brick(cougar_UDs_summer$UD$year_2018_summer, cougar_UDs_summer$UD$year_2019_summer, cougar_UDs_summer$UD$year_2020_summer,
                         cougar_UDs_winter$UD$year_2018_winter, cougar_UDs_winter$UD$year_2019_winter, cougar_UDs_winter$UD$year_2020_winter)
names(cougar_UD_stack) <- c("yr_2018_summer", "yr_2019_summer", "yr_2020_summer","yr_2018_winter", "yr_2019_winter", "yr_2020_winter")
raster::plot(cougar_UD_stack)

# scaling UDs according to min/max across all years
minCoug <- min(cellStats(cougar_UD_stack, stat = "min"))
maxCoug <- NULL
for(i in nlayers(cougar_UD_stack)) {
  maxCoug <- c(maxCoug, max(quantile(cougar_UD_stack[[i]], probs = 0.999)))}
cougar_UD_stack_scale <- ((cougar_UD_stack-minCoug)/(maxCoug-minCoug))
for(i in 1:nlayers(cougar_UD_stack_scale)) {
  r <- cougar_UD_stack_scale[[i]]
  r[r > 1] <- 1
  cougar_UD_stack_scale[[i]] <- r - 0.5
}

# save and read raster stack of cougar UDs
#writeRaster(cougar_UD_stack_scale, "Analysis/CougarUDStack")
cougar_UD_stack_scale <- stack("Analysis/CougarUDStack")

# sanity checks
raster::plot(cougar_UD_stack_scale)
zero_coug <- cougar_UD_stack_scale$yr_2018_summer
zero_coug[zero_coug < -0.49] <- NA
raster::plot(zero_coug, colNA = "blue")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 3.2 Cougar Resource Selection Function ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# defining study area polygons using MCP around all cougar locations at each site
NE_poly_coug <- cougar %>% filter(StudyArea == "NE") %>% make_track(.x = Long, .y = Lat, .t = date_time, crs = 4326) %>% hr_mcp(levels = 1)
MV_poly_coug <- cougar %>% filter(StudyArea == "MV") %>% make_track(.x = Long, .y = Lat, .t = date_time, crs = 4326) %>% hr_mcp(levels = 1)

# for each animal, distribute 10 random points across study site
coug_rsf_dat1 <- list()
for(i in unique(cougar$ID)){
  print(i)
  dat <- cougar[cougar$ID == i,]
  study_area <- dat$StudyArea[1]
  tracks <- dat %>% make_track(.x = Long, .y = Lat, .t = date_time, crs = 4326) 
  coug_rsf_dat1[[i]] <- random_points(get(paste(study_area, "poly_coug", sep = "_")), n = nrow(dat)*10, presence = tracks) %>%
    mutate(ID = i, studyArea = study_area)
}

# extract human footprint
coug_rsf_dat2 <- coug_rsf_dat1 %>% bind_rows() %>%
  mutate(human_footprint = raster::extract(hfi, .[,c("x_","y_")]),
         case_binary = ifelse(case_ == T, 1, 0), ID = factor(paste(ID))) %>%
  # add sample size to use as filter
  group_by(ID) %>% mutate(n = n()/11) %>% ungroup() %>%
  # filter for animals with at least 100 locs
  filter(n >= 100) %>%
  # add large weight for available points, as recommended for RSFs by Muff et al (2020) JAE
  mutate(weight = ifelse(case_ == TRUE, 1, 1000))

# fit RSF model with glmmTMB 
glmmTMB_coug <- glmmTMB(case_ ~ human_footprint + I(human_footprint^2) + 
                         # random intercepts and slopes
                         (1 | ID) + (0 + human_footprint | ID) + (0 + I(human_footprint^2) | ID),
                       family=binomial(), doFit=T,
                       data = coug_rsf_dat2, weights = weight) 
#saveRDS(glmmTMB_coug, "Analysis/saved models/glmmTMB_coug_20230215.rds")
#glmmTMB_coug <- readRDS("Analysis/saved models/glmmTMB_coug_20230215.rds")
summary(glmmTMB_coug)

# create new data to predict log-RSS
coug_x1 <- data.frame(human_footprint = seq(-0.5, 0.5, by = 0.01)) %>% mutate(ID = NA, weight = 1)
coug_x2 <- coug_x1 %>% filter(human_footprint == min(human_footprint))

# predict population-level log-RSS by excluding random effects; see functions script for details
coug_logRSS <- calc_logRSS_glmmTMB(x1 = coug_x1, x2 = coug_x2, model = glmmTMB_coug, ci_level = 0.95, model_type = "RSF")
#saveRDS(coug_logRSS, "Analysis/saved models/coug_logRSS_20230215.rds")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 3.3 export figure of wolf and cougar RSS ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

cougar_orange <- image_read("Analysis/cougar_orange1.png")
wolf_grey <- image_read("Analysis/wolf_grey.png")

dat_wolf_coug <- bind_rows(wolf_logRSS %>% mutate(Species = "Wolf"),  
                           coug_logRSS %>% mutate(Species = "Cougar"))

wolf_coug_rsf <- ggplot() +
  geom_line(dat_wolf_coug, mapping = aes(human_footprint, logRSS, 
                                         colour = Species), size = 0.75) +
  geom_ribbon(data = dat_wolf_coug,
              mapping = aes(x = human_footprint, ymin = lwr, ymax = upr, fill = Species), alpha = 0.3) +
  labs(y = "Log relative selection strength") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c("chocolate1","grey40")) +
  scale_fill_manual(values = c("chocolate1","grey40")) +
  labs(x = "Human footprint index") +
  theme(legend.position = c(0.25,0.25)) +
  guides(fill=guide_legend(title=element_blank()), colour=guide_legend(title=element_blank())) 
wolf_coug_rsf

#pdf(file = "Analysis/wolf_coug_rsf_Fig1_logRSS_20230216.pdf", width = 2.5, height = 2.4) 
wolf_coug_rsf +
  draw_image(cougar_orange, scale = 1, x = -0.75, y = -2) +
  draw_image(cougar_orange, scale = 1, x = -0.75, y = -2) +
  draw_image(wolf_grey, scale = 1.2, x = -0.12, y = 0.33)
dev.off()

#jpeg("Analysis/wolf_coug_rsf_Fig1_logRSS_20230216.jpg",width = 2.5, height = 2.4, units = "in", res = 1000)
wolf_coug_rsf +
  draw_image(cougar_orange, scale = 1, x = -0.75, y = -2) +
  draw_image(cougar_orange, scale = 1, x = -0.75, y = -2) +
  draw_image(wolf_grey, scale = 1.2, x = -0.12, y = 0.33)
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 4.0 Format mesopredator data for SSF ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# load and format coyote and bobcat GPS locations
meso <- readRDS("final_project/data_raw/Prugh_provided_files/WPPP_mesopredator_GPS_locations.rds") %>%
  rename(date_time = Acquisition.Time, Lat = GPS.Latitude, Long = GPS.Longitude, ID = AnimalID) %>%
  group_by(ID) %>%
  arrange(date_time) %>%
  ungroup() %>%
  mutate(Sp = ifelse(Species == "BOB", "Bobcat", "Coyote")) %>%
  # excluding MVBOB71M after dispersal date
  filter(ID != "MVBOB71M" | date_time < as.POSIXct("2019-09-24 00:00:00")) 

# plot a random subset of the locations (to reduce processing time)
meso_for_plot <- meso %>% mutate(rand = round(runif(nrow(meso), 0, 5)), mock = 1) %>% filter(rand == 1) 
leaflet() %>% 
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
  addCircles(data = meso_for_plot , 
             color = ~ colorFactor(cols25(), domain = factor(ID))(factor(ID)),
             label = paste(meso_for_plot$ID,
                           meso_for_plot$date_time), opacity = 1) 

# select variables, nest by individual animal, and make 'amt' track
meso1 <- meso %>% 
  data.frame() %>%
  dplyr::select(ID, StudyArea, Species, Sex, date_time, Lat, Long) %>% 
  nest(.by=c(ID, StudyArea, Species)) %>%
  mutate(trk = lapply(data, function(d){
    make_track(d, .x = Long, .y = Lat, .t = date_time, crs = 4326)}))

# summarise sampling rate
mesoSummary <- meso1 %>% mutate(sr = lapply(trk, summarize_sampling_rate, time_unit = "hour")) %>%
  dplyr::select(ID, sr) %>% unnest %>% 
  left_join(distinct(data.frame(meso)[,c("ID", "Sex", "StudyArea", "Species")])) %>%
  arrange(Species, median)
print(mesoSummary, n = 100)

# separate by species
coyote <- meso1 %>% filter(grepl("COY", ID))
bobcat <- meso1 %>% filter(grepl("BOB", ID))

# resample and filter tracks for temporal consistency 
coyote1 <- coyote[-8,] %>% # omitting coyote in row 8 because tripping up function due to too few consecutive data points 
  mutate(stp = lapply(trk, function(x) {
  x %>% track_resample(rate = hours(4), tolerance = minutes(10)) %>%
    filter_min_n_burst(min_n = 3) %>% 
    steps_by_burst() %>%  # calc step lengths
    random_steps(sl_distr = fit_distr(.$sl_, "gamma"), n_control = 10) %>% # Distribute 10 available points for each used location, based on the gamma distribution
    mutate(log_sl_ = log(sl_ + .1), cos_ta = cos(ta_))})) %>%
  dplyr::select(-data, -trk) %>%
  unnest(cols = stp) %>%
  left_join(distinct(data.frame(meso)[,c("ID", "Sex", "StudyArea", "Species")])) %>%
  mutate(case_binary = ifelse(case_ == TRUE, 1, 0))

# resample and filter tracks for temporal consistency 
bobcat1 <- bobcat[bobcat$ID != "MVBOB87M",] %>% # omitting MVBOB87M because tripping up function due to too few consecutive data points 
  # drop individuals with fewer than 50 points  
  filter(map_int(data, nrow) > 50) %>%
  mutate(stp = lapply(trk, function(x) {
    x %>% track_resample(rate = hours(4), tolerance = minutes(10)) %>%
      filter_min_n_burst(min_n = 3) %>%
      steps_by_burst()})) %>%
  mutate(stp = lapply(stp, function(x) {
    x %>% random_steps(n_control = 10)})) %>%  # Distribute 10 available points for each used location
  dplyr::select(-data, -trk) %>%
  unnest(cols = stp) %>%
  left_join(distinct(data.frame(meso)[,c("ID", "Sex", "StudyArea", "Species")])) %>%
  mutate(case_binary = ifelse(case_ == TRUE, 1, 0))

#saveRDS(coyote1, "Analysis/coyote_extracted_20220917.rds")
#saveRDS(bobcat1, "Analysis/bobcat_extracted_20220917.rds")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 4.1 Extract predictors ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Read processed meso data 
coyote_extracted <- readRDS("Analysis/coyote_extracted_20220917.rds") %>%
  dplyr::select(-log_sl_, -sl_) %>%
  mutate(# step lengths from amt are in degrees; overwrite to meters
         sl_ = distGeo(.[,c("x1_","y1_")],.[,c("x2_","y2_")]),
         log_sl = log(sl_),
         # instead of using calendar years, we start the year on 1 Dec of the prior year
         analysis_year = ifelse(month(t2_) == 12, year(t2_) + 1, year(t2_)),
         season = ifelse(month(t2_) %in% 4:11, "summer", "winter")) 
bobcat_extracted <- readRDS("Analysis/bobcat_extracted_20220917.rds") %>%
  dplyr::select(-sl_) %>%
  mutate(# step lengths from amt are in degrees; overwrite to meters
         sl_ = distGeo(.[,c("x1_","y1_")],.[,c("x2_","y2_")]),
         log_sl = log(sl_),
         # instead of using calendar years, we start the year on 1 Dec of the prior year
         analysis_year = ifelse(month(t2_) == 12, year(t2_) + 1, year(t2_)),
         season = ifelse(month(t2_) %in% 4:11, "summer", "winter")) 

# COYOTES: extract time-varying predictors
coyote_extracted1 <- data.frame()
for(i in 2018:2020) { # for each year
  print(i)
  dat <- data.frame()
  for(j in unique(coyote_extracted$season)) { # and for each season within each year
    print(j)
    dat0 <- coyote_extracted %>% filter(analysis_year == i, season == j) %>% 
      # wolf UD
      mutate(wolfUD = raster::extract(UD_stack_scale[[c(grep(paste(i, j, sep = "_"), names(UD_stack_scale)))]], .[,c("x2_","y2_")])) %>%
      # cougar UD
      mutate(cougarUD = raster::extract(cougar_UD_stack_scale[[c(grep(paste(i, j, sep = "_"), names(cougar_UD_stack_scale)))]], .[,c("x2_","y2_")])) %>%
      # % forest cover
      mutate(forest_pc = raster::extract(forest_pc[[c(grep(i, names(forest_pc)))]], .[,c("x2_","y2_")])) %>%
      # human footprint index
      mutate(human_footprint = raster::extract(human_footprint_scaled[[c(grep(i, names(human_footprint_scaled)))]], .[,c("x2_","y2_")]))
    dat <- rbind(dat, dat0)
    }
  coyote_extracted1 <- rbind(coyote_extracted1, dat)
}
# save fully processed data for SSF
#saveRDS(coyote_extracted1, "Analysis/coyote_extracted1_20230215.rds")

coyote_extracted2 <- coyote_extracted1 %>% 
  # create unique step ID for each animal
  mutate(step_id_ = paste(ID, step_id_, sep = "_"),
         log_sl_ = log_sl,
         human_footprint = human_footprint,
         forest_pc = forest_pc - 0.5, #subtracting 0.5 to get forest cover on same -0.5 to 0.5 scale as the other variables  
         wolfUD = wolfUD,
         cougarUD = cougarUD) %>%
  group_by(ID) %>% 
  # calculate sample size for each animal (dividing by 11 because of available points)
  mutate(n = n()/11) %>%
  ungroup()
# save file
# saveRDS(coyote_extracted2, "Analysis/coyote_extracted2_20230215.rds")


# BOBCATS: extract time-varying predictors
bobcat_extracted1 <- data.frame()
for(i in 2018:2020) { # for each year
  print(i)
  dat <- data.frame()
  for(j in unique(bobcat_extracted$season)) { # and for each season within each year
    print(j)
    dat0 <- bobcat_extracted %>% filter(analysis_year == i, season == j) %>% 
      # wolf UD
      mutate(wolfUD = raster::extract(UD_stack_scale[[c(grep(paste(i, j, sep = "_"), names(UD_stack_scale)))]], .[,c("x2_","y2_")])) %>%
      # cougar UD
      mutate(cougarUD = raster::extract(cougar_UD_stack_scale[[c(grep(paste(i, j, sep = "_"), names(cougar_UD_stack_scale)))]], .[,c("x2_","y2_")])) %>%
      # % forest cover
      mutate(forest_pc = raster::extract(forest_pc[[c(grep(i, names(forest_pc)))]], .[,c("x2_","y2_")])) %>%
      # human footprint index
      mutate(human_footprint = raster::extract(human_footprint_scaled[[c(grep(i, names(human_footprint_scaled)))]], .[,c("x2_","y2_")]))
    dat <- rbind(dat, dat0)
  }
  bobcat_extracted1 <- rbind(bobcat_extracted1, dat)
}
# save fully processed data for SSF
#saveRDS(bobcat_extracted1, "bobcat_extracted1_20230215.rds")

# read fully pre-processed data for SSF
bobcat_extracted2 <- bobcat_extracted1 %>% 
  # create unique step ID for each animal
  mutate(step_id_ = paste(ID, step_id_, sep = "_"),
         log_sl_ = log_sl,
         human_footprint = human_footprint,
         forest_pc = forest_pc - 0.5, #subtracting 0.5 to get forest cover on same -0.5 to 0.5 scale as the other variables
         wolfUD = wolfUD,
         cougarUD = cougarUD) %>%
  group_by(ID) %>% 
  # calculate sample size for each animal (dividing by 11 because of available points)
  mutate(n = n()/11) %>%
  ungroup()
# save fully processed data for SSF
#saveRDS(bobcat_extracted2, "Analysis/bobcat_extracted2_20230215.rds")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 5.0 Coyote Step Selection Function ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# read fully processed data for SSFs
coy_ssf_dat <- readRDS("Analysis/coyote_extracted2_20230215.rds")  %>% 
  # select coyotes with more than 100 fixes
  filter(n <= 100)

# plots showing bi-variate density of data, with respect to humans and large carnivores, showing that it is rare to have both high human and high large carnivore
# function to create bivariate density plots
hex_plot_fun <- function(dat, yvar, ylab, title_) {
  ggplot(dat , aes(human_footprint, get(yvar))) +
    geom_bin_2d(aes(fill = stat(log(density)))) +
    scale_fill_continuous(type = "viridis", limits=c(-14,-1)) +
    theme_bw() +
    labs(x = "human footprint", y = ylab, title = title_)+
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
    stat_density_2d(colour = "red", breaks = c(0.05), n = 15, size = 1)+
    geom_vline(xintercept = 0, colour = "gray70", size = 1) +
    geom_hline(yintercept = 0, colour = "gray70", size = 1)
}
hex1 <- hex_plot_fun(dat = coy_ssf_dat, yvar = "wolfUD", ylab = "Wolf UD", title_ = "Human and wolf"); hex1
hex2 <- hex_plot_fun(dat = coy_ssf_dat, yvar = "cougarUD", ylab = "Cougar UD", title_ = "Human and cougar"); hex2

# fit SSF with glmmTMB following Muff et al (2020) JAE
nt <- min(parallel::detectCores()) - 2
glmmTMB_coy <-  glmmTMB(case_ ~ -1 + 
                          wolfUD + (0 + wolfUD | ID) +
                          cougarUD + (0 + cougarUD | ID) + 
                          human_footprint + I(human_footprint^2) + (0 + human_footprint | ID) + (0 + I(human_footprint^2) | ID) +
                          human_footprint*wolfUD + (0 + human_footprint:wolfUD | ID) +
                          human_footprint*cougarUD + (0 + human_footprint:cougarUD | ID) +
                          forest_pc + I(forest_pc^2) + (0 + forest_pc | ID) + (0 + I(forest_pc^2) | ID) +
                          # control for step length
                          log_sl_ + (0 + log_sl_ | ID) + 
                          # strata
                          (1|step_id_),
                          family=poisson, doFit=T,
                          data = coy_ssf_dat, 
                          # Tell glmmTMB not to change the last standard deviation, all other values are estimated freely
                          map = list(theta = factor(c(1:9, NA))),
                          # Set the value of the standard deviation of the strata (the last ranef) to large constant value. See Muff et al (2020) JAE
                          start = list(theta = c(rep(0, times = 9),log(1000))),
                          control = glmmTMBControl(parallel = nt)) 
#saveRDS(glmmTMB_coy, "Analysis/glmmTMB_coy_20230215_minusToPlusScale.rds")
#glmmTMB_coy <- readRDS("Analysis/glmmTMB_coy_20230215_minusToPlusScale.rds")
summary(glmmTMB_coy)

# Graphs of coyote selection
# avg-effects plot, following the method of Avgar et al 2017 Ecol and Evol
# First step is to predict from model in order to produce data needed to make avg effect plots.
# Only predicting to random subset of points from each strata to keep prediction time manageable 
# as this prediction step takes a *very* long time.
coy_ssf_dat1 <-  coy_ssf_dat %>%
  group_by(step_id_) %>%
  mutate(samp = sample(n())) %>%
  filter(samp <= 5, case_ == F) %>% mutate(ID = NA) 
coy_pred <- predict(glmmTMB_coy, coy_ssf_dat1, re.form = NA, se.fit = T)
coy_ssf_dat1$fit <- coy_pred$fit
coy_ssf_dat1$se <- coy_pred$se
coy_ssf_dat1 <- coy_ssf_dat1 %>% ungroup()
#saveRDS(coy_ssf_dat1, "Analysis/coy_ssf_dat1_20230215_scaled.rds")
#coy_ssf_dat1 <- readRDS("Analysis/coy_ssf_dat1_20230215_scaled.rds")

# average effect plots
avg_eff_coy_human_wolf <- avg_eff_plot(fittedResponse = coy_ssf_dat1 %>% slice_sample(n = 100000), 
                                       xvar = "human_footprint", predator = "wolfUD", predatorThreshold = -0.49,
                                       predictorSpName = "wolf", nsim = 1000, showPeakValue = T)
avg_eff_coy_human_cougar <- avg_eff_plot(fittedResponse = coy_ssf_dat1 %>% slice_sample(n = 100000), 
                                         xvar = "human_footprint", predator = "cougarUD", predatorThreshold = -0.49, 
                                         predictorSpName = "cougar", nsim = 1000, showPeakValue = T)
# report peaks
avg_eff_coy_human_wolf$peaks
avg_eff_coy_human_cougar$peaks


# Relative Selection Strength (RSS); see details in functions script
# wolves
coy_wolf_RSS_dat <- find_observed_range_singleVar(fitted_data = coy_ssf_dat1, xvar = "wolfUD", xvar_point_of_comparison = min, model = glmmTMB_coy,  predatorSp = "wolf") 
coy_wolf_RSS <- calc_logRSS_glmmTMB(x1 = coy_wolf_RSS_dat$x1, x2 = coy_wolf_RSS_dat$x2, model = glmmTMB_coy, ci_level = 0.95, model_type = "SSF")
# cougars
coy_cougar_RSS_dat <- find_observed_range_singleVar(fitted_data = coy_ssf_dat1, xvar = "cougarUD", xvar_point_of_comparison = min, model = glmmTMB_coy, predatorSp = "cougar") 
coy_cougar_RSS <- calc_logRSS_glmmTMB(x1 = coy_cougar_RSS_dat$x1, x2 = coy_cougar_RSS_dat$x2, model = glmmTMB_coy, ci_level = 0.95, model_type = "SSF")

# join large predators together and plot
coy_largePreds <- bind_rows(
  coy_wolf_RSS %>% rename(largePredatorUD = wolfUD),  
  coy_cougar_RSS %>% rename(largePredatorUD = cougarUD)) %>%
  mutate(predator = factor(predator, levels = c("wolf", "cougar")))

coy_vs_largePred_plot <- ggplot(coy_largePreds, aes(largePredatorUD, exp(logRSS), colour = predator, fill = predator)) +
  theme_minimal() +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "Large predator utilization", y = "Relative selection")+ 
  scale_colour_manual(values = c("grey40","chocolate1")) +
  scale_fill_manual(values = c("grey40","chocolate1")) +
  guides(fill=guide_legend(title=element_blank()), colour=guide_legend(title=element_blank()))  +
  theme(legend.position = "bottom", text=element_text(size=9),legend.title=element_text(size=7)) 
coy_vs_largePred_plot


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 6.0 Bobcat Step Selection Function ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# read fully pre-processed data for SSF 
bob_ssf_dat <- readRDS("Analysis/bobcat_extracted2_20230215.rds") %>%
  # select bobcats with at least 100 fixes
  filter(n >= 100)

# plots to show distribution of bi-variate data
hex3 <- hex_plot_fun(dat = bob_ssf_dat, yvar = "wolfUD", ylab = "Wolf UD", title_ = "Human and wolf"); hex3
hex4 <- hex_plot_fun(dat = bob_ssf_dat, yvar = "cougarUD", ylab = "Cougar UD", title_ = "Human and cougar"); hex4

# load silhouettes
thumbnail_coyote <- image_read("Analysis/coyote1.png")
thumbnail_bobcat <- image_read("Analysis/lynx.png")

#jpeg("Analysis/hexPlot_20230215.jpg",width=9,height=6, units = "in", res = 1000)
plot_grid(NULL, hex1, hex2, NULL, hex3, hex4, rel_widths = c(0.5,1,1.4,0.5,1,1.4), 
          label_fontface = "plain", label_size = 10, labels = c("", "A", "B", "", "C", "D")) +
  draw_image(thumbnail_coyote, scale = 0.17, x = -.42, y = 0.25) +
  draw_image(thumbnail_bobcat, scale = 0.15, x = -.42, y = -0.2)
dev.off()  

# Fit SSF with glmmTMB 
nt <- min(parallel::detectCores()) - 2
glmmTMB_bob <- glmmTMB(case_binary ~ -1 + 
                         wolfUD + (0 + wolfUD | ID) +
                         cougarUD + (0 + cougarUD | ID) +
                         human_footprint + I(human_footprint^2) + (0 + human_footprint | ID) +
                         human_footprint*wolfUD + (0 + human_footprint:wolfUD | ID) + 
                         human_footprint*cougarUD + (0 + human_footprint:cougarUD | ID) + 
                         forest_pc + (0 + forest_pc | ID) +
                         # control for step length
                         log_sl_ + (0 + log_sl_ | ID)  + 
                         # strata
                         (1|step_id_),
                       family=poisson, doFit=T,
                       data = bob_ssf_dat, 
                       # Tell glmmTMB not to change the last standard deviation, all other values are freely estimated 
                       map = list(theta = factor(c(1:7, NA))),
                       # Set the value of the standard deviation of the strata (the last ranef) to large value 
                       start = list(theta = c(rep(0,times=7), log(1000))),
                       control = glmmTMBControl(parallel = nt)) 
summary(glmmTMB_bob)
#saveRDS(glmmTMB_bob, "Analysis/glmmTMB_bob_20230215_minusToPlusScale.rds")
#glmmTMB_bob <- readRDS("Analysis/glmmTMB_bob_20230215_minusToPlusScale.rds")

# avg-effects plot, following the method of Avgar et al 2017 Ecol and Evol
# First step is to predict from model in order to produce data needed to make avg effect plots
# only predicting to random subset of points from each strata to keep prediction time manageable 
# as this prediction step takes a *very* long time.
bob_ssf_dat1 <-  bob_ssf_dat %>%
  group_by(step_id_) %>%
  mutate(samp = sample(n())) %>%
  filter(samp <= 5, case_ == F) %>% ungroup()
bob_pred <- predict(glmmTMB_bob, bob_ssf_dat1, re.form = NA, se.fit = T)
bob_ssf_dat1$fit <- bob_pred$fit
bob_ssf_dat1$se <- bob_pred$se
#saveRDS(bob_ssf_dat1, "Analysis/bob_ssf_dat1_2023015.rds")
# bob_ssf_dat1 <- readRDS("Analysis/bob_ssf_dat1_2023015.rds")

# Graphs of bobcat selection
# avg effects plots; see functions script for more details
avg_eff_bob_human_wolf <- avg_eff_plot(fittedResponse = bob_ssf_dat1 %>% slice_sample(n = 100000),  xvar = "human_footprint", predator = "wolfUD", 
                                       predatorThreshold = -0.49, predictorSpName = "wolf", nsim = 1000, showPeakValue = T) 
avg_eff_bob_human_cougar <- avg_eff_plot(fittedResponse = bob_ssf_dat1 %>% slice_sample(n = 100000),  xvar = "human_footprint", predator = "cougarUD", 
                                         predatorThreshold = -0.49, predictorSpName = "cougar", nsim = 1000, showPeakValue = F) 
# because the curve is flat for predator-present in the range of HFI -0.5 to ~0.1, we need to define the point 
# at which selection starts to decline. Thus, getting the value at which selection starts declining.
maxVal_bob_coug <- avg_eff_bob_human_cougar$plot$data %>% filter(predatorPresAbs == "absent" | xvar > 0) %>% 
  filter(est == max(est)) %>% slice(rep(1:n(), times = 2)) 
maxVal_bob_coug[c(2,4),]$est <- -Inf 
# update graph
avg_eff_bob_human_cougar$plot <- avg_eff_bob_human_cougar$plot + 
  geom_line(data = maxVal_bob_coug, aes(colour = predatorPresAbs), linetype = "dashed", size = 0.5, alpha = 0.4)

# wolves
bob_wolf_RSS_dat <- find_observed_range_singleVar(fitted_data = bob_ssf_dat1, xvar = "wolfUD", 
                                                  xvar_point_of_comparison = min, model = glmmTMB_bob, predatorSp = "wolf") 
bob_wolf_RSS <- calc_logRSS_glmmTMB(x1 = bob_wolf_RSS_dat$x1, x2 = bob_wolf_RSS_dat$x2, model = glmmTMB_bob, ci_level = 0.95, model_type = "SSF")

# cougars
bob_cougar_RSS_dat <- find_observed_range_singleVar(fitted_data = bob_ssf_dat1, xvar = "cougarUD", 
                                                    xvar_point_of_comparison = min, model = glmmTMB_bob, predatorSp = "cougar") 
bob_cougar_RSS <- calc_logRSS_glmmTMB(x1 = bob_cougar_RSS_dat$x1, x2 = bob_cougar_RSS_dat$x2, model = glmmTMB_bob, ci_level = 0.95, model_type = "SSF")

# plot large predators together
bob_largePreds <- bind_rows(
  bob_wolf_RSS %>% rename(largePredatorUD = wolfUD),  
  bob_cougar_RSS %>% rename(largePredatorUD = cougarUD)) %>%
  mutate(predator = factor(predator, levels = c("wolf", "cougar")))

bob_vs_largePred_plot <- ggplot(bob_largePreds, aes(largePredatorUD, exp(logRSS), colour = predator, fill = predator)) +
  theme_minimal() +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.3, colour = NA) +
  coord_cartesian(ylim=c(0, 2)) +
  geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(x = "Large predator utilization", y = "Relative selection")+ 
  scale_colour_manual(values = c("grey40","chocolate1")) +
  scale_fill_manual(values = c("grey40","chocolate1")) +
  guides(fill=guide_legend(title=element_blank()), colour=guide_legend(title=element_blank()))  +
  theme(legend.position = "bottom", text=element_text(size=9),legend.title=element_text(size=7)) 
bob_vs_largePred_plot


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 7.0 export multi-panel meopredator plot ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# load silhouettes
thumbnail_cougar <- image_read("Analysis/cougar.png")
thumbnail_wolf <- image_read("Analysis/wolf.png")
thumbnail_coyote <- image_read("Analysis/coyote1.png")
thumbnail_bobcat <- image_read("Analysis/lynx.png")

meso_multi <- plot_grid(
  nrow = 2, align = "hv", axis = "lr", 
  labels = c("", "A", "C", "E", "", "B", "D", "F"),
  label_size = 8, label_fontface = "plain", label_x = 0.0, label_y = c(1.05,1,1,1,1.05,1.05,1.05,1.05),
  rel_widths = c(0.45,1,1,1,0.65,1,1,1), 
  rel_heights = c(0.7,1),
  
  # padding on left for silhouette
  NULL,
  
  # top left
  coy_vs_largePred_plot + theme(legend.position = "none", text=element_text(size=8)) + 
    lims(y = c(0,2.5))  +
    labs(x = element_blank()) +
    theme(panel.grid.major = element_line(size = 0.3), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 6), axis.title = element_text(size = 8)), 
  
  # top middle
  avg_eff_coy_human_wolf$plot + 
    theme(legend.position = "none", text=element_text(size=8), axis.title.x = element_blank()) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim=c(-1.5, 0.6))  +
    theme(panel.grid.major = element_line(size = 0.3), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 6), axis.title = element_text(size = 8)),
  
  # top right
  avg_eff_coy_human_cougar$plot + theme(legend.position = "none", text=element_text(size=8), axis.title.x = element_blank()) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim=c(-1.5, 0.6))  +
    theme(panel.grid.major = element_line(size = 0.3), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 6), axis.title = element_text(size = 8)), 

  # padding on left for silhouette
  NULL,
  
  # bottom left
  bob_vs_largePred_plot + theme(legend.position = "bottom", text=element_text(size=8),legend.key.size = unit(0.3, 'cm')) +
    coord_cartesian(ylim=c(0, 2.5))  +
    theme(panel.grid.major = element_line(size = 0.3), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 6), axis.title = element_text(size = 8), 
          legend.text = element_text(size = 7), legend.title = element_text(size = 7)), 

  # bottom middle
  avg_eff_bob_human_wolf$plot + 
    theme(legend.position = "bottom", text=element_text(size=7), legend.title = element_text(size=7), legend.key.size = unit(0.3, 'cm')) +
    labs(x = "Human footprint index") +
    coord_cartesian(ylim=c(-0.6, 0.2))  +
    theme(panel.grid.major = element_line(size = 0.3), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 6), axis.title = element_text(size = 8), 
          legend.text = element_text(size = 7), legend.title = element_text(size = 7)),

  # bottom right
  avg_eff_bob_human_cougar$plot + 
    theme(legend.position = "bottom", text=element_text(size=7), legend.title = element_text(size=7), legend.key.size = unit(0.3, 'cm')) +
    labs(x = "Human footprint index")  +
    coord_cartesian(ylim=c(-0.6, 0.2))  +
    theme(panel.grid.major = element_line(size = 0.3), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 6), axis.title = element_text(size = 8), 
          legend.text = element_text(size = 7), legend.title = element_text(size = 7))
  
  ) +
    draw_image(thumbnail_coyote, scale = 0.12, x = -.43, y = 0.29) +
    draw_image(thumbnail_bobcat, scale = 0.12, x = -.43, y = -0.08)

#jpeg("Analysis/meso_relative_selection_20230215.jpg",width=6.3,height=3.35, units = "in", res = 1000)
meso_multi
dev.off()

#pdf(file = "Analysis/meso_relative_selection_20230215.pdf", width = 6.3, height = 3.25) 
meso_multi
dev.off()
