# Functions for analysis of lethal human shields 
# Code written by C. Cunningham 2023-02-17
# Accompanies the paper: LR Prugh, CX Cunningham, RM Windell, BN Kertson, 
# R Ganz, SL Walker, AJ Wirsing. The paradox of the lethal human shield outside protected areas. 

# These functions accompany the main script and are used to calculate: 1) utilization distributions of 
# wolves and cougars, 2) relative selection strength, 3) average effects plots, 4) and help with plotting

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Wolf UDs ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

#create UD for each year*season for each relevant wolf pack
# wolfPackInfo is a table that contains the relevant information about pack sizes, and what year to get data from in the 
# case that a pack was known to exist but missing a collar.
# min_locs sets the threshold for minimum number of locations needed to construct a seasonal UD from a wolf's locations
wolf_UD_function <- function(wolfPackInfo, wolf_gps_data, min_locs = 50){
  UD_list <- list()
  KDE_list <- list()
  isopleth_list <- list()
  UD_summary <- data.frame()
  
  for(i in 1:nrow(wolfPackInfo)) {
    print(i)
    season_ <- wolfPackInfo[i,"season"]
    # select year of data to be used for each UD; when no data available in given year, choosing nearby years specified in wolfPackInfo 
    yr <- wolfPackInfo[i,"Year"]
    data_for_UD <- wolfPackInfo[i,"Year_for_UD"]
    packName <- wolfPackInfo[i,"pack"]
    
    # select seasonal data from specific time window
    dat <- wolf_gps_data %>% 
      filter(pack == packName) %>% # filter for pack of interest
      filter(analysisYear == data_for_UD & season == season_) # select data within pre-specified year and season
    
    # keep individuals with more than min_locs data points
    indiv <- dat %>% group_by(ID) %>% count() %>% filter(n >= min_locs)
    
    # calculate UDs, nesting by individual
    UDs <- dat %>% 
      filter(ID %in% indiv$ID) %>%
      nest(track_data = -ID) %>%
      mutate(hr_kde_ = map(track_data, hr_kde, trast = raster_template)) %>%
      mutate(ud = map(hr_kde_, hr_ud)) %>%
      mutate(isopleth = map(ud, hr_isopleths, level = 0.95)) # also get isopleths, for creating a map
    
    # when multiple pack members collared, average the UDs
    n_UDs <- nrow(UDs)
    if(n_UDs > 1) {
      # print what year we are using for each pack
      print(paste(packName, ": Using", data_for_UD, "data for", yr, "UD; Averaging", n_UDs, "animals" ))
      # create blank UD into which we will average all pack member's UDs
      UD_avg <- UDs[1,]$ud[[1]]
      UD_avg[!is.na(UD_avg)] <- 0 
      for(j in 1:n_UDs) {
        UD_avg <- UD_avg + UDs[j,]$ud[[1]]}
      UD <- UD_avg / n_UDs}
    
    # otherwise, if only one member collared, use the single collared animal's UD
    if(n_UDs == 1) {
      print(paste(packName, ": Using", data_for_UD, "data for", yr, "UD; single animal" ))
      UD <- UDs$ud[[1]]}
    
    # multiply each pack's UD by the minimum known pack size, as reported by Washington Department of Fish and Wildlife
    UD <- UD * wolfPackInfo[i,"minPackSize"]
    
    # add products to lists
    UD_list[[paste(packName, yr, season_, sep = "_")]] <- UD
    KDE_list[[paste(packName, yr, season_, sep = "_")]] <- UDs$hr_kde_
    isopleth_list[[paste(packName, yr, season_, sep = "_")]] <- UDs$isopleth
    dat1 <- UDs %>% unnest(track_data) %>% group_by(ID) %>% count() %>%
      mutate(Year = yr, Season = season_, Pack = packName, numberIndividuals = n_UDs,
             DataFromCorrectYear = ifelse(yr == data_for_UD, "Yes", paste("No:", data_for_UD)))
    
    UD_summary <- rbind(UD_summary, dat1)
  }
  return(list(UD = UD_list, KDE = KDE_list, isopleths = isopleth_list, UD_summary = UD_summary))
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Cougar UDs ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# create UD for each year*season, averaged across all cougars
cougar_UD_function <- function(cougar_gps_data, min_locs = 50){
  UD_list_coug <- list()
  UD_summary_coug <- data.frame()
  for(i in 2018:2020) {
    print(i)
    
    # select year of data to be used for each UD, choosing nearby years when no data available in a given year 
    coug_for_UD <- cougar_gps_data[cougar_gps_data$analysis_year == i,]
    
    # keep individuals with more than min_locs data points
    indiv <- coug_for_UD %>% data.frame() %>% group_by(ID) %>% count() %>% filter(n >= min_locs)
    
    # separately for the two study areas, and nested by individual cougar, calc UDs
    UDs_NE <- coug_for_UD %>% 
      filter(ID %in% indiv$ID, StudyArea == "NE") %>%
      nest(track_data = -ID) %>%
      mutate(hr_kde_ = map(track_data, hr_kde, trast = raster_template1)) %>%
      mutate(ud = map(hr_kde_, hr_ud)) 
    UDs_MV <- coug_for_UD %>% 
      filter(ID %in% indiv$ID, StudyArea == "MV") %>%
      nest(track_data = -ID) %>%
      mutate(hr_kde_ = map(track_data, hr_kde, trast = raster_template1)) %>%
      mutate(ud = map(hr_kde_, hr_ud)) 
    
    # create blank UD into which we will average all cougars
    UD_avg_NE <- UDs_NE[,]$ud[[1]]
    UD_avg_NE[!is.na(UD_avg_NE)] <- 0 
    UD_avg_MV <- UDs_MV[,]$ud[[1]]
    UD_avg_MV[!is.na(UD_avg_MV)] <- 0 
    n_UDs_NE <- nrow(UDs_NE)
    n_UDs_MV <- nrow(UDs_MV)
    
    # average the UDs at each site
    for(j in 1:n_UDs_NE) {
      UD_avg_NE <- UD_avg_NE + UDs_NE[j,]$ud[[1]]}
    UD_NE <- UD_avg_NE / n_UDs_NE
    
    for(k in 1:n_UDs_MV) {
      UD_avg_MV <- UD_avg_MV + UDs_MV[k,]$ud[[1]]}
    UD_MV <- UD_avg_MV / n_UDs_MV
    
    # add products to list
    UD_list_coug[[paste("year", i, sep = "_")]] <- (UD_NE + UD_MV)/2
  }
  return(list(UD = UD_list_coug))
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# create new data frame for use in calculating log-RSS ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# create new data frame based on the observed values of the data. Output is for use in logRSS functions
# fitted data = data used to fit SSF/RSF model
# xvar = the explanatory variable of interest
# xvar_point_of_comparison = a function to determine what value we are are setting as the point of comparison, i.e. what relative selection strength is relative to
# byAnimal = whether we calculate separately by individual (TRUE) or by population (FALSE)
# predatorSp = character string specifying the species of large predator
find_observed_range_singleVar <- function(fitted_data, xvar, xvar_point_of_comparison = min, model, byAnimal = F, predatorSp) {
  
  if(byAnimal == T) {
    newdata <- data.frame()
    # find the range of xvar observed for each animal
    for(i in unique(fitted_data$ID)) {
      indiv <- fitted_data %>% filter(ID == i) %>% 
        summarise(minDat = min(get(xvar)), 
                  maxDat = quantile(get(xvar), probs = 0.9995))
      dat <- data.frame(ID = i, xvar = seq(indiv$minDat, indiv$maxDat, length = 101))
      newdata <- bind_rows(newdata, dat)
    }
    names(newdata)[2] <- xvar
  }
  
  if(byAnimal == F) {
    minMax <- fitted_data %>% summarise(minDat = min(get(xvar)), maxDat = quantile(get(xvar), probs = 0.9995))
    newdata <- data.frame(ID = NA, xvar = seq(minMax$minDat, minMax$maxDat, length = 101))                               
  }
  
  # find all variables used in the model
  # separate interactions into individual terms
  n = names(model$frame)
  nn <- strsplit(n, ":")
  
  vars <- c()
  for(i in 1:length(nn)) {
    nnn <- nn[[i]]
    for(j in 1:length(nnn)) {
      vars <- c(vars, nnn[j])
    }
  }
  
  # remove quadratic terms
  vars <- unique(vars)
  vars <- vars[c(grepl("2)", vars) == FALSE)]
  vars <- vars[c(grepl("3)", vars) == FALSE)]
  vars <- vars[!vars %in% c(xvar, "step_id_", "ID", "case_", "(weights)")]
  
  # add means for all other variables
  for(i in unique(vars)) {
    newdata[,i] <- mean(fitted_data[[i]], na.rm = T)
  }
  
  # add in step_id_ and animal ID, even though it won't be used in prediction
  newdata$step_id_ <- NA
  newdata$predator = predatorSp
  
  newdata <- newdata %>% mutate(xvar = plyr::round_any((xvar), accuracy = 0.01))
  names(newdata)[2] <- paste(xvar)
  
  # create data frame against which Relative Selection Strength is calculated for each animal (i.e. the baseline)
  if(byAnimal == T) {
    newdat_baseline <- data.frame()
    for(i in unique(newdata$ID)) {
      n <- newdata %>% filter(ID == i) %>% slice_head(n = 1)
      newdat_baseline <- bind_rows(newdat_baseline, n)
    }}
  
  if(byAnimal == F) {
    newdat_baseline <- newdata %>% filter(get(xvar) == min(get(xvar)))
  }
  
  newdat_baseline[,xvar] <- xvar_point_of_comparison(fitted_data[[xvar]])
  return(list(x1 = newdata, x2 = newdat_baseline))
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# calculate log-RSS ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Based on Brian Smith's log_rss function from the amt package. 
# See here for his code: https://bsmity13.github.io/log_rss/ 
# Below is an adaptation for use with glmmTMB

# x1 and x2 come from the function "find_observed_range_singleVar" 
# model type can either be SSF (conditional logistic regression) or RSF (logistic regression)
calc_logRSS_glmmTMB <- function(x1, x2, model, ci_level = 0.95, model_type, calc_CIs = T) {
  
  # Predict from model, with re.form = NA in order to produce population-level inference
  x1_pred <- predict(model, newdata = x1, se.fit = F, re.form = NA)
  x2_pred <- predict(model, newdata = x2, se.fit = F, re.form = NA)
  # calc logRSS by subtracting predictions (x1) from the point of comparison (x2)
  x1$logRSS <- x1_pred - x2_pred
  
  # calculate confidence intervals
  if(calc_CIs == T) {
    # this section is adapted for glmmTMB from on amt's log_rss function (which only works on glm and clogit)
    # create model matrix for x1 and x2 manually (model.matrix wasn't working with glmmtmb in same way it would with glm)
    vars <- fixef(model)
    mod.form <- paste(names(vars$cond), collapse = " + ")
  
    if(model_type == "SSF") {
      mod.form <- paste("~ -1 + ", mod.form, sep = "")
      x1_mm <- model.matrix(formula(mod.form), data = x1)
      x2_mm <- model.matrix(formula(mod.form), data = x2)}
  
    if(model_type == "RSF"){
      mod.form <- paste("~ ", mod.form, sep = "")
      x1$Intercept <- 1
      x2$Intercept <- 1
      x1_mm <- model.matrix(formula(mod.form), data = x1)
      x2_mm <- model.matrix(formula(mod.form), data = x2)
      # dropping Intercept column because it's in there twice
      x1_mm <- x1_mm[,colnames(x1_mm)!="Intercept"]
      x2_mm <- x2_mm[,colnames(x2_mm)!="Intercept"]
    }
  
    #Subtract x2 model matrix from each row of x1
    delta_mm <- sweep(data.matrix(x1_mm), 2, data.matrix(x2_mm))
    #Get model variance-covariance matrix
    m_vcov <- stats::vcov(model, full = F)
    #Get variance of log-RSS prediction
    var_pred <- diag(delta_mm %*% m_vcov$cond %*% t(delta_mm))
    #Get standard error of prediction
    logrss_se <- unname(sqrt(var_pred))
    #Get critical value
    p <- 1 - ((1 - ci_level)/2)
    zstar <- qnorm(p)
    #Compute bounds
    x1$lwr <- x1$logRSS - logrss_se * zstar  
    x1$upr <- x1$logRSS + logrss_se * zstar  
  }
  return(x1)
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# calculate average effects plots ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# this function is an adaptation of average effect plots described by Avgar et al 2017 Ecol and Evol,
# except here we fit the curves with mgcv::gam instead of ksmooth
# fittedResponse = the predicted fits of the SSF for the "available" data used to fit the model
# xvar = explanatory variable of interest
# predator = species of large predator (wolf or cougar)
# predatorThreshold = the threshold of UD above which we deem a pixel to be occupied by a predator during a given season
# k = max smoothness parameter used the generalized additive model
# predictorSpName = string name for the large predator used for plotting 
# nsim = number of times to take random draws from model fits + se
# showPeakValue = whether or not plot should identify the peak of selection (i.e., max value)

avg_eff_plot <- function(fittedResponse, xvar, predator, predatorThreshold, k = 7, predictorSpName, 
                         nsim = 20, showPeakValue) {
  
  # define predator presence/absence based on utilization distribution threshold
  fittedResponse <- fittedResponse %>% 
    # if UD is greater than predatorThreshold, then define those areas as having predators
    mutate(predatorPresAbs = ifelse(get(predator) > predatorThreshold, "present", "absent")) %>%
    # set values between 0.005 and predatorThreshold as NAs, as these represent areas where perhaps predators occurred in very low numbers
    mutate(predatorPresAbs = as.factor(ifelse(predatorPresAbs == "absent" & get(predator) > 0.005, NA, predatorPresAbs)))

  # use mgcv::gam to fit smooth function of the predicted values following the logic of the average effect plots presented by Avgar et al. 
  # however, we are additionally propagating error using a monte carlo error propagation approach We chose to do this because the data we
  # are fitting are from model predictions (+-SE) themselves, so we wish to also propagate the SSF prediction uncertainty. 
  # Fitting the gam smooth to nsim random draws from the SSF predictions and then summarising mean CIs.
  set.seed(123)
  fit_sample <- sapply(rep(1, nsim), function (x){ 
    # take random sample based on each fitted value and it's associated standard error
    rnorm(n = nrow(fittedResponse), mean = fittedResponse$fit, sd = fittedResponse$se)})
  # end up with a matrix with number of rows equal to input data, and number of columns equal to nsim
  fit_sample1 <- matrix(fit_sample, nrow = nrow(fittedResponse), ncol = nsim)

  smooth_list <- list()
  for (j in 1:nsim) {
    print(j)
    # fit smooth function
    smooth_list[[j]] <- bam(fit_sample1[,j] ~ s(get(paste(xvar)), by = predatorPresAbs, bs = "ts", k = k) + predatorPresAbs, 
                         data = fittedResponse, select = T, discrete = T, nthreads = 7) %>%
      # extract smooths
      smooth_estimates(unconditional = F, overall_uncertainty = T) %>% add_confint() %>% rename(xvar = `get(paste(xvar))`) 
  }
  
  # take average of the many runs
  avgEff <- bind_rows(smooth_list) %>% group_by(smooth, by,predatorPresAbs, xvar) %>%
    summarise(est = mean(est), upper_ci = mean(upper_ci), lower_ci = mean(lower_ci))%>%
    mutate(predatorPresAbs = forcats::fct_relevel(predatorPresAbs, "present", "absent"))

  # plot
  p <- ggplot(avgEff, aes(x = xvar, y = est, colour = predatorPresAbs, fill = predatorPresAbs, linetype = predatorPresAbs)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.4) +
    geom_line(size = 0.75) +
    theme_minimal() +
    labs(x = xvar, y = "Relative use") +
    theme(legend.position = c(0.2,0.2)) +
    scale_fill_discrete(name = predator) +
    scale_colour_discrete(name = predator) +
    scale_colour_manual(values = c("darkred","steelblue4")) +
    scale_fill_manual(values = c("darkred","steelblue4")) +
    guides(fill=guide_legend(title=predictorSpName), colour=guide_legend(title=predictorSpName), linetype=guide_legend(title=predictorSpName))

  returnlist <- list(plot = p)
  
  # if we wish to plot vertical lines to identify the maximum selection values, then find these values and add vertical line to graphs
  if(showPeakValue == T) {
    max_val <- avgEff %>% group_by(predatorPresAbs) %>% filter(est == max(est)) %>% ungroup() %>% slice(rep(1:n(), times = 2)) 
    max_val[3:4,]$est <- -Inf 
    p <- p + geom_line(data = max_val, aes(colour = predatorPresAbs), linetype = "dashed", size = 0.5, alpha = 0.4)
    returnlist <- list(plot = p, peaks = max_val)
  }
  
  print(p)
  return(returnlist)
}

