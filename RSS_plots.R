# RSS plots
coy_ssf <- readRDS("outputs/coyote_ssf.rds")
bob_ssf_sum <- readRDS("outputs/bobcat_ssf_summer.rds")
bob_ssf_win <- readRDS("outputs/bobcat_ssf_winter.rds")

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


## Bobcat Summer ----
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
  geom_line(aes(y=y), color="#de2c2c", lwd = 1.5) +
  geom_line(aes(y=ciu), color="#de2c2c", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#de2c2c", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Human Footprint", y = "RSS", title = "Bobcat") +
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
  geom_line(aes(y=y), color="#99795b", lwd = 1.5) +
  geom_line(aes(y=ciu), color="#99795b", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#99795b", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,26) +
  labs(x = "Trail Proximity", y = "RSS", title = "Bobcat") +
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
  geom_line(aes(y=y), color="#53c153", lwd = 1.5) +
  geom_line(aes(y=ciu), color="#53c153", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#53c153", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Percent Canopy Cover", y = "RSS", title = "Bobcat") +
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
  geom_line(aes(y=y), color="cornflowerblue", lwd = 1.5) +
  geom_line(aes(y=ciu), color="cornflowerblue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="cornflowerblue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Elevation", y = "RSS", title = "Bobcat") +
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

## Bobcat Winter ----
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
  geom_line(aes(y=y), color="#771818", lwd = 1.5) +
  geom_line(aes(y=ciu), color="#771818", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#771818", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Human Footprint", y = "RSS", title = "Bobcat") +
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
  geom_line(aes(y=y), color="#55371b", lwd = 1.5) +
  geom_line(aes(y=ciu), color="#55371b", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#55371b", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,.0012) +
  labs(x = "Trail Proximity", y = "RSS", title = "Bobcat") +
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
  geom_line(aes(y=y), color="#334e33", lwd = 1.5) +
  geom_line(aes(y=ciu), color="#334e33", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="#334e33", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,0.0012) +
  labs(x = "Percent Canopy Cover", y = "RSS", title = "Bobcat") +
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
  geom_line(aes(y=y), color="darkblue", lwd = 1.5) +
  geom_line(aes(y=ciu), color="darkblue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="darkblue", lwd = 1,
            alpha = 0.6, linetype = 2) +
  ylim(0,0.5) +
  labs(x = "Elevation", y = "RSS", title = "Bobcat") +
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
## Coyote ----
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
  geom_line(aes(y=y), color="firebrick", lwd = 1.5) +
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
  geom_line(aes(y=y), color="tan4", lwd = 1.5) +
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
  geom_line(aes(y=y), color="darkgreen", lwd = 1.5) +
  geom_line(aes(y=ciu), color="darkgreen", lwd = 1,
            alpha = 0.6, linetype = 2) +
  geom_line(aes(y=cil), color="darkgreen", lwd = 1,
            alpha = 0.6, linetype = 2) +
  labs(x = "Canopy Cover", y = "RSS", title = "Coyote") +
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
  geom_line(aes(y=y), color="blue", lwd = 1.5) +
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