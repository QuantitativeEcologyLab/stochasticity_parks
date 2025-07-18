library('ggplot2') #for fancy figures
library('ggpubr') #arrange multiple figures
library('lubridate') #for decimal_date function
library('mgcv') #for predict.gam
library('sf') #for importing ecozone data
library('terra') #for basemaps
#library('basemaps') #for figure basemaps
#library('tidyterra')
library('ggspatial')
library('dplyr') # for data wrangling
library('khroma') # for color palettes

# color palette for mean NDVI
NDVI_cols <- c("#e7dbce", "#e3d6c6", "#cdbf9f", "#c7b995", "#c2b58c",
               "#b5a975", "#aaa263", "#8d8e37", "#868924", "#69761f",
               "#536321", "#485921", "#3d4f21", "#2e4511", "#274009",
               "#193401", "#1d3900", "#0f2902")

#import data -----------------------------------------------
canada_albers <- 'ESRI:102001' # albers conic projection for canada

#canada shapefile
canada <- rgeoboundaries::geoboundaries("Canada") %>%
  st_transform(canada_albers)

#ecozones
ecozones <-
  st_read("Data/ecodistricts/Canada_Ecodistricts.shp", quiet = TRUE) %>%
  st_transform(canada_albers) %>%
  group_by(zone_name) %>%
  summarize(geometry = st_union(geometry))
plot(ecozones)

#protected areas
pas <-
  st_read('Data/protected-areas/BDCAPC_CPCAD_2024.shp', quiet = TRUE) %>%
  st_transform(canada_albers) %>%
  filter(BIOME != 'Marin | Marine') %>%
  st_geometry() %>%
  st_as_sf()
plot(pas)

# raster of protected areas
r_pas <- rast('Data/protected-areas/protected-areas-0-1.tif') %>%
  project(canada_albers, method = 'mode')
plot(r_pas)

# add predicted values and residuals for each observation to dataset
#' ******************UPDATE******************************************
if(file.exists(d_sum, 'Data_annotated/ndvi-data-spatial-stats.rds')){
  d_sum <- readRDS('Data_annotated/ndvi-data-spatial-stats.rds')
} else {
  m_beta <- readRDS('Models/canada-mean-ndvi-aggr-2-2025-07-18-beta.rds')
  
  d <- readRDS("Data_annotated/ndvi-data-with-fitted-and-e.rds")
  
  d_sum <- d %>%
    summarize(mean_NDVI = mean(NDVI, na.rm = TRUE),
              var_NDVI = var(NDVI, na.rm = TRUE),
              mu_hat = mean(mu_hat, na.rm = TRUE),
              mean_e = mean(e, na.rm = TRUE),
              var_e = var(e, na.rm = TRUE),
              .by = c(x, y)) %>%
    rast() %>%
    `crs<-`('EPSG:4326') %>%
    project(canada_albers) %>%
    mask(canada) %>%
    as.data.frame(xy = TRUE) %>%
    as_tibble() %>%
    mutate(prop = extract(project(rast('Data/proportion-water.tif'),
                                  canada_albers),
                          select(., x, y))[, 2],
           elev_m = extract(project(rast('Data/elevation.tif'),
                                    canada_albers),
                            select(., x, y))[, 2],
           pa = extract(project(r_pas, canada_albers),
                        select(., x, y))[, 2]) %>%
    filter(! is.na(elev_m)) %>%
    mutate(ecozone = extract(rasterize(ecozones, r_pas,
                                       field = 'zone_name'),
                             select(., x, y))[, 2],
           ecozone = case_when(ecozone == 'Boreal PLain' ~ 'Boreal Plain',
                               ecozone == 'MixedWood Plain' ~ 'Mixedwood Plain',
                               .default = ecozone))
  
  saveRDS(d_sum, 'Data_annotated/ndvi-data-spatial-stats.rds')
  
  # test
  ggplot(d_sum) +
    coord_sf(crs = canada_albers) +
    geom_raster(aes(x, y, fill = mean_NDVI)) +
    scale_fill_gradientn('Mean NDVI', colours = NDVI_cols)
}

#some data wrangling -------------------------------------------------------------
#calculate mean +  variance spatially 
# 
# VAR <- r %>%
#   group_by(x, y, park) %>%
#   summarize(mean = mean(preds),
#             var = var(res),
#             mean_res = mean(res),
#             cv = cv(preds, aszero = T),
#             cv.2 = sd(preds)/mean(preds),
#             mean.2 = mean(NDVI)) %>%
#   mutate(median_ndvi = median(NDVI), .by = ecozone) %>%
#   arrange(desc(median_ndvi)) %>%
#   mutate(ecozone = factor(ecozone, levels = unique(ecozone)))
# 
# VAR <- readRDS('Data_annotated/summarized-spatial-stats-albers.rds') %>%
#   mutate(pa = extract(r_pas, select(., x, y))[, 2],
#          ecozone = extract(rasterize(ecozones, r_spat_stats,
#                                      field = 'zone_name'),
#                            select(., x, y))[, 2])
# 
#create rasters for this object to plot easier 
# r_spat_stats <- select(VAR, - c(pa, ecozone)) %>%
#   rast() %>%
#   `crs<-`('EPSG:4326') %>%
#   project(canada_albers)
# plot(r_spat_stats)
# 
# mean.rast <- rast$mean
# terra::writeCDF(mean.rast, 'Canada/mean_predNDVI_raster.nc', overwrite = TRUE)
# 
# var.rast <- rast$var
# terra::writeCDF(var.rast, 'Canada/var_residuals_raster.nc', overwrite = TRUE)
# 
# meanres.rast <- rast$mean_res
# terra::writeCDF(meanres.rast, 'Canada/mean_residuals_raster.nc', overwrite = TRUE)
# 
# cv.rast <- rast$cv
# terra::writeCDF(cv.rast, 'Canada/coeff_variation_raster.nc', overwrite = TRUE)
# 
# meanNDVI.rast <- rast$mean.2
# terra::writeCDF(meanNDVI.rast, 'Canada/mean_rawNDVI_raster.nc', overwrite = TRUE)
# 
# RES <- r %>%
#   group_by(year, doy) %>%
#   summarize(mean = mean(res))
# 
# saveRDS(RES, "Canada/results/RES.rds")
# 
# ELEV <- r %>%
#   group_by(elevation) %>%
#   summarize(mean = mean(res))
# 
# saveRDS(ELEV, "Canada/results/ELEV.rds")
# 
# MEAN.DOY <- r %>%
#   group_by(doy, park) %>%
#   summarize(mean = mean(preds))
# 
# saveRDS(MEAN.DOY, "Canada/results/MEAN.DOY.rds")
# 
# MEAN.YEAR <- r %>%
#   group_by(year, park) %>%
#   summarize(mean = mean(preds))
# 
# saveRDS(MEAN.YEAR, "Canada/results/MEAN.YEAR.rds")
# 
# VAR.DOY <- r %>%
#   group_by(doy, park) %>%
#   summarize(var = var(res))
# 
# saveRDS(VAR.DOY, "Canada/results/VAR.DOY.rds")
# 
# VAR.YEAR <- r %>%
#   group_by(year, park) %>%
#   summarize(var = var(res))
# 
# saveRDS(VAR.YEAR, "Canada/results/VAR.YEAR.rds")
# 
# CV.YEAR <- r %>%
#   group_by(year, park) %>%
#   summarize(cv = cv(preds, aszero = T))
# 
# saveRDS(CV.YEAR, "Canada/results/CV.YEAR.rds")
# 
# CV.DOY <- r %>%
#   group_by(doy, park) %>%
#   summarize(cv = cv(preds, aszero = T))
# 
# saveRDS(CV.DOY, "Canada/results/CV.DOY.rds")
# 
# ECOZ <- r %>%
#   group_by(ecozone) %>%
#   summarize(mean = mean(NDVI),
#             var = var(preds),
#             cv = cv(preds, aszero = T))

#plot data (figure 1) -----------------------------------------------------
p_obs_mean <-
  #' *REMOVE*
  filter(d_sum, mean_NDVI < 0.75) %>% # drop two cells with extreme values
  #' **********************************************************************
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = 'transparent') + 
  geom_raster(aes(x, y, fill = mean_NDVI)) +
  scale_fill_gradientn('Mean NDVI', colours = NDVI_cols) +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, -0.5, -0.2, -0.5), "cm"), #top, right, bottom, left
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5, face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 10, barheight = 0.2,
                               title="Mean NDVI",direction = "horizontal"))

parks.plot <- p_obs_mean +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) +
  geom_sf(data = pas, fill = '#00000040', color = "black", lwd = 0.1)

ecoz.plot <-
  p_obs_mean +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) +
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.2)

fig_1 <- ggarrange(parks.plot, ecoz.plot, nrow = 1, common.legend = T,
                   labels = "AUTO")

ggsave("Figures/figure1.png", fig_1, units = "in", width = 6.46,
       height = 3, bg = "white", dpi = 600, scale = 2); beepr::beep(4)

#plot mean trends (figure 2) -----------------------------------------------------
#spatial mean trends
mean <- 
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = 'transparent') + 
  geom_raster(aes(x, y, fill = mu_hat), d_sum) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn('Mean', colours = NDVI_cols) + 
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 20, barheight = 0.5,
                               title = "Mean NDVI",
                               direction = "horizontal"))

#plot smooths

#mean trends by year
if(file.exists('Outputs/smooth-estimates-year.rds') &
   file.exists('Outputs/smooth-estimates-doy.rds')) {
  preds_year <- readRDS('Outputs/smooth-estimates-year.rds')
  preds_doy <- readRDS('Outputs/smooth-estimates-doy.rds')
} else {
  preds_year <-
    expand.grid(year = seq(1981, 2025, length.out = 500), pa = 0:1,
                doy = 0, x = 0, y = 0, prop_water = 0, elev_m = 0) %>%
    bind_cols(.,
              predict(m_beta, newdata = ., type = 'link',se.fit=TRUE,
                      terms = c('(Intercept)', 's(year)', 's(year):pa')) %>%
                as.data.frame() %>%
                transmute(mu_hat = fit,
                          lwr_95 = fit - 1.96 * se.fit,
                          upr_95 = fit + 1.96 * se.fit) %>%
                mutate(across(everything(), m_beta$family$linkinv))) %>%
    mutate(pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
             factor())
  saveRDS(preds_year, 'Outputs/smooth-estimates-year.rds')
  
  preds_doy <-
    expand.grid(doy = seq(0.5, 366.5, length.out = 500), pa = 0:1,
                year = 0, x = 0, y = 0, prop_water = 0, elev_m = 0) %>%
    bind_cols(.,
              predict(m_beta, newdata = ., type = 'link',se.fit=TRUE,
                      terms = c('(Intercept)', 's(doy)', 's(doy):pa')) %>%
                as.data.frame() %>%
                transmute(mu_hat = fit,
                          lwr_95 = fit - 1.96 * se.fit,
                          upr_95 = fit + 1.96 * se.fit) %>%
                mutate(across(everything(), m_beta$family$linkinv))) %>%
    mutate(pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
             factor())
  saveRDS(preds_doy, 'Outputs/smooth-estimates-doy.rds')
}

parkymean <- 
  ggplot(preds_year) +
  geom_ribbon(aes(year, ymin = lwr_95, ymax = upr_95, fill = pa),
              alpha = 0.3) +
  geom_line(aes(year, mu_hat, color = pa)) +
  scale_colour_manual(name = NULL, values=c("#A7C957", "darkgreen"),
                      aesthetics = c('color', 'fill')) +
  labs(x = "Year", y = "Mean NDVI") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = 'inside', legend.position.inside = c(0.8, 0.9))

#mean trends by day of year
parkdoymean <- 
  ggplot(preds_doy) +
  geom_ribbon(aes(doy, ymin = lwr_95, ymax = upr_95, fill = pa),
              alpha = 0.3) +
  geom_line(aes(doy, mu_hat, color = pa)) +
  scale_colour_manual(name = NULL, values=c("#A7C957", "darkgreen"),
                      aesthetics = c('color', 'fill')) +
  xlab("Day of year") +
  ylab("Mean NDVI") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

#boxplot to plot mean in different ecozones
boxmean <- 
  d_sum %>%
  filter(! is.na(ecozone)) %>%
  mutate(median_ndvi = median(mean_NDVI), .by = ecozone) %>%
  arrange(median_ndvi) %>%
  mutate(ecozone = factor(ecozone, levels = unique(ecozone)),
         pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
           factor()) %>%
  ggplot(aes(mu_hat, ecozone, group = paste(ecozone, pa), fill = pa)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.2, outlier.alpha = 0.1) +
  scale_fill_manual(name = "", labels = c("Outside PAs", "Within PAs"),
                    values=c("#A7C957", "darkgreen")) +
  xlab("Mean NDVI") +
  ylab("Ecozone") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.7, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
        legend.position = 'inside', legend.position.inside = c(0.1, 0.9))

fig_2 <-
  ggarrange(parkymean, parkdoymean, nrow = 2, labels = c("B", "C")) %>%
  ggarrange(mean, ., ncol = 2, labels = c("A", ""), widths = c(0.75, 0.5)) %>%
  ggarrange(boxmean, nrow = 2, labels = c("", "D"), heights = c(0.75, 0.5))

ggsave('Figures/figure2.png', fig_2, units = "in", width = 6.86, height = 8.52,
       bg = "white", dpi = 600)

#plot variance trends (figure 3) -------------------------------------------------

#spatial variance trends
#' *REMOVE*
#' ensure there is no banding
filter(d_sum,
       y > data.frame(x = -90, y = 65) %>%
         st_as_sf(coords = c('x','y')) %>%
         st_set_crs('EPSG:4326') %>%
         st_transform(canada_albers) %>%
         st_coordinates() %>%
         `[`(1, 'Y')) %>%
  mutate(var_e = if_else(var_e < 0.02, var_e, NA_real_)) %>%
  ggplot() +
  geom_raster(aes(x, y, fill = var_e)) +
  scale_fill_devon(name = 'Variance in NDVI', limits = c(0, NA),
                   reverse = TRUE, na.value = 'red') +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  # add a reference line at 68N
  geom_path(aes(X, Y, group = lat),
            expand.grid(x = -145:-63, y = c(68, 70)) %>%
              mutate(lat = y) %>%
              st_as_sf(coords = c('x','y')) %>%
              st_set_crs('EPSG:4326') %>%
              st_transform(canada_albers) %>%
              bind_cols(., st_coordinates(.)),
            linewidth = 0.5, lty = 'dashed')

filter(d_sum,
       y > data.frame(x = -90, y = 65) %>%
         st_as_sf(coords = c('x','y')) %>%
         st_set_crs('EPSG:4326') %>%
         st_transform(canada_albers) %>%
         st_coordinates() %>%
         `[`(1, 'Y')) %>%
  select(x, y, mean_NDVI, var_NDVI, mean_e, var_e) %>%
  tidyr::pivot_longer(c(mean_NDVI, var_NDVI, mean_e, var_e)) %>%
  group_by(name) %>%
  mutate(value = if_else(value > quantile(value, 0.999),
                         quantile(value, 0.999), value),
         value = (value - min(value)),
         value = value / max(value)) %>%
  ungroup() %>%
  ggplot() +
  facet_wrap(~ name) +
  geom_raster(aes(x, y, fill = value)) +
  scale_fill_smoothrainbow() +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  # add a reference line at 68N
  geom_path(aes(X, Y, group = lat),
            expand.grid(x = -145:-63, y = c(69.5, 70)) %>%
              mutate(lat = y) %>%
              st_as_sf(coords = c('x','y')) %>%
              st_set_crs('EPSG:4326') %>%
              st_transform(canada_albers) %>%
              bind_cols(., st_coordinates(.)),
            linewidth = 0.5, lty = 'dashed')

var <-
  ggplot() +
  geom_raster(aes(x, y, fill = var_e)) +
  geom_sf(data = canada, fill = NA, color = "black") +
  scale_fill_devon(name = 'Variance in NDVI', limits = c(0, NA),
                   reverse = TRUE) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(
    title.position = "top", ticks.colour = NA, barwidth = 18,
    barheight = 0.5, title = "Variance in NDVI",  direction = "horizontal"))

#plot smooths

#separate data by within/outside parks
d %>%
  group_by(pa, doy) %>%
  
  doyparkyes <- VAR.DOY[VAR.DOY$park == 1,]
doyparkno <- VAR.DOY[VAR.DOY$park == 0,]

yearparkyes <- VAR.YEAR[VAR.YEAR$park == 1,]
yearparkno <- VAR.YEAR[VAR.YEAR$park == 0,]

parkyvar <- 
  ggplot() +
  geom_point(yearparkyes, mapping = aes(year, var), size = 0.3, color = "dodgerblue3") +
  geom_point(yearparkno, mapping = aes(year, var), size = 0.3, color = "lightskyblue2") +
  geom_smooth(yearparkyes, mapping = aes(year, var, colour = "Within Parks"), span = 0.25, se = FALSE) +
  geom_smooth(yearparkno, mapping = aes(year, var, colour = "Outside Parks"),  span = 0.25, se = FALSE) +
  scale_colour_manual(name = "", values=c('lightskyblue2', 'dodgerblue3')) +
  xlab("Year") +
  ylab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.25, 0.9))

#residual trends by day of year
parkdoyvar <- 
  ggplot() +
  geom_point(doyparkyes, mapping = aes(doy, var), size = 0.3, color = "dodgerblue3") +
  geom_point(doyparkno, mapping = aes(doy, var), size = 0.3, color = "lightskyblue2") +
  geom_smooth(doyparkyes, mapping = aes(doy, var, colour = "Within Parks"), span = 0.15, se = FALSE) +
  geom_smooth(doyparkno, mapping = aes(doy, var, colour = "Outside Parks"),  span = 0.15, se = FALSE) +
  scale_colour_manual(name = "", values=c('lightskyblue2', 'dodgerblue3')) +
  xlab("Day of year") +
  ylab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

#plot mean trends by ecozone

#boxplot to plot mean in different ecozones
boxvar <- 
  ggplot(VAR, aes(x = layer, y = var, fill = park)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.2) +
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values=c('lightskyblue2', 'dodgerblue3')) +
  scale_x_discrete(labels=c("2" = "Northern Arctic", "1" = "Arctic Cordillera", "3" = "Southern Arctic", "9" = "Boreal Plain",
                            "6" = "Boreal Shield", "15" = "Hudson Plain", "14" = "Montane Cordillera", "8" = "Mixedwood Plain",
                            "13" = "Pacific Maritime", "10" = "Prairies", "11" = "Taiga Cordillera", "4" = "Taiga Plain",
                            "5" = "Taiga Shield", "7" = "Atlantic Maritime", "12" = "Boreal Cordillera")) +
  xlab("Ecozone") +
  ylab("Variance in NDVI") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1, size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.7, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.1, 0.9))

bc <- ggarrange(parkyvar, parkdoyvar, nrow = 2, labels = c("b", "c"))
abc <- ggarrange(var, bc, ncol = 2, labels = c("a", ""), widths = c(0.75, 0.5))
ggarrange(abc, boxvar, nrow = 2, labels = c("", "d"), heights = c(0.75, 0.5))

ggsave('Figures/figure3.png', units = "in", width = 6.86, height = 8.52,
       bg = "white", dpi = 600)

#plot quantiles (figure 4) ---------------------------------------------------------

#calculate quantiles
quantile(VAR$mean, probs = 0.7) #0.1127625 
quantile(VAR$var, probs = 0.3, na.rm = TRUE) #0.002437498 
quantile(VAR$cv, probs = 0.3) 

VAR$mean.quant <- with(VAR, ifelse(mean >= quantile(VAR$mean, probs = 0.7), 1, 0))
VAR$var.quant <- with(VAR, ifelse(var <= quantile(VAR$var, probs = 0.3, na.rm = TRUE), 1, 0))
VAR$cv.quant <- with(VAR, ifelse(cv <= quantile(VAR$cv, probs = 0.3), 1, 0))

VAR$mean.quant <- as.factor(VAR$mean.quant)
VAR$var.quant <- as.factor(VAR$var.quant)
VAR$cv.quant <- as.factor(VAR$cv.quant)

#save quantiles as rasters to speed up subsequent analyses
rast <- rast(VAR, type = 'xyz')

mean.quant.rast <- rast$mean.quant
terra::writeCDF(mean.quant.rast, 'Canada/mean_predNDVI_quantile_raster.nc', overwrite = TRUE)

var.quant.rast <- rast$var.quant
terra::writeCDF(var.quant.rast, 'Canada/mean_predNDVI_quantile_raster.nc', overwrite = TRUE)

cv.quant.rast <- rast$cv.quant
terra::writeCDF(cv.quant.rast, 'Canada/mean_predNDVI_quantile_raster.nc', overwrite = TRUE)

mean.quant <- 
  ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = mean.quant, alpha = mean.quant)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_manual(values = c('grey80',"darkgreen"), labels = c('bottom 70th quantile', 'top 30th quantile')) + 
  scale_alpha_discrete(guide = "none", range = c(0.6, 1), labels = c('bottom 70th quantile', 'top 30th quantile')) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.15, 0.95),
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

var.quant <- ggplot() +
  geom_raster(subset(VAR, !is.na(var.quant)), mapping = aes(x, y, fill = var.quant, alpha = var.quant)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_manual(values = c('grey80', 'dodgerblue3'), labels = c('top 70th quantile', 'bottom 30th quantile')) + 
  scale_alpha_discrete(guide = "none", range = c(0.6, 1), labels = c('bottom 70th quantile', 'top 30th quantile')) +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.15, 0.95),
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

cv.quant <- ggplot() +
  geom_raster(subset(VAR, !is.na(cv.quant)), mapping = aes(x, y, fill = cv.quant, alpha = cv.quant)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_manual(values = c("grey80", "#5b3f64"), labels = c('top 70th quantile', 'bottom 30th quantile')) + 
  scale_alpha_discrete(guide = "none", range = c(0.6, 1), labels = c('bottom 70th quantile', 'top 30th quantile')) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.15, 0.95),
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

#some data wrangling

VAR2 <- VAR
VAR2$park <- "2"

VAR.q <- rbind(VAR, VAR2)
VAR.q$park <- as.factor(VAR.q$park)

ggarrange(mean.quant, var.quant, cv.quant, p_spp_richness,
          nrow = 3, ncol = 2, labels = "AUTO")

ggsave('Figures/figure4.png', units = "in", bg = "white", width = 6.86,
       height = 8.5, dpi = 600)

#plot model residuals (appendix) ----------------------------------------------

#spatial log mean of residuals
logmeanres <- ggplot() +
  geom_raster(subset(VAR, !is.na(mean_res)), mapping = aes(x, y, fill = log(mean_res + 2*abs(min(mean_res))))) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn("log Mean Residuals", colours = colors) + #limits = c(0, 0.15)) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 18,
                               barheight = 0.5, title = "log Mean Residuals",  direction = "horizontal"))

#mean of residuals for each area - plots faster than all residuals
mrdoy <- ggplot(RES) + 
  geom_point(aes(doy, mean, alpha = 0.1)) +
  xlab("Day of year") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

mryear <- ggplot(RES) + 
  geom_point(aes(year, mean, alpha = 0.1)) +
  xlab("Year") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

mrlat <- ggplot(VAR) + 
  geom_point(aes(y, mean, alpha = 0.1)) +
  xlab("Latitude") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

mrlong <- ggplot(VAR) + 
  geom_point(aes(x, mean, alpha = 0.1)) +
  xlab("Longitude") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

mrelev <- ggplot(ELEV) + 
  geom_point(aes(elevation, mean, alpha = 0.1)) +
  xlab("Elevation (m)") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

bcf <- ggarrange(mrlat, mrlong, mrelev, nrow = 3, labels = c("b", "c", "f"))
de <- ggarrange(mrdoy, mryear, ncol = 2, labels = c("d", "e"))
ade <- ggarrange(logmeanres, de, nrow = 2, labels = c("a", ""), heights = c(1, 0.5))
ggarrange(ade, bcf, ncol = 2, widths = c(1, 0.5))

ggsave('residuals.png', units = "in", bg = "white", width = 6.86,
       height = 6.86, dpi = 600)

#scatterplot of mean vs var
ggplot() +
  geom_point

#histogram of model residuals
hist(r$res, main = "", xlab = "Model Residuals")

#coefficient of variation (appendix) ---------------------------------------------------------

#spatial variance trends
cv <- 
  ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = cv/100)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_acton(name = 'Coefficient of Variation') +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 18,
                               barheight = 0.5, title = "Coefficient of Variation",  direction = "horizontal"))

#plot smooths

#separate data by within/outside parks
doyparkyes <- CV.DOY[CV.DOY$park == 1,]
doyparkno <- CV.DOY[CV.DOY$park == 0,]

yearparkyes <- CV.YEAR[CV.YEAR$park == 1,]
yearparkno <- CV.YEAR[CV.YEAR$park == 0,]

parkycv <- 
  ggplot() +
  geom_point(yearparkyes, mapping = aes(year, cv/100), size = 0.3, color = "#5b3f64") +
  geom_point(yearparkno, mapping = aes(year, cv/100), size = 0.3, color = "#ccb3d1") +
  geom_smooth(yearparkyes, mapping = aes(year, cv/100, colour = "Within Parks"), span = 0.25, se = FALSE) +
  geom_smooth(yearparkno, mapping = aes(year, cv/100, colour = "Outside Parks"),  span = 0.25, se = FALSE) +
  scale_colour_manual(name = "", values = c("#ccb3d1", "#5b3f64")) +
  xlab("Year") +
  ylab("Coefficient of Variation") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.25, 0.2))

#residual trends by day of year
parkdoycv <- 
  ggplot() +
  geom_point(doyparkyes, mapping = aes(doy, cv/100), size = 0.3, color = "#5b3f64") +
  geom_point(doyparkno, mapping = aes(doy, cv/100), size = 0.3, color = "#ccb3d1") +
  geom_smooth(doyparkyes, mapping = aes(doy, cv/100, colour = "Within Parks"), span = 0.15, se = FALSE) +
  geom_smooth(doyparkno, mapping = aes(doy, cv/100, colour = "Outside Parks"),  span = 0.15, se = FALSE) +
  scale_colour_manual(name = "", values = c("#ccb3d1", "#5b3f64")) +
  xlab("Day of year") +
  ylab("Coefficient of Variation") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

#plot mean trends by ecozone

#boxplot to plot mean in different ecozones
boxcv <- 
  ggplot(edata, aes(x = layer, y = cv/100, fill = park)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.2) +
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values = c("#ccb3d1", "#5b3f64")) +
  scale_x_discrete(labels=c("2" = "Northern Arctic", "1" = "Arctic Cordillera", "3" = "Southern Arctic", "9" = "Boreal Plain",
                            "6" = "Boreal Shield", "15" = "Hudson Plain", "14" = "Montane Cordillera", "8" = "Mixedwood Plain",
                            "13" = "Pacific Maritime", "10" = "Prairies", "11" = "Taiga Cordillera", "4" = "Taiga Plain",
                            "5" = "Taiga Shield", "7" = "Atlantic Maritime", "12" = "Boreal Cordillera")) +
  xlab("Ecozone") +
  ylab("Coefficient of Variation") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1, size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.7, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.1, 0.9))

bc <- ggarrange(parkycv, parkdoycv, nrow = 2, labels = c("b", "c"))
abc <- ggarrange(cv, bc, ncol = 2, labels = c("a", ""), widths = c(0.75, 0.5))
ggarrange(abc, boxcv, nrow = 2, labels = c("", "d"), heights = c(0.75, 0.5))

ggsave('cv.png', units = "in", width = 6.86, height = 8.52, bg = "white",
       dpi = 600)

#correlation between mean and variance (appendix) ---------------------------------

FIG <- 
  ggplot(data=DATA, aes(x=mean, y=var)) +
  geom_point(size = 1, alpha = 0.3,stroke = 0,shape=16) +
  ylab("Variance in NDVI") +
  xlab("Mean NDVI") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_x_continuous(limits = c(-0.1,0.5), expand = c(0.01,0)) +
  scale_y_continuous(limits = c(0,0.032), expand = c(0.01,0))



#Save the figure
ggsave(FIG, width = 6.86, height = 4.5, units = "in", dpi = 600,
       bg = "white", file="Mean_Var_Correlation.png")
