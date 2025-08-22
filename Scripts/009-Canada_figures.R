library('ggplot2') #for fancy figures
library('ggpubr') #arrange multiple figures
library('lubridate') #for decimal_date function
library('mgcv') #for predict.gam
library('sf') #for importing ecozone data
library('terra') #for basemaps
library('ggspatial')
library('dplyr') # for data wrangling
library('khroma') # for color palettes
library('cowplot') # for multi-panel figures
source('Functions/scale-ndvi.R')

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
  read_sf("Data/ecodistricts/Canada_Ecodistricts.shp") %>%
  st_transform(canada_albers) %>%
  summarize(geometry = st_union(geometry), .by = zone_name)

#protected areas
pas <-
  read_sf('Data/protected-areas/BDCAPC_CPCAD_2024.shp') %>%
  filter(BIOME != 'Marin | Marine') %>%
  st_transform(canada_albers) %>%
  st_geometry() %>%
  st_as_sf()

# raster of protected areas
r_pas <- rast('Data/protected-areas/protected-areas-0-1.tif') %>%
  project(canada_albers, method = 'mode')

# import data with residuals and projected spatial summary
d <- readRDS('Data_annotated/ndvi-data-with-fitted-and-e.rds')
est_albers <- readRDS('Data_annotated/summarized-spatial-stats-albers.rds')

#plot data (figure 1) -----------------------------------------------------
p_obs_mean <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = 'transparent') + 
  geom_raster(aes(x, y, fill = unmodeled_mean), est_albers) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) +
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
  geom_sf(data = pas, fill = '#00000040', color = "black", lwd = 0.1)

ecoz.plot <-
  p_obs_mean +
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.1)

fig_1 <- ggarrange(parks.plot, ecoz.plot, nrow = 1, common.legend = T,
                   labels = "AUTO")

ggsave("Figures/figure-1.png", fig_1, units = "in", width = 6.46,
       height = 3, bg = "white", dpi = 600, scale = 2)

#plot mean trends (figure 2) -----------------------------------------------------
#spatial mean trends
p_mean_canada <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = 'transparent') + 
  geom_raster(aes(x, y, fill = mu_hat), est_albers) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_gradientn('Mean', colours = NDVI_cols) + 
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
                               barwidth = 15, barheight = 0.5,
                               title = "Mean NDVI",
                               direction = "horizontal"))

#plot smooths

#mean trends by year
if(file.exists('Outputs/smooth-estimates-year.rds') &
   file.exists('Outputs/smooth-estimates-doy.rds')) {
  preds_year <- readRDS('Outputs/smooth-estimates-year.rds')
  preds_doy <- readRDS('Outputs/smooth-estimates-doy.rds')
} else {
  m_beta <- readRDS('Models/canada-mean-ndvi-aggr-2-2025-07-31-beta.rds')
  
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
                mutate(across(everything(), m_beta$family$linkinv),
                       across(everything(), ndvi_to_11))) %>%
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
                mutate(across(everything(), m_beta$family$linkinv),
                       across(everything(), ndvi_to_11))) %>%
    mutate(pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
             factor())
  saveRDS(preds_doy, 'Outputs/smooth-estimates-doy.rds')
  rm(m_beta)
  gc()
}

# find maximum differences between PA and non-PA areas throughout the years
preds_year %>%
  select(year, pa, mu_hat) %>%
  tidyr::pivot_wider(names_from = pa, values_from = mu_hat) %>%
  mutate(diff = `Within PAs` - `Outside PAs`) %>%
  pull(diff) %>%
  abs() %>%
  max() # 0.002655057

parkymean <-
  ggplot(preds_year) +
  geom_vline(xintercept = 2014, lwd = 0.25, lty = 'dashed') +
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
        legend.position = 'inside', legend.position.inside = c(0.3, 0.2))

#mean trends by day of year
parkdoymean <-
  ggplot(preds_doy) +
  geom_ribbon(aes(doy, ymin = lwr_95, ymax = upr_95, fill = pa),
              alpha = 0.3) +
  geom_line(aes(doy, mu_hat, color = pa)) +
  scale_colour_manual(name = NULL, values = c("#A7C957", "darkgreen"),
                      aesthetics = c('color', 'fill')) +
  scale_x_continuous(NULL, breaks = c(79, 171, 263, 354),
                     labels = \(x) format(as.Date('2024-12-31') + x, '%b %d')) +
  ylab("Mean NDVI") +
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
  est_albers %>%
  filter(! is.na(ecozone)) %>%
  mutate(median_ndvi = median(mu_hat), .by = ecozone) %>%
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
  ggarrange(p_mean_canada, ., ncol = 2, labels = c("A", ""),
            widths = c(0.75, 0.5)) %>%
  ggarrange(boxmean, nrow = 2, labels = c("", "D"), heights = c(0.75, 0.5))

ggsave('Figures/figure-2.png', fig_2, units = "in", width = 6.86,
       height = 8.52, bg = "white", dpi = 600)

#plot variance trends (figure 3) -------------------------------------------------
#spatial variance trends
hist(est_albers$s2_hat)
mean(est_albers$s2_hat > 0.04) # about 0.2% of data has variance > 0.04

p_var_canada <-
  est_albers %>%
  mutate(s2_hat = if_else(s2_hat > 0.04, 0.04, s2_hat)) %>%
  ggplot() +
  geom_sf(data = canada, fill = 'grey') +
  geom_raster(aes(x, y, fill = s2_hat)) +
  geom_sf(data = canada, fill = 'transparent', color = "black", lwd = 0.2) +
  scale_fill_devon(name = 'Variance in NDVI residuals', limits = c(0, NA),
                   reverse = TRUE) +
  coord_sf(datum = canada_albers) +
  theme_void() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(
    title.position = "top", ticks.colour = NA,
    barwidth = 15, barheight = 0.5, direction = "horizontal"))

# calculate variance (i.e., mean squared residuals) over the years and doy
# while differentiating between outside and inside PAs
if(all(file.exists(c('Outputs/variance-estimates-year.rds',
                     'Outputs/variance-estimates-doy.rds')))) {
  s2_year <- readRDS('Outputs/variance-estimates-year.rds')
  s2_doy <- readRDS('Outputs/variance-estimates-doy.rds')
} else {
  s2_year <- d %>%
    summarize(s2 = var(e, na.rm = TRUE), .by = c(pa, year, x, y)) %>%
    summarize(s2_se = sqrt(var(s2, na.rm = TRUE) / n()),
              s2 = mean(s2, na.rm = TRUE),
              .by = c(pa, year)) %>%
    mutate(pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
             factor())
  
  s2_doy <- d %>%
    summarize(s2 = var(e, na.rm = TRUE), .by = c(pa, doy, x, y)) %>%
    summarize(s2_se = sqrt(var(s2, na.rm = TRUE) / n()),
              s2 = mean(s2, na.rm = TRUE),
              .by = c(pa, doy)) %>%
    mutate(pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
             factor())
  
  saveRDS(s2_year, 'Outputs/variance-estimates-year.rds')
  saveRDS(s2_doy, 'Outputs/variance-estimates-doy.rds')
}

# figure 3B
parkyvar <-
  ggplot(s2_year) +
  geom_vline(xintercept = 2014, lwd = 0.25, lty = 'dashed') +
  geom_errorbar(aes(year, ymin = s2 - s2_se, ymax = s2 + s2_se,
                    color = pa), linewidth = 0.1, width = 1) +
  geom_point(aes(year, s2, color = pa), size = 0.25) +
  scale_colour_manual(name = NULL, aesthetics = c('color', 'fill'),
                      values = c('lightskyblue2', 'dodgerblue3')) +
  xlab("Year") +
  ylab("Variance in NDVI residuals") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = 'inside', legend.position.inside = c(0.25, 0.9),
        legend.key.width = rel(0.4))

# figure 3C: residual trends by day of year
parkdoyvar <-
  ggplot(s2_doy) +
  geom_errorbar(aes(doy, ymin = s2 - s2_se, ymax = s2 + s2_se,
                    color = pa), linewidth = 0.1, width = 7.5) +
  geom_point(aes(doy, s2, color = pa), size = 0.25) +
  scale_colour_manual(name = NULL, aesthetics = c('color', 'fill'),
                      values = c('lightskyblue2', 'dodgerblue3')) +
  scale_x_continuous(NULL, breaks = c(79, 171, 263, 354),
                     labels = \(x) format(as.Date('2024-12-31') + x, '%b %d')) +
  ylab("Variance in NDVI residuals") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none",
        legend.key.width = rel(0.4))

# figure 3D: boxplots of variance in different ecozones
boxvar <-
  est_albers %>%
  filter(! is.na(ecozone)) %>%
  mutate(median_ndvi = median(s2_hat), .by = ecozone) %>%
  arrange(median_ndvi) %>%
  mutate(ecozone = factor(ecozone, levels = unique(ecozone)),
         pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
           factor()) %>%
  ggplot(aes(s2_hat, ecozone, group = paste(ecozone, pa), fill = pa)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.2, outlier.alpha = 0.1) +
  scale_fill_manual(name = "", labels = c("Outside PAs", "Within PAs"),
                    values=c("lightskyblue2", "dodgerblue3")) +
  xlab("Variance in NDVI residuals") +
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
        legend.position = 'inside', legend.position.inside = c(0.9, 0.2))

fig_3 <-
  ggarrange(parkyvar, parkdoyvar, nrow = 2, labels = c("B", "C")) %>%
  ggarrange(p_var_canada, ., ncol = 2, labels = c("A", ""),
            widths = c(0.75, 0.5)) %>%
  ggarrange(boxvar, nrow = 2, labels = c("", "D"), heights = c(0.75, 0.5))

ggsave('Figures/figure-3.png', fig_3, units = "in", width = 6.86,
       height = 8.52, bg = "white", dpi = 600)

# check years with oddly high variance
filter(s2_year, s2 > 0.015)

#plot quantiles (figure 6) ---------------------------------------------------------
#' *NOTE: converting negative CV to Inf since CV --> Inf as NDVI --> 0*
plot(seq(-0.1, 1, by = 1e-3), 0.01 / seq(-0.1, 1, by = 1e-3), type = 'l',
     col = 'red', lwd = 5, xlab = 'Mean NDVI', ylab = 'CV of NDVI')

# add columns of species richness and 
est_albers <-
  mutate(est_albers,
         rich = extract(rast('Data/species_richness.tif'),
                        tibble(x, y))[, 2],
         extr = extract(rast('Outputs/n-extreme-months-1981-2024.tif'),
                        tibble(x, y))[, 2])

#calculate top and bottom 30th percentiles
cutoffs <- est_albers %>%
  mutate(cv_hat = if_else(cv_hat < 0, Inf, cv_hat)) %>%
  summarize(top_mean = quantile(mu_hat, probs = 0.7),
            bottom_var = quantile(s2_hat, probs = 0.3),
            bottom_cv = quantile(cv_hat, probs = 0.3),
            top_rich = quantile(rich, probs = 0.7, na.rm = TRUE),
            bottom_extr = quantile(extr, probs = 0.3, na.rm = TRUE))

# identify areas in each quantile group
areas <-
  est_albers %>%
  transmute(x, y,
            top_mean = if_else(mu_hat >= cutoffs$top_mean, mu_hat, NA_real_),
            bottom_var = if_else(s2_hat <= cutoffs$bottom_var, s2_hat, NA_real_),
            bottom_cv = if_else(cv_hat <= cutoffs$bottom_cv, cv_hat, NA_real_),
            top_rich = if_else(rich >= cutoffs$top_rich, rich, NA_real_),
            bottom_extr = if_else(extr <= cutoffs$bottom_extr, extr, NA_real_))

# save quantiles as rasters for ease of access
r_areas <- rast(areas, crs = canada_albers)
plot(r_areas)

mask(as.numeric(! is.na(r_areas$top_mean)) +
       as.numeric(! is.na(r_areas$bottom_cv)),
     st_transform(canada, crs(r_areas))) %>%
  plot()

writeRaster(r_areas$top_mean, 'Outputs/top-30-percent-mean-ndvi.tif')
writeRaster(r_areas$bottom_var, 'Outputs/bottom-30-percent-variance-ndvi.tif')
writeRaster(r_areas$bottom_cv, 'Outputs/bottom-30-percent-cv-ndvi.tif')

# make figure 6
theme_6 <-
  ggplot() +
  theme_void() +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.1, 0.95),
        legend.justification = 'left',
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) +
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 6, barheight = 0.5,
                               direction = "horizontal"))
mean.quant <-
  theme_6 +
  geom_sf(data = canada, color = "black", fill = 'grey50', lwd = 0.2) + 
  geom_raster(aes(x, y, fill = top_mean), areas) +
  geom_sf(data = canada, color = "black", fill = NA, lwd = 0.2) + 
  scale_fill_gradientn(name = 'Mean NDVI', colours = NDVI_cols,
                       limits = range(est_albers$mu_hat))

var.quant <-
  theme_6 +
  geom_sf(data = canada, color = "black", fill = 'grey50', lwd = 0.2) + 
  geom_raster(aes(x, y, fill = bottom_var), areas) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_devon(name = 'Variance in NDVI', reverse = TRUE,
                   limits = c(0, 0.04))

cv.quant <-
  theme_6 +
  geom_sf(data = canada, color = "black", fill = 'grey50', lwd = 0.2) + 
  geom_raster(aes(x, y, fill = bottom_cv), areas) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_acton(name = 'CV in NDVI', reverse = TRUE,
                   limits = c(0, 0.5))

p_spp_richness <-
  theme_6 +
  geom_sf(data = canada, color = "black", fill = 'grey50', lwd = 0.2) + 
  geom_raster(aes(x, y, fill = top_rich), areas) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) +
  scale_fill_batlow(name = 'Species richness', reverse = TRUE,
                    limits = c(0, NA))

p_extr <-
  theme_6 +
  geom_sf(data = canada, color = "black", fill = 'grey50', lwd = 0.2) + 
  geom_raster(aes(x, y, fill = bottom_extr), areas) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) +
  scale_fill_lajolla(name = 'Months with extreme temperatures',
                     limits = c(24, 164))

p_pas <- ggplot() +
  theme_void() +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.1, 0.95),
        legend.justification = 'left',
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) +
  geom_raster(aes(x, y, fill = pa),
              as.data.frame(r_pas, xy = TRUE) %>%
                mutate(pa = if_else(layer == 1, 'Yes', 'No') %>%
                         factor(levels = c('Yes', 'No')))) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) +
  scale_fill_manual('Protected areas', values = c('#255016', 'grey50'))

fig_6 <- plot_grid(p_pas, mean.quant, var.quant,
                   cv.quant, p_spp_richness, p_extr,
                   nrow = 2, labels = "AUTO")

ggsave('Figures/figure-6.png', fig_6, units = "in", bg = "white",
       width = 12.75, height = 7, dpi = 600)

#plot model residuals (Fig S1 in appendix) --------------------------------

#spatial mean of residuals
p_res_canada <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = mean_e), est_albers) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_distiller("Mean NDVI residuals", type = 'div',
                       palette = 5, limits = c(-0.2, 0.2)) +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 18, barheight = 0.5,
                               direction = "horizontal"))

# hex plots of residuals vs other variables (using only a subset)
if(file.exists('Data_annotated/ndvi-data-with-fitted-and-e-1-percent-subset.rds')) {
  d_sub <- readRDS('Data_annotated/ndvi-data-with-fitted-and-e-1-percent-subset.rds')
} else {
  d_sub <- slice_sample(d, prop = 0.01)
  saveRDS(d_sub, 'Data_annotated/ndvi-data-with-fitted-and-e-1-percent-subset.rds')
}
nrow(d_sub)

p_hex <-
  ggplot(d_sub) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(face = 'bold'),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "top") +
  ylab("NDVI residuals") +
  scale_fill_lapaz(name = 'Count', reverse = TRUE,
                   labels = function(x) {
                     .e <- floor(log10(x))
                     paste0(round(x / 10^.e), 'e', .e)
                   }) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 6, barheight = 0.5,
                               direction = "horizontal"))

p_res_doy <-
  p_hex +
  geom_hex(aes(doy, e)) +
  geom_hline(yintercept = 0, lwd = 0.1) +
  xlab("Day of year")

p_res_year <-
  p_hex +
  geom_hex(aes(year, e)) +
  geom_hline(yintercept = 0, lwd = 0.1) +
  geom_vline(xintercept = 2014, lwd = 0.1, lty = 'dashed') +
  xlab("Year")

p_res_lat <-
  p_hex +
  geom_hex(aes(y, e)) +
  geom_hline(yintercept = 0, lwd = 0.1) +
  xlab("Latitude")

p_res_long <-
  p_hex +
  geom_hex(aes(x, e)) +
  geom_hline(yintercept = 0, lwd = 0.1) +
  xlab("Longitude")

p_res_elev <-
  p_hex +
  geom_hex(aes(elev_m, e)) +
  geom_hline(yintercept = 0, lwd = 0.1) +
  xlab("Elevation (m)")

p_s1 <-
  ggarrange(
    ggarrange(p_res_canada,
              ggarrange(p_res_doy, p_res_year, ncol = 2,
                        labels = c("D", "E")),
              nrow = 2, labels = c("A", ""), heights = c(1, 0.5)),
    ggarrange(p_res_lat, p_res_long, p_res_elev, nrow = 3,
              labels = c("B", "C", "F")),
    ncol = 2, widths = c(1, 0.5))

ggsave('Figures/figure-s1-residuals.png', units = "in", bg = "white", width = 6.86,
       height = 6.86, dpi = 600)

# histogram of model residuals (Fig. S2 in appendix) ----
p_s2 <-
  d_sub %>%
  ggplot() +
  geom_histogram(aes(e), fill = 'grey', color = 'black', binwidth = 0.025,
                 center = 0) +
  xlab("NDVI Residuals") +
  ylab('Count') +
  scale_y_continuous(expand = c(0.025, 0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"))

ggsave('Figures/figure-s2-residuals-hist.png', p_s2, units = "in",
       bg = "white", width = 6.86, height = 4.5, dpi = 600)

#correlation between mean and variance (Fig. S3 in appendix) --------------
p_s3 <-
  ggplot(est_albers, aes(x = mu_hat, y = s2_hat)) +
  geom_hex(bins = 75) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=8, family = "sans", face = "bold"),
        axis.text = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12,
                                  family = "sans", face = "bold"),
        legend.position = "top",
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"),
        legend.title = element_text(face = "bold")) +
  scale_x_continuous('Mean NDVI') +
  scale_y_continuous('Variance in NDVI') +
  khroma::scale_fill_lapaz(name = 'Count', reverse = TRUE,
                           breaks = c(1, 500, 1000, 1500)) +
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 15, barheight = 0.5,
                               direction = "horizontal"))

ggsave("Figures/figure-s3-mean-variance-correlation.png", p_s3,
       width = 6.86, height = 6.86, units = "in", dpi = 600, bg = "white")
