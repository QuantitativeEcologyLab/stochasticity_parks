library('ggplot2') #for fancy figures
library('ggpubr') #arrange multiple figures
library('lubridate') #for decimal_date function
library('mgcv') #for predict.gam
library('sf') #for importing ecozone data
library('terra') #for basemaps
library('ggspatial')
library('dplyr') # for data wrangling
library('khroma') # for color palettes
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
plot(ecozones)

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
plot(r_pas)

# import model, data with residuals, and projected spatial summary
m_beta <- readRDS('Models/canada-mean-ndvi-aggr-2-2025-07-31-beta.rds')
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
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.2)

fig_1 <- ggarrange(parks.plot, ecoz.plot, nrow = 1, common.legend = T,
                   labels = "AUTO")

ggsave("Figures/figure1.png", fig_1, units = "in", width = 6.46,
       height = 3, bg = "white", dpi = 600, scale = 2); beepr::beep(4)

#plot mean trends (figure 2) -----------------------------------------------------
#spatial mean trends
p_mean_canada <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = 'transparent') + 
  geom_raster(aes(x, y, fill = mu_hat), est_albers) +
  geom_sf(data = canada, fill = NA, color = "black") + 
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

ggsave('Figures/figure2.png', fig_2, units = "in", width = 6.86,
       height = 8.52, bg = "white", dpi = 600)

#plot variance trends (figure 3) -------------------------------------------------
#spatial variance trends
#' **HERE**
hist(est_albers$s2_hat)

p_var_canada <-
  est_albers %>%
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "black") +
  geom_raster(aes(x, y, fill = s2_hat)) +
  geom_sf(data = canada, fill = 'transparent', color = "black") +
  scale_fill_devon(name = 'Variance in NDVI', limits = c(0, NA),
                   reverse = TRUE) +
  coord_sf(datum = canada_albers) +
  theme_void() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(
    title.position = "top", ticks.colour = NA, barwidth = 18,
    barheight = 0.5, title = "Variance in NDVI",  direction = "horizontal"))

#separate data by within/outside parks
parkyvar <-
  d %>%
  summarize(s2 = mean(e^2, na.rm = TRUE), .by = c(pa, year, x, y)) %>%
  summarize(s2 = mean(s2, na.rm = TRUE), .by = c(pa, year)) %>%
  mutate(pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
           factor()) %>%
  ggplot() +
  geom_hline(yintercept = INT_S2, color = 'grey', lty = 'dashed') +
  geom_point(aes(year, s2, color = pa)) +
  geom_smooth(aes(year, s2, colour = pa, fill = pa), method = 'gam') +
  scale_colour_manual(name = NULL, aesthetics = c('color', 'fill'),
                      values = c('lightskyblue2', 'dodgerblue3')) +
  xlab("Year") +
  ylab("Variance in NDVI") +
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
        legend.position = 'inside', legend.position.inside = c(0.25, 0.9))

#residual trends by day of year
parkdoyvar <- 
  d %>%
  summarize(s2 = mean(e^2, na.rm = TRUE), .by = c(pa, doy, x, y)) %>%
  summarize(s2 = mean(s2, na.rm = TRUE), .by = c(pa, doy)) %>%
  mutate(pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
           factor()) %>%
  ggplot() +
  geom_hline(yintercept = INT_S2, color = 'grey', lty = 'dashed') +
  geom_point(aes(doy, s2, color = pa)) +
  geom_smooth(aes(doy, s2, colour = pa, fill = pa), method = 'gam') +
  scale_colour_manual(name = NULL, aesthetics = c('color', 'fill'),
                      values = c('lightskyblue2', 'dodgerblue3')) +
  xlab("Day of year") +
  ylab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

#plot mean trends by ecozone

#boxplot to plot variance in different ecozones
boxvar <-
  est_albers %>%
  filter(! is.na(ecozone)) %>%
  mutate(median_ndvi = median(s2_hat), .by = ecozone) %>%
  arrange(median_ndvi) %>%
  mutate(ecozone = factor(ecozone, levels = unique(ecozone)),
         pa = if_else(pa == 0, "Outside PAs", "Within PAs") %>%
           factor()) %>%
  ggplot(aes(s2_hat, ecozone, group = paste(ecozone, pa), fill = pa)) +
  geom_vline(xintercept = INT_S2, color = 'grey', lty = 'dashed') +
  geom_boxplot(outlier.size = 0.1, lwd = 0.2, outlier.alpha = 0.1) +
  scale_fill_manual(name = "", labels = c("Outside PAs", "Within PAs"),
                    values=c("lightskyblue2", "dodgerblue3")) +
  xlab("Variance in response residuals") +
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
        legend.position = 'inside', legend.position.inside = c(0.9, 0.9))

fig_3 <-
  ggarrange(parkyvar, parkdoyvar, nrow = 2, labels = c("B", "C")) %>%
  ggarrange(p_var_canada, ., ncol = 2, labels = c("A", ""),
            widths = c(0.75, 0.5)) %>%
  ggarrange(boxvar, nrow = 2, labels = c("", "D"), heights = c(0.75, 0.5))

ggsave('Figures/figure3.png', fig_3, units = "in", width = 6.86,
       height = 8.52, bg = "white", dpi = 600)

#plot quantiles (figure 4) ---------------------------------------------------------
#' *NOTE: using absolute CV to avoid issues with NDVI < 0 and negative CV*
#' *not sure if this is the best choice*
plot(seq(-0.1, 1, by = 1e-3), 0.01 / seq(-0.1, 1, by = 1e-3), type = 'l',
     col = 'red', lwd = 5, xlab = 'Mean NDVI', ylab = 'CV of NDVI')
plot(seq(-0.1, 1, by = 1e-3), 0.01 / ((seq(-0.1, 1, by = 1e-3) + 0.1) / 1.1),
     type = 'l', col = 'red', lwd = 5, xlab = 'Mean NDVI',
     ylab = 'CV of transformed NDVI')

#calculate top and bottom 30th percentiles
cutoffs <- est_albers %>%
  mutate(cv = if_else(mu_hat < 0, Inf, cv)) %>%
  summarize(top_mean = quantile(mu_hat, probs = 0.7),
            bottom_var = quantile(s2_hat, probs = 0.3, na.rm = TRUE),
            bottom_cv = quantile(abs(cv), probs = 0.3))

# identify areas in each quantile group
areas <-
  est_albers %>%
  transmute(x, y,
            top_mean = mu_hat >= cutoffs$top_mean,
            bottom_var = s2_hat <= cutoffs$bottom_var,
            bottom_abs_cv = abs(cv) <= cutoffs$bottom_cv)

# save quantiles as rasters for ease of access
r_areas <- rast(areas)
plot(r_areas)

plot(r_areas$top_mean + r_areas$bottom_abs_cv)

# writeRaster(r_areas$top_mean, 'Canada/top-30-percent-mean-ndvi.tif')
# writeRaster(r_areas$bottom_var, 'Canada/bottom-30-percent-variance.tif')
# writeRaster(r_areas$bottom_abs_cv, 'Canada/bottom-30-percent-cv.tif')

mean.quant <-
  ggplot() +
  geom_sf(data = canada, color = "black", fill = 'grey50') + 
  geom_raster(aes(x, y, fill = mu_hat),
              filter(est_albers, mu_hat >= cutoffs$top_mean)) +
  geom_sf(data = canada, color = "black", fill = NA) + 
  scale_fill_gradientn(name = 'Mean NDVI', colours = NDVI_cols,
                       limits = range(est_albers$mu_hat)) + 
  theme_void() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0.1, 0.95),
        legend.justification = 'left',
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

var.quant <-
  ggplot() +
  geom_sf(data = canada, color = "black", fill = 'grey50') + 
  geom_raster(aes(x, y, fill = s2_hat),
              filter(est_albers, s2_hat <= cutoffs$bottom_var)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_devon(name = 'Variance in NDVI', reverse = TRUE,
                   limits = c(0, max(est_albers$s2_hat))) + 
  theme_void() +
  theme(legend.position = c(0.1, 0.95),
        legend.justification = 'left',
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

cv.quant <-
  est_albers %>%
  mutate(cv = if_else(mu_hat < 0, Inf, cv)) %>%
  filter(cv <= cutoffs$bottom_cv) %>%
  ggplot() +
  geom_sf(data = canada, color = "black", fill = 'grey50') + 
  geom_raster(aes(x, y, fill = cv)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_acton(name = 'Absolute CV in NDVI', reverse = TRUE,
                   limits = c(0, max(abs(est_albers$cv)))) +
  theme_void() +
  theme(legend.position = c(0.1, 0.95),
        legend.justification = 'left',
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

p_spp_richness <-
  ggplot() +
  geom_sf(data = canada, color = "black", fill = 'grey50') + 
  geom_raster(aes(x, y, fill = rich),
              mutate(est_albers, rich = rpois(n(), 10)) %>%
                filter(rich >= quantile(rich, 0.7))) +
  geom_sf(data = canada, fill = NA, color = "black") +
  scale_fill_batlow(name = 'Species richness', reverse = TRUE,
                    limits = c(0, NA)) +
  theme_void() +
  theme(legend.position = c(0.1, 0.95),
        legend.justification = 'left',
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

p_quantiles <- ggarrange(mean.quant, var.quant, cv.quant, p_spp_richness,
                         nrow = 2, ncol = 2, labels = "AUTO")

ggsave('Figures/figure4.png', p_quantiles, units = "in", bg = "white",
       width = 8.5, height = 8.5, dpi = 600)

#plot model residuals (Fig S1 in appendix) --------------------------------

#spatial mean of residuals
p_res_canada <-
  ggplot() +
  geom_raster(aes(x, y, fill = mean_e), est_albers) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_distiller("Mean response residuals", type = 'div',
                       palette = 5, limits = c(-0.2, 0.2)) +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 18, barheight = 0.5,
                               direction = "horizontal"))

#mean of residuals for each area - plots faster than all residuals
p_res_doy <-
  ggplot(d, aes(doy, e)) +
  geom_hex() +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), color = 'black') +
  xlab("Day of year") +
  ylab("Response residuals") +
  scale_fill_lapaz(name = 'Count', reverse = TRUE, limits = c(0, NA)) + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(face = 'bold'),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "top")

p_res_year <-
  ggplot(d, aes(year, e)) +
  geom_hex() +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), color = 'black') +
  xlab("Year") +
  ylab("Response residuals") + 
  scale_fill_lapaz(name = 'Count', reverse = TRUE, limits = c(0, NA)) + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(face = 'bold'),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "top")

p_res_lat <-
  ggplot(d, aes(y, e)) +
  geom_hex() +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), color = 'black') +
  xlab("Latitude") +
  ylab("Response residuals") + 
  scale_fill_lapaz(name = 'Count', reverse = TRUE, limits = c(0, NA)) + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(face = 'bold'),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "top")

p_res_long <-
  ggplot(d, aes(x, e)) +
  geom_hex() +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), color = 'black') +
  xlab("Longitude") +
  ylab("Response residuals") + 
  scale_fill_lapaz(name = 'Count', reverse = TRUE, limits = c(0, NA)) + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(face = 'bold'),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "top")

p_res_elev <-
  ggplot(d, aes(elev_m, e)) +
  geom_hex() +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5), color = 'black') +
  xlab("Elevation (m)") +
  ylab("Response residuals") + 
  scale_fill_lapaz(name = 'Count', reverse = TRUE, limits = c(0, NA)) + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(face = 'bold'),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "top")

ggarrange(
  ggarrange(p_res_canada,
            ggarrange(p_res_doy, p_res_year, ncol = 2,
                      labels = c("D", "E")),
            nrow = 2, labels = c("A", ""), heights = c(1, 0.5)),
  ggarrange(p_res_lat, p_res_long, p_res_elev, nrow = 3,
            labels = c("B", "C", "F")),
  ncol = 2, widths = c(1, 0.5))

ggsave('Figures/residuals.png', units = "in", bg = "white", width = 6.86,
       height = 6.86, dpi = 600)

# histogram of model residuals (Fig. S2 in appendix) ----
p_res_hist <-
  tibble(res = rnorm(10)) %>%
  ggplot() +
  geom_histogram(aes(res), fill = 'grey', color = 'black', binwidth = 0.25,
                 center = 0) +
  xlab("Response Residuals") +
  ylab('Count') +
  scale_y_continuous(expand = c(0.025, 0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"))

ggsave('Figures/residuals-hist.png', p_res_hist, units = "in",
       bg = "white", width = 6.86, height = 4.5, dpi = 600)

#coefficient of variation (currently not in appendix) ---------------------

# #spatial variance trends
# p_cv_canada <- 
#   ggplot() +
#   geom_raster(est_albers, mapping = aes(x, y, fill = cv)) +
#   geom_sf(data = canada, fill = NA, color = "black") + 
#   scale_fill_acton(name = 'Coefficient of Variation') +
#   theme_void() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
#         legend.position = "bottom",
#         legend.title = element_text(face = "bold")) + 
#   guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 18,
#                                barheight = 0.5, title = "Coefficient of Variation",  direction = "horizontal"))
# 
# #plot smooths
# 
# #separate data by within/outside parks
# parkycv <- 
#   ggplot() +
#   geom_smooth(aes(year, cv, colour = pa), span = 0.25, se = FALSE) +
#   geom_smooth(aes(year, cv, colour = "Outside Parks"),  span = 0.25, se = FALSE) +
#   scale_colour_manual(name = "", values = c("#ccb3d1", "#5b3f64")) +
#   xlab("Year") +
#   ylab("Coefficient of Variation") +
#   scale_x_continuous(expand = c(0,0)) +
#   theme_classic() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
#         axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
#         axis.text.y = element_text(size = 8, family = "sans"),
#         axis.text.x = element_text(size = 8, family = "sans"),
#         legend.background = element_rect(fill = "transparent"),
#         legend.key.size = unit(0.5, 'cm'),
#         legend.spacing.y = unit(0.1, 'cm'),
#         legend.text = element_text(size = 6, family = "sans", face = "bold"),
#         panel.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
#         legend.position = c(0.25, 0.2))
# 
# #residual trends by day of year
# parkdoycv <- 
#   ggplot() +
#   geom_point(doyparkyes, mapping = aes(doy, cv), size = 0.3, color = "#5b3f64") +
#   geom_point(doyparkno, mapping = aes(doy, cv), size = 0.3, color = "#ccb3d1") +
#   geom_smooth(doyparkyes, mapping = aes(doy, cv, colour = "Within Parks"), span = 0.15, se = FALSE) +
#   geom_smooth(doyparkno, mapping = aes(doy, cv, colour = "Outside Parks"),  span = 0.15, se = FALSE) +
#   scale_colour_manual(name = "", values = c("#ccb3d1", "#5b3f64")) +
#   xlab("Day of year") +
#   ylab("Coefficient of Variation") +
#   scale_x_continuous(expand = c(0,0)) +
#   theme_classic() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
#         axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
#         axis.text.y = element_text(size = 8, family = "sans"),
#         axis.text.x = element_text(size = 8, family = "sans"),
#         panel.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
#         legend.position = "none")
# 
# #plot mean trends by ecozone
# 
# #boxplot to plot mean in different ecozones
# boxcv <- 
#   ggplot(edata, aes(x = layer, y = cv, fill = park)) +
#   geom_boxplot(outlier.size = 0.3, lwd = 0.2) +
#   scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values = c("#ccb3d1", "#5b3f64")) +
#   scale_x_discrete(labels=c("2" = "Northern Arctic", "1" = "Arctic Cordillera", "3" = "Southern Arctic", "9" = "Boreal Plain",
#                             "6" = "Boreal Shield", "15" = "Hudson Plain", "14" = "Montane Cordillera", "8" = "Mixedwood Plain",
#                             "13" = "Pacific Maritime", "10" = "Prairies", "11" = "Taiga Cordillera", "4" = "Taiga Plain",
#                             "5" = "Taiga Shield", "7" = "Atlantic Maritime", "12" = "Boreal Cordillera")) +
#   xlab("Ecozone") +
#   ylab("Coefficient of Variation") +
#   theme_classic() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
#         axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
#         axis.text.y = element_text(size = 8, family = "sans"),
#         axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1, size = 8, family = "sans"),
#         legend.background = element_rect(fill = "transparent"),
#         legend.key.size = unit(0.7, 'cm'),
#         legend.spacing.y = unit(0.1, 'cm'),
#         legend.text = element_text(size=6, family = "sans", face = "bold"),
#         panel.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
#         legend.position = c(0.1, 0.9))
# 
# ggarrange(
#   ggarrange(
#     cv, ggarrange(parkycv, parkdoycv, nrow = 2, labels = c("B", "C")),
#     ncol = 2, labels = c("A", ""), widths = c(0.75, 0.5)),
#   boxcv, nrow = 2, labels = c("", "D"), heights = c(0.75, 0.5))
# 
# ggsave('cv.png', units = "in", width = 6.86, height = 8.52, bg = "white",
#        dpi = 600)

#correlation between mean and variance (Fig. S3 in appendix) --------------

p_correl <- 
  ggplot(est_albers, aes(x = mu_hat, y = s2_hat)) +
  geom_hex() +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5),
              method.args = list(family = Gamma('log')), color = 'black') +
  ylab("Variance in NDVI") +
  xlab("Mean NDVI") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=8, family = "sans", face = "bold"),
        axis.text = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "top",
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_y_continuous(expand = c(0.01,0)) +
  khroma::scale_fill_lapaz(name = 'Count', reverse = TRUE)

ggsave("Figures/Mean_Var_Correlation.png", p_correl,
       width = 6.86, height = 4.5, units = "in", dpi = 600, bg = "white")
