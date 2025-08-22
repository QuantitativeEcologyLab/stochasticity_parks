library('dplyr')      # for data wrangling
library('sf')         # for spatial data
library('climatenaR') # for downloading historical climate data
library('purrr')      # for functional programming
library('tidyr')      # for data wrangling
library('terra')      # for raster data (masks tidyr::extract())
library('ggplot2')    # for fancy plots
library('cowplot')    # for fancy plots in grids
library('ggspatial')  # for maps in ggplot2
library('elevatr')    # for digital elevation models
library('khroma')     # for color palettes
library('mgcv')       # for GAMs
library('gratia')       # for GAMs

theme_set(
  theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size = 9, family = 'sans', face = 'bold'),
          axis.title.x = element_text(size = 9, family = 'sans', face = 'bold'),
          axis.text.y = element_text(size = 8, family = 'sans'),
          axis.text.x = element_text(size = 8, family = 'sans'),
          legend.background = element_rect(fill = 'transparent'),
          legend.key.size = unit(0.5, 'cm'),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.text = element_text(size = 6, family = 'sans', face = 'bold'),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), 'cm')))

canada_albers <- 'ESRI:102001'
canada <- rgeoboundaries::geoboundaries('Canada') %>%
  st_geometry() %>%
  st_as_sf() %>%
  st_transform(canada_albers)
plot(canada)

# rasters of estimated statistics of NDVI
r_mu <- rast('Outputs/estimated-mean.tif')
r_s2 <- rast('Outputs/estimated-variance.tif')
r_cv <- rast('Outputs/coefficient-of-variation.tif')

plot(r_mu)
plot(canada, add = TRUE)

# import/create a dataset of the number of months with extreme temperature
if(file.exists('Data/extreme-temperature-months-1981-2024.rds')) {
  extremes_longlat <- readRDS('Data/extreme-temperature-months-1981-2024.rds')
} else {
  # download a raster based on the raster of the spatially aggregated data
  r_0 <- list.files('Data/ndvi-data/AVHRR-Land_v005',
                    pattern = '.nc', full.names = TRUE) %>%
    first() %>%
    rast() %>%
    terra::aggregate(2)
  
  # using a coarse resolution since it's just a supporting analysis
  dem <- get_elev_raster(r_0, z = 1) %>%
    rast() %>%
    crop(., st_transform(canada, crs(.)), mask = TRUE)
  res(dem)
  plot(dem)
  plot(st_transform(canada, crs(dem)), add = TRUE)
  
  #' change the working directory as required by `climatenaR`
  #' this directory should have all the files for climateNA to run
  setwd('ClimateNA_v760')
  writeRaster(dem, 'can-dem.tif')
  plot(rast('can-dem.tif')) # check
  plot(st_transform(canada, 'EPSG:4326'), add = TRUE)
  
  #' convert the can DEM to a csv as required by `climatenaR`
  if(! file.exists('can-dem.csv')) {
    demToCSV(file = 'can-dem.tif',
             outdir = '.', # save in current folder
             srs = NULL) # keep NULL if in lat/long
    
    # check the csv
    print(read.csv('can-dem.csv', nrows = 5))
    
    read.csv('can-dem.csv') %>%
      ggplot() +
      geom_raster(aes(-long, lat, fill = el)) +
      coord_sf(crs = 'EPSG:4326')
  }
  
  if(! dir.exists('can-dem')) dir.create('can-dem')
  
  # downloading all possible historical data (2025-08-04)
  map(1901:1980,
      \(y) {
        cat(paste0('Downloading estimated historical data for ', y, '...\n'))
        histClimateNA(file = 'can-dem.csv',
                      dateR = as.character(y), # year
                      tFrame = 'M', # monthly averages are the finest scale
                      exe = 'ClimateNA_v7.60.exe', # must be in wd
                      outdir = 'can-dem')
      })
  
  # estimating extremes using a single mean and sd for each pixel to get
  # the number of events outside mean +/- 1.96 * SE
  setwd('..') # go back to main folder
  
  extremes_longlat <-
    map(1901:2024, function(.y) {
      .fn <- paste0('ClimateNA_v760/can-dem/can-dem_', .y, '.csv')
      data.table::fread(.fn, na.strings = '-9999',
                        select = c('Latitude', 'Longitude', 'Elevation',
                                   'Tave01', 'Tave02', 'Tave03', 'Tave04',
                                   'Tave05', 'Tave06', 'Tave07', 'Tave08',
                                   'Tave09', 'Tave10', 'Tave11', 'Tave12')) %>%
        mutate(year = .y) %>%
        return()
    }, .progress = TRUE) %>%
    bind_rows() %>%
    filter(! is.na(Tave01)) %>% # drop rows without temperature data
    pivot_longer(Tave01:Tave12, values_to = 'temp_C', names_to = 'month') %>%
    # get number of months with extreme temperatures for each year
    # using 1901-1980 as reference years
    group_by(Latitude, Longitude, Elevation, month) %>%
    mutate(
      # find lower and upper bounds using only data prior to 1981
      q_0.025 = if_else(year < 1981, temp_C, NA_real_) %>%
        quantile(probs = 0.025, na.rm = TRUE),
      q_0.975 = if_else(year < 1981, temp_C, NA_real_) %>%
        quantile(probs = 0.975, na.rm = TRUE),
      extreme = year >= 1981 & (temp_C < q_0.025 | temp_C > q_0.975)) %>%
    group_by(Latitude, Longitude, Elevation) %>%
    summarize(n_extr = sum(extreme), .groups = 'drop')
  
  saveRDS(extremes_longlat, 'Data/extreme-temperature-months-1981-2024.rds')
}

total_extremes <- sum(extremes_longlat$n_extr)
total_extremes # 1,747,106

# percentage of months with extreme temperatures (44 complete years of data)
n_locations <- nrow(read.csv('ClimateNA_v760/can-dem.csv'))
paste0(round(total_extremes / ((2025 - 1981) * 12 * n_locations) * 100, 1),
       '%')

range(extremes_longlat$n_extr) # 24 to 164
paste0(round(range(extremes_longlat$n_extr) / (12*(2024-1980)) * 100, 1),
       '%')

# project to Albers CRS
extremes_rast <-
  extremes_longlat %>%
  select(Longitude, Latitude, n_extr) %>%
  rast() %>%
  `crs<-`('EPSG:4326') %>%
  project(canada_albers) %>%
  round()

plot(extremes_rast)
plot(canada, add = TRUE)
writeRaster(extremes_rast, 'Outputs/n-extreme-months-1981-2024.tif')

extremes <-
  as.data.frame(extremes_rast, xy = TRUE, na.rm = TRUE) %>%
  mutate(mu = extract(r_mu, data.frame(x, y))[, 2],
         s2 = extract(r_s2, data.frame(x, y))[, 2]) %>%
  na.omit()

# drop variance values with too much leverage
hist(extremes$s2)
sum(extremes$s2 > 0.04, na.rm = TRUE)
paste(round(mean(extremes$s2 > 0.04, na.rm = TRUE) * 100, 2), '%')

extremes <- filter(extremes, s2 <= 0.04)

map <-
  ggplot(extremes) +
  geom_sf(data = canada, fill = 'grey', color = 'transparent') +
  geom_raster(aes(x, y, fill = n_extr)) +
  geom_sf(data = canada, fill = 'transparent', color = 'black', lwd = 0.1) +
  scale_fill_lajolla(name = 'Number of extreme-temperature months',
                     limits = range(extremes_longlat$n_extr)) +
  theme_map() +
  guides(fill = guide_colorbar(title.position = 'top', ticks.colour = NA,
                               barwidth = 20, barheight = 0.5,
                               direction = 'horizontal')) +
  theme(legend.position = 'top', legend.justification = 'center',
        panel.grid = element_blank(), text = element_text(face = 'bold'))

# scatterplots of extreme events by mean and variance
# ndvi can be < 0 and cv can be too large, so cv gives odd relationships
hexes_events <-
  ggplot() +
  stat_summary_hex(aes(mu, s2, z = n_extr), extremes, na.rm = TRUE,
                   bins = 75, fun = 'mean') +
  scale_fill_lajolla(name = 'Number of extreme-temperature months',
                     limits = range(extremes_longlat$n_extr)) +
  labs(x = 'Mean NDVI', y = 'Variance in NDVI') +
  theme(legend.position = 'top', legend.justification = 'center',
        panel.grid = element_blank(), text = element_text(face = 'bold')) +
  guides(fill = guide_colorbar(title.position = 'top', ticks.colour = NA,
                               barwidth = 20, barheight = 0.5,
                               direction = 'horizontal'))

hexes_n_obs <-
  ggplot() +
  stat_summary_hex(aes(mu, s2, z = n_extr), extremes, na.rm = TRUE,
                   bins = 75, fun = 'length') +
  scale_fill_lapaz(name = 'Number of observations', reverse = TRUE) +
  labs(x = 'Mean NDVI', y = 'Variance in NDVI') +
  theme(legend.position = 'top', legend.justification = 'center',
        panel.grid = element_blank(), text = element_text(face = 'bold')) +
  guides(fill = guide_colorbar(title.position = 'top', ticks.colour = NA,
                               barwidth = 20, barheight = 0.5,
                               direction = 'horizontal'))

# find the coordinates of the main peak
plot_grid(hexes_events +
            geom_hline(yintercept = 0.01) +
            geom_vline(xintercept = 0.05),
          hexes_n_obs +
            geom_hline(yintercept = 0.01) +
            geom_vline(xintercept = 0.05))

# the secondary peak has little data to support it
plot_grid(hexes_events +
            geom_hline(yintercept = 0.005) +
            geom_vline(xintercept = 0.15),
          hexes_n_obs +
            geom_hline(yintercept = 0.005) +
            geom_vline(xintercept = 0.15))

# plot the panels together
fig_s4 <- plot_grid(map,
                    plot_grid(hexes_events, hexes_n_obs, nrow = 1,
                              labels = c('B', 'C')),
                    ncol = 1, labels = c('A', ''))

ggsave('Figures/figure-s4-canada-extreme-events.png', fig_s4,
       width = 12.75, height = 11, dpi = 600, bg = 'white')

# fit a model with mean and variance in NDVI ----
m <-
  bam(
    n_extr ~ s(mu, k = 10) + s(sqrt(s2), k = 10) + ti(mu, sqrt(s2), k = 3),
    family = nb(link = 'log'),
    data = extremes,
    discrete = TRUE,
    method = 'fREML')
draw(m, n = 250, rug = FALSE, nrow = 1)
appraise(m, point_alpha = 0.3)
summary(m)

# predict from the model ----
preds_mu <-
  tibble(mu = seq(min(extremes$mu), max(extremes$mu), length.out = 500),
         s2 = 0) %>%
  bind_cols(.,
            predict(m, newdata = ., type = 'link', se = TRUE,
                    terms = c('(Intercept)', 's(mu)')) %>%
              as.data.frame()) %>%
  mutate(est = exp(fit),
         lwr_95 = exp(fit - 1.96 * se.fit),
         upr_95 = exp(fit + 1.96 * se.fit))

preds_s2 <-
  tibble(s2 = seq(min(extremes$s2), max(extremes$s2), length.out = 500),
         mu = 0) %>%
  bind_cols(.,
            predict(m, newdata = ., type = 'link', se = TRUE,
                    terms = c('(Intercept)', 's(sqrt(s2))')) %>%
              as.data.frame()) %>%
  mutate(est = exp(fit),
         lwr_95 = exp(fit - 1.96 * se.fit),
         upr_95 = exp(fit + 1.96 * se.fit))

preds_int <-
  expand_grid(mu = seq(min(extremes$mu), max(extremes$mu), length.out = 500),
              s2 = seq(min(extremes$s2), max(extremes$s2), length.out = 500)) %>%
  filter(! too_far(mu, s2, extremes$mu, extremes$s2, dist = 0.05)) %>%
  mutate(est = predict(m, newdata = ., type = 'response', se = FALSE))

# create the figure ----
# mean NDVI
fig_5A <-
  ggplot(preds_mu) +
  geom_ribbon(aes(mu, ymin = lwr_95, ymax = upr_95),
              fill = 'darkgreen',
              alpha = 0.3) +
  geom_line(data = preds_mu, aes(mu, est),
            colour = 'darkgreen', linewidth = 1) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = 'sans', face = 'bold'),
        axis.title.x = element_text(size=10, family = 'sans', face = 'bold'),
        axis.text.y  = element_text(size=10, family = 'sans', face = 'bold', color = 'black'),
        axis.text.x  = element_text(size=10, family = 'sans', face = 'bold', color = 'black'),
        plot.title = element_text(hjust = -0.05, size = 12, family = 'sans', face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = 'sans', face = 'bold'),
        legend.background = element_rect(fill = 'transparent'),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), 'cm')) +
  scale_x_continuous('Mean NDVI', expand = c(0,0.005)) +
  scale_y_continuous('Number of extreme-temperature months',
                     limits = c(0, NA))

# NDVI variance
fig_5B <- 
  ggplot(preds_s2) +
  geom_ribbon(aes(s2, ymin = lwr_95, ymax = upr_95),
              fill = 'dodgerblue3', alpha = 0.3) +
  geom_line(aes(s2, est), colour = 'dodgerblue3',
            linewidth = 1) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = 'sans', face = 'bold'),
        axis.title.x = element_text(size=10, family = 'sans', face = 'bold'),
        axis.text.y  = element_text(size=10, family = 'sans', face = 'bold', color = 'black'),
        axis.text.x  = element_text(size=10, family = 'sans', face = 'bold', color = 'black'),
        plot.title = element_text(hjust = -0.05, size = 12, family = 'sans', face = 'bold'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = 'sans', face = 'bold'),
        legend.background = element_rect(fill = 'transparent'),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), 'cm')) +
  scale_x_continuous('Variance in NDVI', expand = c(0,0.0004)) +
  scale_y_continuous('Number of extreme-temperature months',
                     limits = c(0, NA))

# interaction
fig_5C <-
  ggplot(preds_int) +
  geom_raster(aes(mu, s2, fill = est)) +
  geom_contour(aes(mu, s2, z = est), color = 'black') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_lajolla(
    name = expression(atop(bold('Number of months'),
                           bold('with extreme temperatures'))),
    limits = range(extremes_longlat$n_extr)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size=10, family = 'sans', face = 'bold'),
    axis.title.x = element_text(size=10, family = 'sans', face = 'bold'),
    axis.text.y = element_text(size=8, family = 'sans'),
    axis.text.x  = element_text(size=8, family = 'sans'),
    plot.title = element_text(hjust = -0.05, size = 12, family = 'sans', face = 'bold'),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'top',
    legend.title = element_text(size=8, family = 'sans', face = 'bold'),
    panel.background = element_rect(fill = 'grey40'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    plot.margin = unit(c(0.2,0.1,0.2,0.2), 'cm')) +
  xlab('Mean NDVI') +
  ylab('Variance in NDVI') +
  guides(fill = guide_colorbar(ticks.colour = NA, barwidth = 10,
                               barheight = 1, direction = 'horizontal'))

slice(preds_int, which.max(est))

# final figure ----
cowplot::plot_grid(fig_5A, fig_5B, fig_5C, labels = 'AUTO', nrow = 1)

ggsave('Figures/figure-5-extreme-events-regression.png',
       width = 6.86, height = (6.86 / 3), units = 'in',
       dpi = 600, bg = 'white', scale = 2)
