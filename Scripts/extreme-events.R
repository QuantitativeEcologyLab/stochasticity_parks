library('dplyr')      # for data wrangling
library('sf')         # for spatial data
library('climatenaR') # for downloading historical climate data
library('purrr')      # for functional programming
library('tidyr')      # for data wrangling
library('terra')      # for raster data (masks tidyr::extract())
library('ggplot2')    # for fancy plots
library('cowplot')    # for fancy plots in grids
library('ggspatial')  # for maps in ggplot2

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
          plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), 'cm'), #top, left, bottom, right
          legend.position = 'inside',
          legend.position.inside = c(0.25, 0.9)))

prov <- canadianmaps::PROV %>%
  st_geometry() %>%
  st_as_sf()

mu_r <- rast('Outputs/mean_predNDVI_raster.nc')
s2_r <- rast('Outputs/var_residuals_raster.nc')
cv_r <- rast('Outputs/coeff_variation_raster.nc')

plot(mu_r)
plot(prov, add = TRUE)

dem <- elevatr::get_elev_raster(locations = mu_r, z = 1) %>%
  rast() %>%
  crop(prov) %>%
  mask(prov)
plot(dem)
plot(prov, add = TRUE)

#' change the working directory as required by `climatenaR`
#' this directory should have all the files for climateNA to run
setwd('climatena')
writeRaster(dem, 'can-dem-z1.tif')
plot(rast('can-dem-z1.tif'))
plot(prov, add = TRUE)

if(! file.exists('can-dem-z1.csv')) {
  #' convert the can DEM to a csv as required by `climatenaR`
  demToCSV(file = 'can-dem-z1.tif',
           outdir = '.', # save in current folder
           srs = NULL) # keep NULL if in lat/long
  
  # check the csv
  read.csv('can-dem-z1.csv', nrows = 5)
  
  read.csv('can-dem-z1.csv') %>%
    ggplot() +
    geom_raster(aes(-long, lat, fill = el))
}

if(! dir.exists('can-dem-z1')) dir.create('can-dem-z1')

# downloading all possible historical data (2024-09-17)
map(2022:1901,
    \(y) {
      cat(paste0('Downloading estimated historical data for ', y, '...\n'))
      histClimateNA(file = 'can-dem-z1.csv',
                    dateR = as.character(y), # year
                    tFrame = 'M', # monthly averages are the finest scale
                    exe = 'ClimateNA_v7.42.exe', # must be in wd
                    outdir = 'can-dem-z1')
    })

# single mean and sd for each pixel to get number of events outside mean +/- 2SE
setwd('..')
COLS <- c('Latitude', 'Longitude', 'Elevation', 'Tave01', 'Tave02',
          'Tave03', 'Tave04', 'Tave05', 'Tave06', 'Tave07', 'Tave08',
          'Tave09', 'Tave10', 'Tave11', 'Tave12')

extremes <-
  map_dfr(
    list.files('climatena/can-dem-z1', full.names = TRUE, pattern = 'csv'),
    \(.f) data.table::fread(.f, na.strings = '-9999', select = c(COLS))) %>%
  filter(! is.na(Tave01)) %>%
  pivot_longer(Tave01:Tave12, values_to = 'temp_C', names_to = 'month') %>%
  group_by(Latitude, Longitude, Elevation, month) %>%
  summarize(mu = mean(temp_C),
            sd = sd(temp_C),
            n_extr = sum(abs((temp_C - mu) / sd) >= qnorm(0.975))) %>%
  summarize(n_extr = sum(n_extr), .groups = 'drop') %>%
  mutate(mu = extract(mu_r, data.frame(Longitude, Latitude))[, 2],
         s2 = extract(s2_r, data.frame(Longitude, Latitude))[, 2],
         cv = extract(cv_r, data.frame(Longitude, Latitude))[, 2])

extremes_rast <- rast(select(extremes, Longitude, Latitude, n_extr)) %>%
  `crs<-`('EPSG:4326') %>%
  project('ESRI:102001')
plot(extremes_rast)
extremes_rast <- as.data.frame(extremes_rast, xy = TRUE)

# get summary statistics
sum(extremes$n_extr)
range(extremes$n_extr)

map <-
  ggplot(extremes_rast) +
  geom_raster(aes(x, y, fill = n_extr)) +
  geom_sf(data = st_transform(prov, 'ESRI:102001'), fill = 'transparent',
          color = 'black', lwd = 0.5) +
  khroma::scale_fill_acton(name = 'Number of extreme temperature events',
                           reverse = TRUE, limits = c(50, 90)) +
  theme_map() +
  theme(legend.position = 'top', legend.justification = 'center',
        panel.grid = element_blank(),
        legend.key.width = rel(2),
        text = element_text(face = 'bold'))

# scatterplots of extreme events by mean, variance, and CV
scatters <-
  extremes %>%
  pivot_longer(c(mu, s2, cv)) %>%
  filter(! is.na(value)) %>%
  mutate(name = case_when(name == 'mu' ~ 'Mean NDVI',
                          name == 's2' ~ 'Variance in NDVI',
                          name == 'cv' ~ 'Coefficient of variation') %>%
           factor(., levels = unique(.))) %>%
  ggplot(aes(value, n_extr)) +
  facet_wrap(~ name, nrow = 1, scales = 'free_x', strip.position = 'bottom') +
  geom_jitter(alpha = 0.05, width = 0, height = 0.25, size = 0.1) +
  geom_smooth(aes(color = name, fill = name),
              method = 'gam', formula = y ~ s(x, k = 4),
              method.args = list(family = poisson(), method = 'REML'),
              show.legend = FALSE) +
  labs(x = NULL, y = 'Number of extreme\ntemperature events') +
  scale_color_manual(values = c('forestgreen', 'dodgerblue3', '#93799a'),
                     aesthetics = c('color', 'fill')) +
  theme(legend.position = 'top', panel.grid = element_blank(),
        strip.placement = 'outside', strip.background = element_blank(),
        strip.text = element_text(size = 11),
        text = element_text(face = 'bold'))

# plot the two together
plot_grid(map, scatters, ncol = 1, rel_heights = c(3, 1),
          labels = c('a', 'b'))

ggsave('Figures/canada-extreme-events.png', width = 8, height = 9,
       scale = 1.1, dpi = 600, bg = 'white')
