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
    summarize(q_0.025 = if_else(year < 1981, NA_real_, temp_C) %>%
                quantile(probs = 0.025, na.rm = TRUE),
              q_0.975 = if_else(year < 1981, NA_real_, temp_C) %>%
                quantile(probs = 0.975, na.rm = TRUE),
              n_extr = sum(temp_C < q_0.025 | temp_C > q_0.975),
              .by = c(Latitude, Longitude, Elevation, month)) %>%
    summarize(n_extr = sum(n_extr), .by = c(Latitude, Longitude, Elevation))
  
  saveRDS(extremes_longlat, 'Data/extreme-temperature-months-1981-2024.rds')
}

total_extremes <- sum(extremes_longlat$n_extr)
total_extremes # 4637257

# percentage of months with extreme temperatures
total_extremes / (12 * (2024 - 1981) * nrow(read.csv('ClimateNA_v760/can-dem.csv')))

range(extremes_longlat$n_extr) # 77 to 308
range(extremes_longlat$n_extr) / (12 * (2024 - 1981))

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
  mutate(s2 = if_else(s2 > 0.04, 0.04, s2))

map <-
  ggplot(extremes) +
  geom_sf(data = canada, fill = 'grey', color = 'transparent') +
  geom_raster(aes(x, y, fill = n_extr)) +
  geom_sf(data = canada, fill = 'transparent', color = 'black', lwd = 0.1) +
  scale_fill_lajolla(name = 'Number of extreme-temperature months',
                     limits = range(extremes_longlat$n_extr)) +
  theme_map() +
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 20, barheight = 0.5,
                               direction = "horizontal")) +
  theme(legend.position = 'top', legend.justification = 'center',
        panel.grid = element_blank(),
        text = element_text(face = 'bold'))

# scatterplots of extreme events by mean and variance
# ndvi can be < 0 and cv can be too large, so cv gives odd relationships
hexes <-
  ggplot() +
  stat_summary_hex(aes(mu, s2, z = n_extr), extremes, na.rm = TRUE,
                   bins = 75) +
  scale_fill_lajolla(limits = range(extremes_longlat$n_extr)) +
  labs(x = 'Mean NDVI', y = 'Variance in NDVI residuals') +
  theme(legend.position = 'none', text = element_text(face = 'bold'))

# plot the two together
fig_s4 <- plot_grid(map, hexes, ncol = 1, labels = c('A', 'B'))

ggsave('Figures/figure-s4-canada-extreme-events.png', fig_s4,
       width = 6.46, height = 11, dpi = 600, bg = 'white')
