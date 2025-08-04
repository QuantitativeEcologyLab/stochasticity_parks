# code modified from Stefano Mezzini's project on estimating variance in
# NDVI globally: https://github.com/QuantitativeEcologyLab/ndvi-stochasticity
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('sf')        # for shapefiles
library('terra')     # for rasters
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('lubridate') # for working with dates
library('elevatr')   # for digital elevation models
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs (to estimate cutoffs)
source('functions/is_flagged.R') # this script sources other script
theme_set(theme_bw())

file_names <-
  list.files('../ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/',
             pattern = '.nc', full.names = TRUE, recursive = FALSE)
length(file_names) / 365 # approximate number of years of data

# all rasters have the same CRS
# check crs for first, last, and a random sample of rasters
if(FALSE) {
  file_names[c(1, length(file_names),
               sample(length(file_names), size = 100))] %>%
    map_chr(function(.fn) { # fast enough that it's not worth parallelizing
      crs(rast(.fn))
    }) %>%
    unique()
}

# import or create rasters of other variables
r_pa <- rast('Data/protected-areas/protected-areas-0-1.tif') # pa 0/1
r_ed <- rast('Data/ecodistricts/ecodistrict-id.tif') # ecodistricts
r_pw <- rast('Data/proportion-water.tif') # proportion water
if(file.exists('Data/elevation.tif')) {
  r_el <- rast('Data/elevation.tif')
} else {
  #' `z = 4` gives the resolution just finer than the NDVI rasters 
  r_el <- rast(get_elev_raster(r_pa, z = 4))
  res(r_el)
  res(r_pa)
  r_el <- project(r_el, r_pa) # move to same resolution
  r_el <- mask(r_el, r_pa) # drop sea pixels
  all(res(r_el) == res(r_pa))
  writeRaster(r_el, 'Data/elevation.tif')
}

if(FALSE) {
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  plot(r_pa); plot(r_ed); plot(r_pw); plot(r_el)
  layout(1)
}

# some coastlines have elevations < 0 m
plot(r_el > 0)

# import a data frame of all dates
if(file.exists('Data/ndvi-raster-metadata.rds')) {
  dates <- readRDS('Data/ndvi-raster-metadata.rds')
} else {
  plan(strategy = multisession,
       workers = availableCores(logical = FALSE) - 2)
  dates <-
    tibble(date = seq(as.Date('1981-06-24'), as.Date('2025-05-07'), by =1),
           file_name = future_map_chr(date, function(.date) {
             .fn <- file_names[grepl(format(.date, '_%Y%m%d_'),file_names)]
             
             if(length(.fn) == 1) {
               return(.fn)
             } else if(length(.fn) == 0) {
               return(NA_character_)
             } else {
               .msg <- paste0('Found ', length(.fn), ' files for ',
                              as.character(.date), '!')
               warning(.msg)
               return(.msg)
             }
             return()
           }, .progress = TRUE),
           n_cells = future_map_int(file_name, \(.fn) {
             if(file.exists(.fn)) {
               r <- rast(.fn)
               
               # drop probably/confidently cloudy pixels
               # bits 1 and 0 are 10 or 11
               r$NDVI <- ifel(is_flagged(r$QA, 1), NA, r$NDVI)
               
               not.na(r$NDVI) %>%
                 values() %>%
                 sum() %>%
                 as.integer() %>%
                 return()
             } else return(NA_integer_)
           }, .progress = TRUE))
  plan(sequential)
  saveRDS(dates, 'Data/ndvi-raster-metadata.rds')
}

# some dates are missing a raster (not available on the server)
all(! is.na(dates$file_name))
mean(! is.na(dates$file_name))
filter(dates, is.na(file_name))

# rasters are not saved in this repo to avoid having two copies
dates <- dates %>%
  filter(! is.na(file_name)) %>% # rm missing rasters
  mutate(file_name = paste0('../ndvi-stochasticity/', file_name))

# create the aggregated dataset ----
# spatRast objects moved across parallelized sessions, but data frames can
# https://stackoverflow.com/questions/67445883/terra-package-returns-error-when-try-to-run-parallel-operations/67449818#67449818
# shapefile of canada
shp <- st_read('Data/ecodistricts/Canada_Ecodistricts.shp') %>%
  st_geometry() %>%
  st_as_sf() %>%
  st_transform(crs(rast(file_names[1])))

# there are some bands of extremely high NDVI for high latitudes
NCORES <- min(availableCores() - 4, 12)
plan(multisession, workers = NCORES)

if(file.exists('Data/monthly-ecdf-summary.rds')) {
  d_ecdf <- readRDS('Data/monthly-ecdf-summary.rds') %>%
    mutate(month = factor(month.name[month], levels = month.name))
} else {
  d_ecdf <-
    dates %>%
    mutate(month = month(date)) %>%
    select(! n_cells) %>%
    nest(month_data = ! month) %>%
    mutate(month_data = future_map(month_data, \(.d) {
      .d <-
        map(.d$file_name, \(.fn) {
          .r <- rast(.fn) %>%
            # crop and mask to canada
            crop(shp) %>%
            mask(shp)
          # drop cloudy pixels
          .r$NDVI <- ifel(is_flagged(.r$QA, 1), NA, .r$NDVI)
          return(as.data.frame(.r$NDVI, xy = TRUE))
        }) %>%
        bind_rows() %>%
        summarize(cdf = list(ecdf(NDVI)), .by = y) %>%
        mutate(probs = map(cdf, \(.fun) {
          tibble(ndvi = seq(-0.1, 1, by = 0.01),
                 p = .fun(ndvi))
        })) %>%
        return()
    }, .progress = TRUE)) %>%
    unnest(month_data) %>%
    unnest(probs) %>%
    month = factor(month.name[month], levels = month.name)
    
    saveRDS(d_ecdf, 'Data/monthly-ecdf-summary.rds')
    
    # saving a version without ECDF to be able to push it to GitHub
    d_ecdf %>%
      select(! cdf) %>%
      saveRDS('Data/monthly-ecdf-summary-without-ecdf.rds')
}

# find NDVI values closest to the 99th percentile for each latitude
est_cutoffs <- d_ecdf %>%
  slice(which.min(abs(p - 0.99)), .by = c(month, y)) %>%
  # find NDVI cutoffs after dropping decimals in latitude
  mutate(month = factor(month.name[month], levels = month.name),
         y_2 = floor(y)) %>%
  mutate(ndvi_2 = min(ndvi), .by = c(y_2, month)) %>%
  #' arrange by latitude to make sure `geom_path()` looks right
  arrange(y)

ggplot(est_cutoffs) +
  facet_wrap(~ month) +
  geom_path(aes(ndvi_2, y)) +
  scale_x_continuous('NDVI', expand = c(0, 0)) +
  theme(panel.grid = element_blank(), legend.position = 'top',
        text = element_text(face = 'bold'))

cutoffs <- expand_grid(month = unique(est_cutoffs$month),
                       y = unique(est_cutoffs$y),
                       ndvi = c(0, 1)) %>%
  filter((month == 'January' & ! y < 60) |
           (month == 'February' & ! y < 65) |
           (month == 'March' & ! y < 70) |
           (month == 'September' & ! y < 80) |
           (month == 'October' & ! y < 69) |
           (month == 'November' & ! y < 62) |
           (month == 'December' & ! y < 59)) %>%
  summarize(ymin = min(y), ymax = max(y), .by = c(month, ndvi))

# filter months based on est_cutoffs (black lines or april ref line)
# bright bands simply won't have much data, but the bands are thin anyway
p_cutoffs <-
  ggplot() +
  facet_wrap(~ month) +
  geom_raster(aes(ndvi, y, fill = p), d_ecdf) +
  geom_ribbon(aes(x = ndvi, ymin = ymin, ymax = ymax), cutoffs,
              fill = '#FF000048') +
  geom_path(aes(ndvi, y), est_cutoffs, color = 'grey', linewidth = 0.5) +
  geom_path(aes(ndvi_2, y), est_cutoffs, color = 'white') +
  geom_vline(xintercept = c(0, 0.25), lty = 'dotted') +
  geom_hline(yintercept = c(60, 65, 70, 80), lwd = 0.1, lty = 'dotted') +
  scale_fill_viridis_c(expression(bold('ECDF(NDVI)')), direction = -1) +
  scale_x_continuous('NDVI', expand = c(0, 0)) +
  scale_y_continuous('Latitude', expand = c(0, 0),
                     breaks = c(50, 60, 65, 70, 80)) +
  theme(panel.grid = element_blank(), legend.position = 'top',
        text = element_text(face = 'bold')); p_cutoffs

# the change in cutoff throughout the year can be estimated by a function
# with very little uncertainty, but using these smooth changes risks
# leaving some extreme values in the data. it may be worth doing the check
# at intervals finer than a month (e.g., 7 or 14 days), but we can probably
# make the change in future versions, if necessary.
doy_cutoffs <- tibble(
  date = as.Date(c('2025-01-15', '2025-02-15', '2025-03-15', '2025-04-15',
                   '2025-05-15', '2025-06-15', '2025-07-15', '2025-08-15',
                   '2025-09-15', '2025-10-15', '2025-11-15', '2025-12-15')),
  doy = yday(date),
  month = month(date),
  y = c(60, 65, 70, rep(NA, 5), 80, 69, 62, 59))

m_cutoffs <- gam(y ~ s(doy, bs = 'cc', k = 7), data = doy_cutoffs,
                 method = 'REML', knots = list(doy = c(0.5, 365.5)))
plot.gam(m_cutoffs, residuals = TRUE, scheme = 1, pch = 19,
         xlim = c(0.5, 365.5), seWithMean = TRUE,
         trans = \(x) x + coef(m_cutoffs)['(Intercept)'])
rm(m_cutoffs)

# create the dataset
NCORES <- min(availableCores() - 4, 60)
plan(multisession, workers = NCORES)

d <-
  #' `tidyr::unnest()` without `aggregate(2)` used > TB of RAM!
  future_map2(dates$file_name, dates$date, \(.fn, .date) {
    .r <- rast(.fn, lyr = c('NDVI', 'QA')) %>%
      crop(shp, mask = TRUE)
    .month <- month(.date)
    
    # drop cloudy pixels
    .r$NDVI <- ifel(is_flagged(.r$QA, flag_position = 1), NA, .r$NDVI)
    
    # remove unrealistically high NDVI values at high latitudes
    cutoff <- doy_cutoffs$y[which(doy_cutoffs$month == .month)]
    
    if(! is.na(cutoff)) {
      .r$NDVI <- ifel(.r$NDVI > 0 & init(.r, 'y') >= cutoff, NA, .r$NDVI)
    }
    
    .r %>%
      aggregate(2, na.rm = TRUE) %>% # to ensure data frame does not exceed max size
      as.data.frame(xy = TRUE, na.rm = TRUE) %>%
      mutate(date = .date) %>%
      return()
  }, .progress = TRUE, .options = furrr_options(seed = NULL)) %>%
  bind_rows() %>%
  mutate(
    pa = extract(r_pa, select(., x, y))[, 2],
    ecodistrict = extract(r_ed, select(., x, y))[, 2],
    elev_m = extract(r_el, select(., x, y))[, 2],
    prop_water = extract(r_pw, select(., x, y))[, 2],
    year = year(date),
    doy = yday(date))

saveRDS(d, 'Data/ndvi-data.rds')

# check dataset size relative to max data frame size
nrow(d) / .Machine$integer.max # 27.3%

# for testing
library('ggplot2')
theme_set(theme_bw())

ref_lines <-
  bind_rows(tibble(x = -145:-63, y = 69.5, group = 1),
            tibble(x = -127, y = 62:90, group = 2)) %>%
  st_as_sf(coords = c('x','y')) %>%
  st_set_crs('EPSG:4326') %>%
  st_transform('ESRI:102001') %>%
  bind_cols(., st_coordinates(.))

# check banding
d_sum <- d %>%
  summarize(mean_NDVI = mean(NDVI, na.rm = TRUE),
            var_NDVI = var(NDVI, na.rm = TRUE),
            .by = c(x, y)) %>%
  rast() %>%
  `crs<-`('EPSG:4326') %>%
  project(crs(shp)) %>%
  mask(shp) %>%
  project('ESRI:102001') %>%
  as.data.frame(xy = TRUE) %>%
  as_tibble()

d_sum %>%
  select(x, y, mean_NDVI, var_NDVI) %>%
  tidyr::pivot_longer(c(mean_NDVI, var_NDVI)) %>%
  group_by(name) %>%
  mutate(value = if_else(value > quantile(value, 0.999, na.rm = TRUE),
                         quantile(value, 0.999, na.rm = TRUE), value),
         value = (value - min(value, na.rm = TRUE)),
         value = value / max(value, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot() +
  facet_wrap(~ name) +
  geom_raster(aes(x, y, fill = value)) +
  khroma::scale_fill_smoothrainbow() +
  coord_sf(datum = "ESRI:102001") +
  theme_void()

subd <- filter(d, date <= date[1] + 10)

cowplot::plot_grid(
  ggplot() +
    geom_sf(data = shp) +
    geom_raster(aes(x, y, fill = NDVI), subd) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(limits = c(-1, 1)),
  ggplot() +
    geom_sf(data = shp) +
    geom_raster(aes(x, y, fill = pa), subd) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(),
  ggplot() +
    geom_sf(data = shp) +
    geom_raster(aes(x, y, fill = ecodistrict), subd) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(),
  ggplot() +
    geom_sf(data = shp) +
    geom_raster(aes(x, y, fill = elev_m), subd) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(),
  ggplot() +
    geom_sf(data = shp) +
    geom_raster(aes(x, y, fill = prop_water), subd) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c())
