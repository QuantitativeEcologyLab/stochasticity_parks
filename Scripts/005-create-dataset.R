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
setwd('../ndvi-stochasticity')
source('functions/is_flagged.R') # this script sources other scripts
setwd('../stochasticity_parks')

file_names <-
  list.files('../ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/',
             pattern = '.nc', full.names = TRUE, recursive = FALSE)
length(file_names) / 365 # approximate number of years of data

# all rasters have the same CRS
# check crs for first, last, and a random sample of rasters
file_names[c(1, length(file_names),
             sample(length(file_names), size = 100))] %>%
  map_chr(function(.fn) { # fast enough that it is not worth parallelizing
    crs(rast(.fn))
  }) %>%
  unique()

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

layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(r_pa)
plot(r_ed)
plot(r_pw)
plot(r_el)
layout(1)

# some coastlines have elevations < 0 m
plot(r_el > 0)

# import a data frame of all dates
# see https://github.com/QuantitativeEcologyLab/ndvi-stochasticity for info
dates <- readRDS('../ndvi-stochasticity/data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds') %>%
  # rasters are not saved in this repo
  mutate(file_name = paste0('../ndvi-stochasticity/', file_name))

# some dates are missing a raster (not available on the server)
all(! is.na(dates$file_name))
mean(! is.na(dates$file_name))
filter(dates, is.na(file_name))

# create the aggregated datasets
# spatRast objects cannot be run in parallel and moved across sessions:
# https://stackoverflow.com/questions/67445883/terra-package-returns-error-when-try-to-run-parallel-operations/67449818#67449818
shp <- st_read('Data/ecodistricts/Canada_Ecodistricts.shp') %>%
  st_geometry() %>%
  st_as_sf() %>%
  st_transform(crs(rast(file_names[1])))

NCORES <- min(availableCores() - 4, 60)
plan(multisession, workers = NCORES)

d <-
  dates %>%
  filter(! is.na(file_name)) %>% # drop missing rasters
  mutate(
    # ndvi data
    ndvi_data = future_map2(file_name, date, \(fn, .date) {
      .r <- rast(fn, lyr = c('NDVI', 'QA')) %>%
        crop(shp, mask = TRUE)
      .month <- month(.date)
      
      # drop cloudy pixels
      .r$NDVI <- ifel(is_flagged(.r$QA, flag_position = 1), NA, .r$NDVI)
      
      # remove unrealistically high NDVI values at high latitudes
      # (are no cells above 70 N, so not filtering in sept or oct)
      if(.month %in% c(1:4, 11:12)) {
        # create a raster with TRUE if above max latitude 
        .lats <- init(.r, 'y') >= 60
        .r$NDVI <- ifel(.r$NDVI > 0.2 & .lats, NA, .r$NDVI)
      }
      
      return(as.data.frame(.r, xy = TRUE, na.rm = TRUE))
    }, .progress = TRUE, .options = furrr_options(seed = NULL))) %>%
  select(! c(file_name, n_cells)) %>%
  unnest(ndvi_data) %>%
  mutate(
    pa = extract(r_pa, select(., x, y))[, 2],
    ecodistrict = extract(r_ed, select(., x, y))[, 2],
    elev_m = extract(r_el, select(., x, y))[, 2],
    prop_water = extract(r_pw, select(., x, y))[, 2],
    year = year(date),
    doy = yday(date))

saveRDS(d, 'Data/ndvi-data.rds')

# for testing
library('ggplot2')
theme_set(theme_bw())

subd <- filter(readRDS('Data/ndvi-data.rds'), date <= date[1] + 10)

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
