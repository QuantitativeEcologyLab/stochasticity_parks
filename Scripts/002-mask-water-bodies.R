# code modified from Stefano Mezzini's project on estimating variance in
# NDVI globally: https://github.com/QuantitativeEcologyLab/ndvi-stochasticity
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('sf')        # for shapefiles
library('terra')     # for rasters
library('elevatr')   # for digital elevation models
library('lubridate') # for working with dates
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('mgcv')      # for GAMs
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids 
library('gratia')    # for fancy plots of GAMs

# raster files are kept in a different project repo folder
# not copied over to save space and avoid corruption during copy
CRS <- crs(rast('../ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc'))

prov <- st_transform(canadianmaps::PROV, CRS) %>%
  st_make_valid() %>%
  st_geometry() %>%
  st_as_sf()

prov_bbox_wide <- st_bbox(prov) %>%
  st_as_sfc() %>%
  st_as_sf() %>%
  st_transform('EPSG:3005') %>% # BC Albers to project to a planar surface
  st_buffer(100e3) %>% # buffer by 100 km
  st_transform(CRS) # project back to lat-long

ggplot() +
  geom_sf(data = prov_bbox_wide, fill = 'red') +
  geom_sf(data = prov)

# import water features
# not dropping seas and oceans to remove problematic coast cells
linear <- read_sf('data/water-shapefiles/linear-water.shp') %>%
  st_transform(CRS) %>%
  st_geometry()
bodies <- read_sf('data/water-shapefiles/water-polygons.shp') %>%
  st_transform(CRS)
unique(warnings())

#' looks like the widest parts of rivers are included in `bodies`, while
#' `linear` only include minor rivers and canals
#' not plotting other `bodies` because some have invalid geometries
filter(bodies, Name1 == 'MacKenzie River') %>%
  st_geometry() %>%
  plot(col = 'red', border = 'red', axes = TRUE)
st_intersection(linear,
                filter(bodies, Name1 == 'MacKenzie River') %>%
                  st_bbox() %>%
                  st_as_sfc() %>%
                  st_as_sf()) %>%
  plot(col = 'cornflowerblue', add = TRUE)

# some bodies only have 2 vertices, so they give errors when plotting
if(FALSE) plot(st_geometry(filter(bodies, Shape_Leng == 0)))

head(bodies) %>%
  mutate(n_vertices = map_int(1:n(), \(i) nrow(st_coordinates(.[i, ]))))

filter(bodies, Shape_Leng == 0) %>%
  mutate(n_vertices = map_int(1:n(), \(i) nrow(st_coordinates(.[i, ]))))

# no bodies with negative length
filter(bodies, Shape_Area < 0)

# drop problematic polygons
bodies <- filter(bodies, Shape_Leng > 0)

# some polygons still have invalid geometries
filter(bodies, ISO_CC == 'CA') %>%
  slice(59) %>%
  st_coordinates() %>%
  data.frame() %>%
  ggplot() +
  geom_path(aes(X, Y))

filter(bodies, ISO_CC == 'CA') %>%
  slice(59) %>%
  st_make_valid() %>%
  st_area()

# find polygons with invalid geometries ----
# runs in ~ 21 minutes
tictoc::tic()
bodies <- bodies %>%
  mutate(valid_geom = map_lgl(1:n(), \(i) {
    ! is.na(sf:::CPL_geos_is_valid(st_geometry(.[i, ])))
  }, .progress = TRUE))
tictoc::toc()

# not many polygons with invalid geometry
sum(! bodies$valid_geom)

filter(bodies, ! valid_geom) %>%
  st_make_valid() %>%
  st_area() %>%
  as.numeric() %>%
  `/`(1e6) %>% # convert to km^2
  hist(., xlab = expression(Polygon~area~(km^2)), breaks = 50, main ='')

filter(bodies, ! valid_geom) %>%
  st_make_valid() %>%
  st_geometry() %>%
  plot()

# make invalid polygons valid
ig <- which(! bodies$valid_geom)
bodies[ig, ] <- st_make_valid(bodies[ig, ])
rm(ig)
bodies <- bodies %>%
  mutate(fixed_geom = ! valid_geom) %>%
  select(! valid_geom)

# calculate a raster of proportion of cell covered by water ----
# not masking to keep to find coastlines
r_bodies <- bodies %>%
  vect() %>%
  rasterize(y = rast(list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
                                pattern = '.nc', full.names = TRUE)[1],
                     lyr = 'QA'),
            cover = TRUE, background = 0) %>%
  crop(prov_bbox_wide) %>%
  mask(prov_bbox_wide)

plot(r_bodies)
plot(st_geometry(canadianmaps::PROV), col = 'transparent',
     border = 'white', add = TRUE)

# save and check the raster
writeRaster(r_bodies, 'Data/water-body-raster-canada.tif')
plot(rast('Data/water-body-raster-canada.tif'))
