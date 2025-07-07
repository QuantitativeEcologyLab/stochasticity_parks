library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('sf')      # for shapefiles
library('terra')   # for rasters
library('ggplot2') # for fancy plots
theme_set(theme_bw() + theme(legend.position = 'top'))

# raster files are kept in a different project repo folder
# not copied over to save space and avoid corruption during copy
r_0 <- rast('../ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc')
CRS <- crs(r_0)

# shapefile of canada
prov <- st_transform(canadianmaps::PROV, CRS) %>%
  st_make_valid() %>%
  st_geometry() %>%
  st_as_sf()

# shapefile of ecodistricts
ed <- read_sf('Data/ecodistricts/Canada_Ecodistricts.shp') %>%
  st_transform(CRS) %>%
  mutate(id = 1:n()) #' `ECODISTRIC` has non-unique values    

# convert the shapefile to a raster of TRUE/FALSE
r_ed <- ed %>%
  vect() %>%
  rasterize(y = r_0, field = 'id') %>%
  crop(prov, mask = TRUE)
plot(r_ed)

writeRaster(r_ed, 'Data/ecodistricts/ecodistrict-id.tif')
