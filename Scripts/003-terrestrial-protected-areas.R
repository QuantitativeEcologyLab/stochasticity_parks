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

# shapefile of protected areas
pas <- read_sf('Data/protected-areas/BDCAPC_CPCAD_2024.shp') %>%
  st_transform(CRS)

table(pas$BIOME, useNA = 'always')
table(pas$TYPE_E, useNA = 'always')


ggplot(pas, aes(fill = BIOME)) +
  geom_sf() +
  scale_fill_manual(values = c('#33658A', '#654321'))

# drop marine PAs
pas <- filter(pas, BIOME != 'Marin | Marine')

# convert the shapefile to a raster of TRUE/FALSE
r_pas <- pas %>%
  vect() %>%
  rasterize(y = r_0, background = 0) %>%
  crop(prov, mask = TRUE)
plot(r_pas)

writeRaster(r_pas, 'Data/protected-areas/protected-areas-0-1.tif')
