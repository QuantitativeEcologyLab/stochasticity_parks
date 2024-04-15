library('dplyr')
library('sf')             
library('MODIStsp') #for downloading NDVI rasters
library('sp')
library('rgdal')
library('rgeos')
library('spData')
library("tidyr")
library('raster')
library("mgcv")    #for wrangling data using gam
library('ggplot2') #for visualizing data
library('rgeoboundaries')

setwd("H:/GitHub/NDVI")

map_boundary <- geoboundaries(country = 'Canada')

#import the shapefiles for the region to be observed
canada <- filter(world, name_long == "Canada")
bbox <- st_bbox(canada)

# download NDVI
MODIStsp(gui = FALSE, 
         out_folder = 'Canada/NDVI', 
         selprod = 'Vegetation Indexes_16Days_250m (M*D13Q1)', 
         prod_version = '061', 
         bandsel = 'NDVI', 
         sensor = 'Terra', 
         user = USERNAME, 
         password = PASSWORD, 
         start_date = '2005.01.01', 
         end_date = '2005.01.20', 
         spatmeth = 'bbox',
         bbox = bbox, 
         out_projsel = 'User Defined', 
         output_proj = '+proj=longlat', 
         resampling = 'bilinear', 
         delete_hdf = TRUE, 
         scale_val = TRUE, 
         ts_format = 'R RasterStack', 
         out_format = 'GTiff', 
         n_retries = 10, 
         verbose = TRUE,
         parallel = TRUE) 


#create an object to save the rasters
rasters <-
  list.files(path = 'Canada/NDVI/VI_16Days_250m_v61/NDVI/',
             pattern = '.tif', full.names = TRUE) %>%
  stack()

canada.rast.1 <- file('Canada/NDVI/VI_16Days_250m_v61/NDVI/MOD13Q1_NDVI_2005_001.tif')

saveRDS(canada.rast.1,'Canada/NDVI/canada.rasters/canada.rast.1.rds' )

# save NDVI data as an rds file of a tibble (save as a data frame)

ras <- readRDS('Canada/NDVI/canada.rasters/canada.rast.1.rds') %>%
  as.data.frame(ras, xy = TRUE) %>%
  filter(layer != 1) %>%
  rename(long = x,
         lat = y,
         ndvi = layer) %>%
  as_tibble()
