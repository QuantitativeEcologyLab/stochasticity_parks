library('foreach') #to run loops in parallel
library('doMC') #to register cores for parallel loops
library('terra') #for extracting raster values
library('sf') #for spatial data wrangling
library('dplyr') #for data wrangling
library('fasterize') #to create rasters quickly
library('lubridate') #for data wrangling
library('tictoc') #to run time estimates
library('purrr') #for map functions
library('furrr') #for map functions

setwd("/run/user/1472454/gvfs/smb-share:server=files.ok.ubc.ca,share=rsmarcus/NDVI")

#import files-----------------------------------------------

canadvi <- list.files(path = 'Canada/Data_5km/',
                      pattern = '.rds', full.names = TRUE)

elevation <- rast('Canada/elevation/wc2.1_30s/wc2.1_30s_elev.tif')
parks <- readRDS('Canada/PAs/parks_fixed.rds')
ecodistricts <- sf::st_read("Canada/ecozones/ecodistrict_shp/Ecodistricts/ecodistricts.shp")
ecodistricts <- st_transform(ecodistricts, crs = "EPSG:4326")

#crop and wrangle NDVI files ---------------------------------

registerDoMC(10) #register how many cores you want to use for this

foreach(i = 1:length(worldvi)) %dopar% {
  r <- rast(worldvi[i])
  r <- project(r, "EPSG:4326")
  c <- crop(r, canada)$NDVI
  c <- mask(c, canada)
  
  dat <- as.data.frame(c, xy = TRUE)
  
  if(!is.null(unique(c$NDVI))) {
    dat$date <- as.Date(substr(worldvi[i],
                               start = nchar(worldvi[i]) - nchar('yyyymmdd_cyyyymmddhhmmss.nc') + 1,
                               stop = nchar(worldvi[i]) - nchar('_cyyyymmddhhmmss.nc')),
                        format = '%Y%m%d')
    
    #save the resulting dataframe
    saveRDS(dat, paste0("Canada/Data_5km/",paste(basename(worldvi[i]),"_Annotated.rds", sep = ""))) #save the dataframe
    rm(dat); rm(r) #clean environment
    print(i) #indicate when each layer has been completed
    
  }
  
}

#create rasters for ecodistricts and park data ---------------

#create raster for ecodistricts
districts <- unique(ecodistricts$DISTRICT_I)
ecodistricts$DISTRICT_I <- as.numeric(factor(ecodistricts$DISTRICT_I))
ecodistricts.rast <- rast(fasterize(ecodistricts, raster::raster(ecodistricts, resolution = 0.05), field = "DISTRICT_I"))

#crashes R session for some reason? 
#writeRaster(ecodistricts.rast, "/run/user/1472454/gvfs/smb-share:server=files.ok.ubc.ca,share=rsmarcus/GitHub/NDVI/Canada/ecozones/ecodistricts_raster.tif")

#merge parks into one raster
parks$true <- 1 
park.rast <- rast(fasterize(parks, raster::raster(parks, resolution = 0.05), field = "true"))
#writeRaster(park.rast, "parks_raster.tif") - crashes R session

#annotate data with parks, ecodistrict and elevation ---------

registerDoMC(25) #register how many cores you want to use for this

foreach(i=1001:length(canadvi)) %dopar% {
  r <- readRDS(canadvi[i])
  
  data <- mutate(r,
                 year = year(date),
                 doy = yday(date))
  
  #convert to formats that take up less space
  data$year <- as.integer(data$year)
  data$doy <- as.integer(data$doy)
  
  coords <- data[,c(1:2)]
  
  park <- terra::extract(park.rast, coords, ID = FALSE)
  data <- cbind(data, park)
  colnames(data)[which(names(data) == "layer")] <- "park"
  data[is.na(data)] <- 0
  
  district.no <- terra::extract(ecodistricts.rast, coords, ID = FALSE)
  data <- cbind(data, district.no)
  colnames(data)[which(names(data) == "layer")] <- "ecodistrict"
  
  elevations <- terra::extract(elevation, coords, ID = FALSE)
  data <- cbind(data, elevations)
  colnames(data)[which(names(data) == "wc2.1_30s_elev")] <- "elevation"
  
  data$park <- as.factor(data$park)
  data$ecodistrict <- as.integer(data$ecodistrict)
  
  saveRDS(data, paste0("Canada/Data_annotated/5km/",
                        paste(basename(canadvi[i]),"_PDE.rds", sep = ""))) #save the dataframe
  
}

#read in all the final data files and bind them into one ----------

#because the final dataset is larger than the number of rows R can handle, 
#the data was imported in subsets, subset again, then merged together to generate
#useable subsets of the data

plan(multicore, workers = 5)

data.pde <- list.files(path = 'Canada/Data_annotated/5km/',
                      pattern = '.rds', full.names = TRUE)

#import data in subsets, then subset further

sub.data <- function(set = 1:1000, sub.by = 1000) {
  data <- furrr::future_map_dfr(data.pde[set], \(.f) readRDS(.f))
  
  slice(data, seq(from = 1, to = nrow(DATA), by = sub.by))
}

DATA1 <- sub.data(set = 1:1000, sub.by = 1000)
DATA2 <- sub.data(set = 1001:2000, sub.by = 1000)
DATA3 <- sub.data(set = 2001:3000, sub.by = 1000)
DATA4 <- sub.data(set = 3001:4000, sub.by = 1000)
DATA5 <- sub.data(set = 4001:5000, sub.by = 1000)
DATA6 <- sub.data(set = 5001:6000, sub.by = 1000)
DATA7 <- sub.data(set = 6001:7000, sub.by = 1000)
DATA8 <- sub.data(set = 7001:8000, sub.by = 1000)
DATA9 <- sub.data(set = 8001:9000, sub.by = 1000)
DATA10 <- sub.data(set = 9001:10000, sub.by = 1000)
DATA11 <- sub.data(set = 10001:11000, sub.by = 1000)
DATA12 <- sub.data(set = 11001:12000, sub.by = 1000)
DATA13 <- sub.data(set = 12001:13000, sub.by = 1000)
DATA14 <- sub.data(set = 13001:14000, sub.by = 1000)
DATA15 <- sub.data(set = 14001:15191, sub.by = 1000)

#bind all subsets together

DATA <- bind_rows(DATA1, DATA2, DATA3, DATA4, DATA5, DATA6, DATA7, DATA8,
                  DATA9, DATA10, DATA11, DATA12, DATA13, DATA14, DATA15)

saveRDS(DATA, 'Canada/Data_annotated/5km_annotated_400MB.rds')

