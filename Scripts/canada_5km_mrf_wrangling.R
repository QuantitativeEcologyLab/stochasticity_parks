library('foreach') #to run loops in parallel
library('doMC') #to register cores for parallel loops
library('terra') #for extracting raster values
library('sf') #for spatial data wrangling
library('dplyr') #for data wrangling
library('fasterize') #to create rasters quickly

#import files-----------------------------------------------

canadvi <- list.files(path = 'smb://files.ok.ubc.ca/rsmarcus/GitHub/NDVI/Canada/Data_5km/',
                      pattern = '.rds', full.names = TRUE)

elevation <- rast('smb://files.ok.ubc.ca/rsmarcus/GitHub/NDVI/Canada/elevation/wc2.1_30s_elev.tif')
parks <- readRDS('smb://files.ok.ubc.ca/rsmarcus/GitHub/NDVI/Canada/PAs/parks_fixeed.rds')
ecodistricts <- sf::st_read("Canada/ecozones/ecodistrict_shp/Ecodistricts/ecodistricts.shp")
ecodistricts <- st_transform(ecodistricts, crs = "EPSG:4326")

#create rasters for ecodistricts and park data ---------------

#create raster for ecodistricts
districts <- unique(ecodistricts$DISTRICT_I)
ecodistricts$DISTRICT_I <- as.numeric(factor(ecodistricts$DISTRICT_I))
ecodistricts.rast <- rast(fasterize(ecodistricts, raster::raster(ecodistricts, resolution = 0.05), field = "DISTRICT_I"))
#writeRaster(ecodistricts.rast, "ecodistricts_raster.tif")

#merge parks into one raster
parks$true <- 1 
park.rast <- rast(fasterize(parks, raster::raster(parks, resolution = 0.05), field = "true"))
#writeRaster(park.rast, "parks_raster.tif") - crashes computer?

#crop and wrangle NDVI files ---------------------------------

registerDoMC(10)

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

#annotate data with parks, ecodistrict and elevation ---------

registerDoMC(25) #register how many cores you want to use for this

#.export = c('parkTF.rast', 'ecodistricts.rast', 'elevation.rast') - test if needed
d <- foreach(i=1:length(canadvi), combine = rbind) %dopar% {
  r <- readRDS(canadvi[i])
  
  data <- mutate(r,
                 year = year(date),
                 doy = yday(date))
  
  #convert to formats that take up less space
  data$year <- as.integer(data$year)
  data$doy <- as.integer(data$doy)
  
  coords <- data[,c(1:2)]
  
  park <- terra::extract(parkTF.rast, coords, ID = FALSE)
  pdata <- cbind(data, park)
  
  district.no <- terra::extract(ecodistricts, coords, ID = FALSE)
  ddata <- cbind(pdata, district.no)
  
  elevations <- terra::extract(elevation, coords, ID = FALSE)
  edata <- cbind(ddata, elevations)
  colnames(edata)[which(names(edata) == "wc2.1_30s_elev")] <- "elevation"
  
  data$park <- as.factor(data$park)
  data$ecodistrict <- as.integer(data$ecodistrict)
  
  saveRDS(edata, paste0("smb://files.ok.ubc.ca/rsmarcus/GitHub/NDVI/Canada/Data_annotated/5km/",
                        paste(basename(canadvi[i]),"_PDE.rds", sep = ""))) #save the dataframe
  
}

