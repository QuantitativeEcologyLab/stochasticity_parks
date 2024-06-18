# this script contains all the steps taken to annotate NOAA AVHRR NDVI data, which was downloaded
# at a 5km resolution but was reprojected to approximately 50km resolution to reduce computational
# load. the updated script 'canada_5km_mrf_wrangling.R` contains the updated version of this script,
# in which the data was kept at the original 5km resolution, as well as more efficient data 
# annotation steps. 

library('tidyr')
library("rlang") 
library('lubridate')
library('dplyr') #data wrangling
library('sf') #spatial data analysis - terra dependency
library('terra') #spatial data analysis
library('rgeoboundaries') # to obtain shapefiles of countries
library('ggplot2') #for plots
library('stringi') #for data analysis involving strings
library('xml2') #for extracting links from the internet
library('rvest') #for formatting links to download from the internet
library('geodata') #for downloading elevation data
library("raster")


#import data 
canada <- geoboundaries("Canada")

worldvi <- unique(list.files(path = 'Canada/NDVI/NOAA_Files/', 
                             pattern = ".nc", full.names = T))

elevation <- raster('Canada/elevation/wc2.1_30s/wc2.1_30s_elev.tif')

parks <- sf::st_read("Canada/PAs/ProtectedConservedArea.gdb")

#remove marine protected areas and proposed protected areas
marine <- parks$TYPE_E[grepl("Marine", parks$TYPE_E)]
proposed <- parks$TYPE_E[grepl("Proposed", parks$TYPE_E)]
st_geometry(parks) <- "GEOMETRY"
terr.parks <- parks[(parks$TYPE_E %in% marine) == F,]
terr.parks <- terr.parks[(terr.parks$TYPE_E %in% proposed) == F,]

lg.parks <- terr.parks[(terr.parks$O_AREA_HA >= 5000),] #limit parks to parks 50 sq km or larger to match raster cells

lg.parks <- st_transform(lg.parks, crs = "EPSG:4326") #reproject parks to the same crs as NDVI rasters

#remove additional MPAs
lg.parks <- lg.parks[-c(885, 910),]

# check what files are corrupt ----

sizes <- file.size(worldvi) / 1e6
hist(sizes, xlab = 'Approximate file size in MB')
worldvi[sizes < 50]
plot(rast(worldvi[sizes < 50][1]))

corrupt <- sapply(worldvi,
                  function(filename) {
                    .r <- tryCatch(rast(filename),
                                   error = function(e) return(as.character(e)))
                    return(is.character(.r))
                  }) %>%
  suppressWarnings()
corrupt
corrupt <- corrupt[which(corrupt)] # only keep TRUE values
corrupt

while(any(corrupt)) {
  # find file names
  files <- substr(x = names(corrupt),
                  start = nchar('Canada/NDVI/NOAA_Files//') + 1,
                  stop = nchar(names(corrupt)))
  
  years <- substr(files,
                  start = nchar(files) - nchar('yyyymmdd_cyyyymmddhhmmss.nc') + 1,
                  stop = nchar(files) - nchar('mmdd_cyyyymmddhhmmss.nc'))
  
  # re-download the corrupt NDVi rasters
  urls <- paste0('https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access/',
                 years, '/', files)
  
  lapply(1:length(urls), function(.i){
    path <- paste0("Canada/NDVI/NOAA_Files/", files[.i])
    try(download.file(urls[.i], destfile = path))
  })
  
  # check again what files are corrupt
  corrupt <- sapply(names(corrupt),
                    function(filename) {
                      .r <- tryCatch(rast(filename),
                                     error = function(e) return(as.character(e)))
                      return(is.character(.r))
                    }) %>%
    suppressWarnings()
  corrupt <- corrupt[which(corrupt)] # only keep TRUE values
}

#crop, wrangle and annotate data ----

for(i in 1250:length(worldvi)){
  r <- rast(worldvi[i])
  c <- crop(r, canada)
  c <- project(c, "EPSG:4326")
  
  #change to 50km resolution (now approximately 5)
  c <- aggregate(c, fact = 0.5/res(c))
  
  #reproject data to canada/Albers projection
  # c <- project(c, "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=50 +lat_2=70 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs") 
  
  if(!is.null(unique(c$NDVI))) {   #if the raster has no data, skip to next file
    
    pa <- list ()
    
    for (j in 1:nrow(lg.parks)){
      
      #create data frame with combined data
      pa.dat <- 
        terra::crop(c, lg.parks[j,], snap = "out") %>%
        as.data.frame(xy = TRUE) %>%
        na.omit()
      
      #drop time of day value
      pa.dat <- pa.dat[which(names(pa.dat) != 'TIMEOFDAY')]
      
      #create a logical column to indicate data is located within a protected area
      if(nrow(pa.dat) > 0){
        pa.dat$park <- TRUE
      } 
      pa[[j]] <- pa.dat
      
    }
    
    pa <- do.call(rbind, pa)
    dat <- as.data.frame(c, xy = TRUE)
    
    if(!is_empty(pa$NDVI)) { #check there is protected area data
      
      dat <- dat[which(names(dat) != 'TIMEOFDAY')]
      dat <- na.omit(dat)
      
      #merge all data into a dataframe with ndvi outside PAs
      dat <-  base::merge(x = dat,
                          y = pa,
                          by.x = c("x","y"),
                          by.y = c("x","y"),
                          all.x = T)
      
      #data carpentry
      names(dat)[3] <- "NDVI"
      dat <- dat[,c(1:4, 7)]
      dat[is.na(dat$park),"park"] <- FALSE
      dat$date <- as.Date(substr(worldvi[i],
                                 start = nchar(worldvi[i]) - nchar('yyyymmdd_cyyyymmddhhmmss.nc') + 1,
                                 stop = nchar(worldvi[i]) - nchar('_cyyyymmddhhmmss.nc')),
                          format = '%Y%m%d')
      
    }
    
    #save the resulting dataframe
    saveRDS(dat, paste0("Canada/reproj_data/",paste(basename(worldvi[i]),"_Annotated.rds", sep = ""))) #save the dataframe
    rm(dat); rm(pa) #clean environment
    print(i) #indicate when each layer has been completed
    
  }
  
}




#debugging test

plot(c$NDVI)
plot(lg.parks, bg = "transparent", add = T)

#merge all data into one large dataframe----------------

#import files
canadvi <- list.files(path = 'Canada/reproj_data/',
                      pattern = '.rds', full.names = TRUE)

d <- list()
for(i in 1:length(canadvi)){
  r <- readRDS(canadvi[i])
  d[[i]] <- r
}

#bind data together into a single dataframe
data <- bind_rows(d, .id = "column_label")
data <- data[,c(2:7)]

#separate dates into year and day of year
data <- mutate(data,
               dec_date = decimal_date(date),
               year = year(date),
               doy = yday(date))

colnames(data) <- c("long", "lat", "NDVI", "QA", "park", "date", "dec_date", "year", "doy")

#convert to formats that take up less space
data$year <- as.integer(data$year)
data$doy <- as.integer(data$doy)
data$park <- as.factor(data$park)
data$ecozone <- as.factor(data$ecozone)

saveRDS(data, "Canada/data_full.rds")

#add in ecozones to the annotated dataset------------------------------------------------

ecozones <- sf::st_read("Canada/ecozones/ecozone_shp/Ecozones/ecozones.shp")
ecozones <- st_transform(ecozones, crs = "EPSG:4326")

data <- readRDS("Canada/data_full.rds")

data.sf <- st_as_sf(data, coords = c("long", "lat"), crs = "EPSG:4326")

missing <- ecozones[c(23,25),]

ECO <- list()

for(j in 1:nrow(ecozones)){
  ed <- data.sf[which(st_within(data.sf, missing[j,], sparse = F)),]
  ed$ecozone <- ecozones$ZONE_NAME[j]
  
  ECO[[j]] <- ed
  
  print(j)
}

dat <- do.call(rbind, ECO)

#data.eco.sf <- rbind(regs23, data.eco) #bind the final regions together for full dataset

#extract coordinates to turn back to dataframe - this step takes a long time
coords <- dat %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as.data.frame()

de <- as.data.frame(dat)
de <- cbind(de, coords)
de <- de[,-9]

colnames(data.eco)[which(names(data.eco) == "X")] <- "long"
colnames(data.eco)[which(names(data.eco) == "Y")] <- "lat"

#rename ecozones to shorter acronyms to save space
arctics <- c("Arctic Cordillera", "Northern Arctic", "Southern Arctic")
data$ecozone <- replace(data$ecozone, data$ecozone %in% arctics, "Arctic")
data$ecozone[data$ecozone == "Arctic"] <- "ARC"
data$ecozone[data$ecozone == "Taiga Cordillera"] <- "TGC"
data$ecozone[data$ecozone == "Taiga Plain"] <- "TGP"
data$ecozone[data$ecozone == "Taiga Shield"] <- "TGS"
data$ecozone[data$ecozone == "Hudson Plain"] <- "HUD"
data$ecozone[data$ecozone == "Boreal Cordillera"] <- "BRC"
data$ecozone[data$ecozone == "Boreal PLain"] <- "BRP"
data$ecozone[data$ecozone == "Boreal Shield"] <- "BRS"
data$ecozone[data$ecozone == "Prairie"] <- "PRA"
data$ecozone[data$ecozone == "Montane Cordillera"] <- "MTC"
data$ecozone[data$ecozone == "Pacific Maritime"] <- "PCM"
data$ecozone[data$ecozone == "Atlantic Maritime"] <- "ATM"
data$ecozone[data$ecozone == "MixedWood Plain"] <- "MXP"

data$ecozone <- as.factor(data$ecozone)

data <- data.eco

saveRDS(data, "Canada/data_eco.rds")

#add elevation to dataset--------------------------------------------------------

#load in final dataset with predicted values
p <- readRDS("Canada/predictions_march6.rds")

coords <- p[,c(8:9)]

elevation_world <- elevation_global(0.5, "Canada/elevation")

elevation <- crop(elevation_world, canada)

#this should work without need for anything else
elevations <- terra::extract(elevation, coords)
edata <- cbind(p, elevations)

