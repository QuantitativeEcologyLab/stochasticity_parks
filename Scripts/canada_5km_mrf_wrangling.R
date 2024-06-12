library('foreach')
library('doMC')
registerDoMC(25) #register how many cores you want to use for this

#import files
canadvi <- list.files(path = 'Canada/Data_5km/',
                      pattern = '.rds', full.names = TRUE)

#annotate data with parks, ecodistrict and elevation

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
  
  saveRDS(edata, paste0("Canada/Data_annotated/5km/",paste(basename(canadvi[i]),"_PDE.rds", sep = ""))) #save the dataframe
  
}