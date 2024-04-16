library('dplyr')    
library('ggplot2')   
library('mgcv')      # for GAMs
library('lubridate')
library('stringi') 
library('tictoc')
library('ggpubr')
library('rgeoboundaries') # to obtain shapefiles of countries

#import dataset
data <- readRDS('Canada/data_eco.rds')
data.thin <- slice(data, seq(from = 1, to = 228587752, by = 5)) %>% #thin the data to a resolution of every 5 days 
  mutate(NDVI_scaled = (NDVI + 1) / 2) #scale ndvi for beta distribution

canada <- geoboundaries("Canada")

d <- slice(data.thin, seq(from = 1, to = 45717551, by = 250))

mean_ndvi <-
  bam(
    NDVI_scaled ~ 
      # fixed effects
      park + ecozone +
      # global smooths
      s(long, lat, bs = 'ds', k = 700) + #area effect
      s(year, bs = 'tp', k = 17) + #year effect
      s(doy, bs = 'cc', k = 10) + #seasonal/day of year effect
      # in/out level smooths
      s(year, park, bs = 'fs', k = 15) + #yearly trends in parks
      s(doy, park, bs = 'fs', k = 12, xt = list(bs = 'cc')) + #seasonal trends in parks
      # tensor interaction terms
      ti(year, doy, bs = c('cr', 'cc'), k = c(15, 20)) + #yearly trends over time
      ti(long, lat, year, bs = c('ds', 'cr'), d = c(2, 1), k = c(20, 10)) + 
      ti(long, lat, doy, bs = c('ds', 'cc'), d = c(2, 1), k = c(20, 15)),
    family = betar(), #beta location scale distribution for the data
    data = d,
    method = 'fREML',
    discrete = TRUE,
    knots = list(doy = c(0.5, 366.5)),
    control = gam.control(nthreads = 1, trace = TRUE))

#check results of model
plot(ndvi_bam, select = 1, scheme = 3, n2 = 100, too.far = 0.02)

saveRDS(ndvi_bam, file = 'Canada/models/ndvi_bam_1GB_march5.rds')

#predict the new data from the model------------------------------

#predict new values
p <- bind_cols(d, 
               as.data.frame(predict(mean_ndvi)))

names(p)[11] <- "mu"

p$back <- p$mu * 2 - 1 #backtransform data
p$resid <- p$NDVI - p$mu #extract residuals 

saveRDS(p, "Canada/predictions_march6.rds")

#for each unique lat long, find the variance of all residuals across time 

latlong <- unique(p[,c(3, 8:10)]) #extract only the spatial data values

residuals <- list()

for(i in 1:nrow(latlong)){
  
  point <- latlong[i,] #extract single point
  
  matches <- p[which(p$lat == point$lat & p$long == point$long),] #find all data related to single point in original dataset
  
  res <- data.frame(lat = point$lat,
                    long = point$long,
                    park = point$park,
                    ecozone = point$ecozone,
                    var = var(matches$resid, na.rm = TRUE), #take variance of all residuals for that point
                    mean_res = mean(matches$resid, na.rm = TRUE),#find mean of the residuals
                    mean = mean(matches$mu, na.rm = TRUE),
                    cv = sd(matches$mu)/mean(matches$mu, na.rm = TRUE)) #calculate coefficient of variation
  
  residuals[[i]] <- res
  
}

results <- do.call(rbind, residuals)

saveRDS(results, "Canada/results_march20.rds")
