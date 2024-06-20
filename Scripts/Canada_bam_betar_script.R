library('dplyr')    
library('ggplot2')   
library('mgcv')      # for GAMs
library('lubridate')
library('stringi') 
library('tictoc')
library('ggpubr')
library('rgeoboundaries') # to obtain shapefiles of countries

#import dataset + transform the data
data <- readRDS('Canada/data_eco.rds')
data.thin <- slice(data, seq(from = 1, to = 228587752, by = 5)) %>% #thin the data to a daily resolution of every 5 days ? 
  mutate(NDVI_scaled = (NDVI + 1) / 2)

canada <- geoboundaries("Canada")

d <- slice(data.thin, seq(from = 1, to = 45717551, by = 250))

ndvi_bam <-
  bam(
    NDVI_scaled ~
      # fixed effects
      park + ecozone +
      # global smooths
      s(long, lat, bs = 'ds', k = 500) + #area effect
      s(year, bs = 'tp', k = 12) + #year effect
      s(doy, bs = 'cc', k = 10) + #seasonal/day of year effect
      # in/out level smooths
      s(year, park, bs = 'fs', k = 12) + #yearly trends in parks
      s(doy, park, bs = 'fs', k = 12, xt = list(bs = 'cc')) + #seasonal trends in parks
      s(doy, ecozone, bs = 'fs', k = 12, xt = list(bs = 'cc')) +
      s(year, ecozone, bs = 'fs', k = 12) +
      # tensor interaction terms
      ti(year, doy, bs = c('cr', 'cc'), k = c(12, 10)) + #yearly trends over time
      ti(long, lat, year, bs = c('ds', 'cr'), d = c(2, 1), k = c(25, 10)) + 
      ti(long, lat, doy, bs = c('ds', 'cc'), d = c(2, 1), k = c(25, 10)),
    family = betar(), #beta location scale distribution for the data
    data = d,
    method = 'fREML',
    discrete = TRUE,
    knots = list(doy = c(0.5, 366.5)),
    control = gam.control(nthreads = 1, trace = TRUE))

plot(ndvi_bam, select = 1, scheme = 3, n2 = 100, too.far = 0.02)

saveRDS(ndvi_bam, file = 'Canada/models/ndvi_bam_25mb_feb1.rds')


#predict the new data from the model

library('foreach')
library('doMC')
registerDoMC(10) #register how many cores you want to use for this

DAYS <- unique(d$date)

RES <- foreach(i=1:length(DAYS), .combine = rbind) %dopar% {
  newd <- filter(ddata, date == DAYS[i])
  
  #add the predicted values back to the data frame with all the data
  preds <-
    bind_cols(newd,
              as.data.frame(predict.gam(tweeds.m, newdata = newd, type = 'response')) %>%
                rename(mu = 'predict.gam(tweeds.m, newdata = newd, type = "response")') %>% #mean
                mutate(mu = mu * 2 - 1)) #rescale ndvi
  
}

r <- mutate(RESULTS,
            year = year(day),
            doy = yday(day))

#calculate mean +  variance spatially
coords <- unique(r[,1:2])

VAR <- foreach(i=1:nrow(coords), .combine = rbind) %dopar% {
  #extract all rows that match that coord pair
  point <- coords[i,]
  matches <- unique(r[r$long == point$long & r$lat == point$lat,])
  
  #calculate mean NDVi and variance of residuals
  var <- bind_cols(point,
                   mu = mean(matches$mu),
                   var = var(matches$residual))
  
}

