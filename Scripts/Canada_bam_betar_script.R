library('dplyr')    
library('ggplot2')   
library('mgcv')      # for GAMs
library('lubridate')
library('stringi') 
library('tictoc')
library('ggpubr')
library('rgeoboundaries') # to obtain shapefiles of countries

#import dataset, transform the data, and save as an rds
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

DAYS <- unique(d$date)
RES <- list()

#create a loop to predict the new data from the model
for(i in 1:length(DAYS)){
  newd <- filter(d, date == DAYS[i])
  
  #add the predicted values back to the data frame with all the data
  preds <-
    bind_cols(newd,
              as.data.frame(predict.gam(ndvi_bam, newdata = newd, type = 'response')) %>%
                rename(mu = 'predict.gam(ndvi_bam, newdata = newd, type = "response")') %>% #mean
                mutate(mu = mu * 2 - 1)) #rescale ndvi
  
  if(!is.na(preds$NDVI[1])){
    RES[[i]] <- data.frame(day = unique(preds$date),
                           mu = mean(preds$mu))
  }
  
}

RESULTS <- do.call(rbind, RES)

r <- mutate(RESULTS,
            dec_date = decimal_date(day),
            year = year(day),
            doy = yday(day))

#plot this for each ecozone ? plot all together or separate?

ggarrange(a, c, ncol = 2)

a <- ggplot(r, aes(year, mu)) +
  geom_smooth(color = 'darkgreen') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        axis.text.x = element_text("Year"),
        axis.text.y = element_text("Mean"))

c <- ggplot(r, aes(doy, mu)) +
  geom_smooth(color = 'darkgreen') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

#create predictions from the model-----------------
#create a new data set for the predictions

latlong <- d[,c(3, 8:10)]
latlong <- distinct(latlong)

latlong$year <- 2010
latlong$doy <- 100
latlong$park <- as.factor(FALSE)

preds <- bind_cols(latlong,
                   as.data.frame(predict(ndvi_bam, newdata = latlong, type = 'response')) %>%
                     rename(mu = 'predict(ndvi_bam, newdata = latlong, type = "response")') %>% # mean parameters
                     mutate(mu = mu * 2 - 1)) # rescale to [-1, 1])


#plot the estimated mean and variance from the models onto a map 

ggplot() +
  geom_raster(data = preds, aes(long, lat, fill = mu)) +
  geom_sf(data = canada, fill = NA, color = "black") +
  scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key.size = unit(1.5, "cm"))


april10 <- rast("Canada/NDVI/NOAA_files/AVHRR-Land_v005_AVH13C1_NOAA-19_20100410_c20170406113543.nc")
c <- crop(april10, canada)
plot(c$NDVI, col = ndvi_pal_lim)
