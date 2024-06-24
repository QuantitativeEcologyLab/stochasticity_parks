library('dplyr')    
library('ggplot2')   
library('mgcv')      # for GAMs
library('lubridate') #for data wrangling
library('tictoc')
library('ggpubr') #for data wrangling
library('spdep') #to run neighbour analysis

setwd("/run/user/1472454/gvfs/smb-share:server=files.ok.ubc.ca,share=rsmarcus/NDVI")

#import dataset + transform the data ----------------------------------
DATA <- readRDS('Canada/Data_annotated/5km_annotated_400MB.rds')
d <- mutate(DATA, NDVI_scaled = (NDVI + 1) / 2)

ecodistricts <- sf::st_read("Canada/ecozones/ecodistrict_shp/Ecodistricts/ecodistricts.shp")
ecodistricts <- st_transform(ecodistricts, crs = "EPSG:4326")

d$ecodistrict <- as.factor(d$ecodistrict)

#d <- slice(data.thin, seq(from = 1, to = 45717551, by = 250)) - subset data further if needed

#create object to hold neighbour structure ---------------------------

canada.ecod <- unique(d$ecodistrict)
ecodistricts.ca <- ecodistricts[ecodistricts$DISTRICT_I %in% canada.ecod,]

#extract neighbour structure
nb <- poly2nb(ecodistricts.ca, row.names = ecodistricts.ca$DISTRICT_I)
names(nb) <- attr(nb, "region.id")

#model data  - mrf only ------------------------------------------------

#smooth only using MRF smooth
canada.mrf <- bam(NDVI_scaled ~ s(ecodistrict, bs = 'mrf', xt = list(nb = nb)),
                  family = betar(), #beta distribution for the data
                  data = d,
                  method = 'fREML',
                  discrete = TRUE,
                  knots = list(doy = c(0.5, 366.5)),
                  control = gam.control(nthreads = 1, trace = TRUE))              

#test full mrf model -----------------------------------------------
canada.mrf.small <-
  bam(
    NDVI_scaled ~
      # fixed effects
      park + #ecozone +
      # global smooths
      s(ecodistrict, bs = 'mrf', xt = list(nb = nb)) + #specify knots?
      s(x, y, bs = 'ds', k = 300) + #area effect
      s(year, bs = 'tp', k = 10) + #year effect
      s(doy, bs = 'cc', k = 15), #seasonal/day of year effect
      # in/out level smooths
      # s(year, park, bs = 'fs', k = 12) + #yearly trends in parks
      # s(doy, park, bs = 'fs', k = 12, xt = list(bs = 'cc')) + #seasonal trends in parks
      # s(doy, ecozone, bs = 'fs', k = 12, xt = list(bs = 'cc')) +
      # s(year, ecozone, bs = 'fs', k = 12) +
      # # tensor interaction terms
      # ti(year, doy, bs = c('cr', 'cc'), k = c(12, 10)) + #yearly trends over time
      # ti(long, lat, year, bs = c('ds', 'cr'), d = c(2, 1), k = c(25, 10)) + 
      # ti(long, lat, doy, bs = c('ds', 'cc'), d = c(2, 1), k = c(25, 10)),
    family = betar(), #beta location scale distribution for the data
    data = d,
    method = 'fREML',
    discrete = TRUE,
    knots = list(doy = c(0.5, 366.5)),
    control = gam.control(nthreads = 1, trace = TRUE))

plot(canada.mrf.small, select = 1, scheme = 3, n2 = 100, too.far = 0.02)

saveRDS(canada.mrf.small, file = 'Canada/models/ndvi_simplemrf_betals_jun20.rds')

canada.mrf.small <- readRDS('Canada/models/ndvi_simplemrf_betals_jun20.rds')

#predict the new data from the model--------------------------------

library('foreach')
library('doMC')
registerDoMC(10) #register how many cores you want to use for this

tic()
  #add the predicted values back to the data frame with all the data
r <- bind_cols(d,
               as.data.frame(predict.gam(canada.mrf.small, newdata = d, type = 'response')) %>%
                rename(mu = 'predict.gam(canada.mrf.small, newdata = d, type = "response")') %>% #mean
                mutate(mu = mu * 2 - 1)) #rescale ndvi

r <- mutate(r,
            res = r$NDVI - r$mu)

#calculate mean +  variance spatially
VAR <- r %>%
  group_by(x, y) %>%
  summarize(mean = mean(mu),
            var = var(res))

#save the resulting datasets
saveRDS(r, "Canada/Data_annotated/5km_annotated_400MB_predict.rds")
saveRDS(VAR, "Canada/Data_annotated/5km_400MB_variance.rds")
