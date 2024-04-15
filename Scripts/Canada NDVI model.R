library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs
library('lubridate')
library('ggplot2') # for figures
library('ggpubr') #for arranging multiple figures
library('stringi') # for working with strings
library('tictoc') #for time keeping
library('rgeoboundaries') #for canada shapefile
source('Functions/betals.r') #beta location scale function for gam
source('Functions/ndvi_pal.r') #ndvi color palette for maps
source('Functions/ndvi_pal_lim.r') #ndvi color palette for maps

# load in all data
data <- readRDS("Canada/data_eco.rds")
data.thin <- slice(data, seq(from = 1, to = 228587752, by = 5)) %>% #thin the data to a daily resolution of every 5 days ? 
  mutate(NDVI_scaled = (NDVI + 1) / 2)

canada <- geoboundaries("Canada")

worldvi <- unique(list.files(path = 'Canada/NDVI/NOAA_Files/', 
                             pattern = ".nc", full.names = T))
parks <- sf::st_read("Canada/PAs/ProtectedConservedArea.gdb")

colors <- c('lightcyan', 'lightskyblue2', 'steelblue1', 'dodgerblue3', 'royalblue4', 'midnightblue')

#thin the data set to be more manageable
d <- slice(data.thin, seq(from = 1, to = 228587752, by = 1000))

#125 kb - 13 hrs - failed
#1.2 mb - 21 hrs - failed to converge
#5.9mb - 

#run the model
tic()
model_ndvi <-
  gam(list(
    # mean predictor
    NDVI_scaled ~
      # fixed effects
      park + ecozone +
      # global smooths
      s(long, lat, bs = 'ds', k = 50) +
      s(year, bs = 'tp', k = 15) +
      s(doy, bs = 'cc', k = 10) +
      # in/out level smooths
      s(year, park, bs = 'fs', k = 12) + #view trends over time
      s(doy, park, bs = 'fs', k = 12, xt = c(bs = 'cc')) + #view yearly trends
      # tensor interaction terms
      ti(year, doy, bs = c('cr', 'cc'), k = c(12, 10)) + #view how yearly trends changed over time
      ti(long, lat, year, bs = c('ds', 'cr'), d = c(2, 1), k = c(25, 10)) +
      ti(long, lat, doy, bs = c('ds', 'cc'), d = c(2, 1), k = c(25, 10)),
    # scale predictor
    ~
      # fixed effects
      park + ecozone +
      # global smooths
      s(long, lat, bs = 'ds', k = 50) +
      s(year, bs = 'tp', k = 15) +
      s(doy, bs = 'cc', k = 10) +
      # in/out level smooths
      s(year, park, bs = 'fs', k = 12) + #view trends over time
      s(doy, park, bs = 'fs', k = 12, xt = c(bs = 'cc')) + #view yearly trends
      # tensor interaction terms
      ti(year, doy, bs = c('cr', 'cc'), k = c(12, 10)) + #view how yearly trends changed over time
      ti(long, lat, year, bs = c('ds', 'cr'), d = c(2, 1), k = c(25, 10)) +
      ti(long, lat, doy, bs = c('ds', 'cc'), d = c(2, 1), k = c(25, 10))),
    family = betals(), #beta location scale distribution for the data
    data = slice_sample(d, n = 1e3),
    method = 'REML',
    knots = list(doy = c(0.5, 366.5)),
    control = gam.control(nthreads = 1, trace = TRUE))
toc()

#save model as rds file
saveRDS(model_ndvi, file = 'Canada/models/ndvi_betals_5mb_jan11.rds')

#plot the model and run a gam.check to check the data and k values
if(FALSE) {
  plot(model_ndvi, pages = 1, scheme = 3, scale = 0, n = 250) # plot smooths
  layout(matrix(1:4, ncol = 2))
  gam.check(model_ndvi)
  layout(1)
}

#plot smooths of predictions
#create objects to hold the unique dates and the results
DAYS <- unique(d$date)
RES <- list()

#create a loop to predict the new data from the model
for(i in 1:length(DAYS)){
  newd <- filter(d, date == DAYS[i])
  
  #add the predicted values back to the data frame with all the data
  preds <-
    bind_cols(newd,
              as.data.frame(predict.gam(model_ndvi, newdata = newd, type = 'response')) %>%
                rename(mu = V1, #mean
                       phi = V2) %>% #scale
                mutate(sigma2 = phi * (1 - mu) * mu, #variance
                       mu = mu * 2 - 1, #rescale ndvi
                       sigma2 = sigma2 * 4)) #scale variance appropriately
  
  if(!is.na(preds$NDVI[1])){
    RES[[i]] <- data.frame(day = unique(preds$date),
                           mu = mean(preds$mu),
                           var = mean(preds$sigma2))
  }
  
  print(i)
  
}

RESULTS <- do.call(rbind, RES)

r <- mutate(RESULTS,
            dec_date = decimal_date(day),
            year = year(day),
            doy = yday(day))


#plot smooths for mean and variance day and year trends

ggarrange(a, b, c, e, 
          nrow = 2,
          ncol = 2)

a <- ggplot(r, aes(year, mu)) +
  geom_smooth(color = 'darkgreen') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        axis.text.x = element_text("Year"),
        axis.text.y = element_text("Mean"))

b <- ggplot(r, aes(year, var)) +
  geom_smooth(color = 'dodgerblue3') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

c <- ggplot(r, aes(doy, mu)) +
  geom_smooth(color = 'darkgreen') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

e <- ggplot(r, aes(doy, var)) +
  geom_smooth(color = 'dodgerblue3') +
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
                   as.data.frame(predict(model_ndvi, newdata = latlong, type = 'response')) %>%
                     rename(mu = V1, # mean parameters
                            phi = V2) %>% # scale parameter
                     mutate(sigma2 = phi * (1 - mu) * mu, # calculate variance
                            mu = mu * 2 - 1, # rescale to [-1, 1]
                            sigma2 = sigma2 * 4)) 

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

ggplot() +
  geom_raster(data = preds, aes(long, lat, fill = sigma2)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn('Variance', colours = colors, limits = c(0, 1)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key.size = unit(1.5, "cm"))

#test

july19 <- rast("Canada/NDVI/NOAA_files/AVHRR-Land_v005_AVH13C1_NOAA-19_20100719_c20170406141657.nc")
j19 <- crop(july19, canada)
plot(j19$NDVI, col = ndvi_pal_lim)
