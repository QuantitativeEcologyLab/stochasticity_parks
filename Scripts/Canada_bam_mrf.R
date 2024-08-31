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
      s(x, y, bs = 'ds', k = 800) + #area effect
      s(year, bs = 'tp', k = 10) + #year effect
      s(doy, bs = 'cc', k = 15) + #seasonal/day of year effect
      s(elevation, bs = "tp", k = 5), #elevation effect
      # in/out level smooths
      s(year, park, bs = "fs", k = 12) + 
      s(doy, park, bs = "fs", k = 12, xt = list(bs = "cc")) + 
      #tensor interaction smooths
      ti(year, doy, bs = c("cr", "cc"), k = c(12, 10)) + 
      ti(x, y, year, bs = c("ds", "cr"), d = c(2, 1), k = c(50, 10)) + 
      ti(x, y, doy, bs = c("ds", "cc"), d = c(2, 1), k = c(50, 10)),
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

#additional statistics ------------------------------------------------------

#find average variance and confidence intervals
model <- lm(var ~ 1, VAR)
confint.lm(model, level = 0.95)

#correlation between mean and variance
cor.test(x = VAR$mean, y = VAR$var, method = "spearman")

#model variance within and outside of parks

var.park.gamma <- gam(var ~ park, data = VAR,  family = "Gamma")
plot(var.park.gamma)

#calculate quantiles

quantile(VAR$mean, probs = 0.7) #0.1127625 
quantile(VAR$var, probs = 0.3, na.rm = TRUE) #0.002437498 
quantile(VAR$cv, probs = 0.3) 

VAR$mean.quant <- with(VAR, ifelse(mean >= quantile(VAR$mean, probs = 0.7), 1, 0))
VAR$var.quant <- with(VAR, ifelse(var <= quantile(VAR$var, probs = 0.3, na.rm = TRUE), 1, 0))
VAR$cv.quant <- with(VAR, ifelse(cv <= quantile(VAR$cv, probs = 0.3), 1, 0))

#statistics for quantiles

#average NDVI of top mean quantile
mean(VAR$mean[VAR$mean.quant == 1])
model <- lm(mean ~ 1, VAR[VAR$mean.quant == 1,])
confint.lm(model, level = 0.95)

#average variance of bottom variance quantile
mean(na.omit(VAR$var[VAR$var.quant == 1]))
model <- lm(na.omit(var) ~ 1, na.omit(VAR[VAR$var.quant == 1,]))
confint.lm(model, level = 0.95)

#percentage of each quantile found in PAs

VAR <- na.omit(VAR)

#amount of high productivity land protected
VAR %>%
  group_by(park, mean.quant) %>%
  summarise(percent = 100 * n() / nrow(VAR))
#3.22% of highest productivity land is protected

(3.22*9984670)/100 #sq km of protected land in this quantile
(26.78*9984670)/100 #sq km of remaining land to protect in this quantile

#amount of low variance land protected
VAR %>%
  group_by(park, var.quant) %>%
  summarise(percent = 100 * n() / nrow(VAR))
#3.56% of lowest variance land is protected

(3.56*9984670)/100 #sq km of protected land in this quantile
(26.44*9984670)/100 #sq km of remaining land to protect in this quantile

#amount of ideal land protected
VAR %>%
  group_by(park, cv.quant) %>%
  summarise(percent = 100 * n() / nrow(VAR))
#3.81% of ideal land is protected

(3.81*9984670)/100 #sq km of protected land in this quantile
(26.19*9984670)/100 #sq km of remaining land to protect in this quantile

#mean and confidence intervals of specific ecozones

mci.eco <- function(ecozone = 1){
  
  mean(VAR$mean[VAR$ecozone == ecozone])
  
  model <- lm(mean ~ 1, VAR[VAR$ecozone == ecozone,])
  confint.lm(model, level = 0.95)
  
  }

mci.eco(15)

