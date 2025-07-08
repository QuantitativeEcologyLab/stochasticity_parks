library('sf')        # for spatial data
library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs
library('lubridate') # for working with dates
library('tictoc')    # for estimating time elapsed
library('ggpubr')    # for data wrangling
library('spdep')     # to run neighbour analysis
library('terra')     # for rasters
source('Functions/scale-ndvi.R') # to scale NDVI to [0, 1] and back

#import dataset + transform the data ----------------------------------
d <- readRDS('Data/ndvi-data.rds') %>%
  filter(prop_water < 0.5) %>% # drop pixels that are mostly water
  select(! c(QA, date, ecodistrict)) %>%
  mutate(across(c(year, doy, pa), .fns = as.integer))

ecodistricts <- st_read("Data/ecodistricts/Canada_Ecodistricts.shp") %>%
  st_geometry() %>%
  st_as_sf() %>%
  st_transform(crs = "EPSG:4326") %>%
  mutate(id = 1:n()) %>%
  filter(id %in% unique(d$ecodistrict))

plot(ecodistricts)
plot(rast('Data/ecodistricts/ecodistrict-id.tif'))
plot(st_geometry(ecodistricts), add = TRUE, lwd = 0.2)

d$ecodistrict <- as.factor(d$ecodistrict)

#create object to hold neighbour structure ---------------------------
#extract neighbour structure
nb <- poly2nb(ecodistricts, row.names = ecodistricts$id)
names(nb) <- attr(nb, "region.id")

# test spatial smooth only (only using first 100 rasters) ----
if(FALSE) {
  d_0 <- filter(d, date <= as.Date('1981-10-01')) %>%
    mutate(NDVI_scaled = ndvi_to_01(NDVI))
  nrow(ecodistricts) #' to find max `k` for MRF smooth
  
  # model only using MRF smooth
  m_mrf_0 <- bam(NDVI_scaled ~ s(ecodistrict, bs = 'mrf', k = 500,
                                 xt = list(nb = nb)),
                 family = gaussian(),
                 data = d_0,
                 method = 'fREML',
                 discrete = TRUE,
                 nthreads = 10,
                 control = gam.control(trace = TRUE))
  summary(m_mrf_0)
  plot_mrf(m_mrf_0, .newdata = d_0, .fun = ndvi_to_11)
  
  # model only using sos smooth
  tictoc::tic() # fits in ~1 minute
  m_sos_0 <- bam(NDVI_scaled ~ s(y, x, bs = 'sos', k = 500),
                 family = gaussian(),
                 data = d_0,
                 method = 'fREML',
                 discrete = TRUE,
                 nthreads = 10,
                 control = gam.control(trace = TRUE))
  tictoc::toc()
  summary(m_sos_0)
  gratia::draw(m_sos_0, rug = FALSE, dist = 0.01)
  
  # beta model takes much longer without appreciably different results
  tictoc::tic() # fits in ~ 800 seconds (13 minutes)
  m_sos_0_beta <- bam(NDVI_scaled ~ s(y, x, bs = 'sos', k = 500),
                      family = betar(), #beta distribution for the data
                      data = d_0,
                      method = 'fREML',
                      discrete = TRUE,
                      nthreads = 10,
                      control = gam.control(trace = TRUE))
  tictoc::toc()
  summary(m_sos_0_beta)
  gratia::draw(m_sos_0_beta, rug = FALSE, dist = 0.01)
}

hist(fitted(m_sos_0) / fitted(m_sos_0_beta), breaks = 100)

# run full model -----------------------------------------------
# not using MRF smooth because it results in visibly discrete areas
tictoc::tic()
m <-
  bam(
    NDVI ~
      # global smooths
      s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
      s(y, x, bs = 'sos', k = 1000) + # smooth of space
      s(year, bs = 'cr', k = 10) + # year effect
      s(doy, bs = 'cc', k = 10) + #seasonal/day of year effect
      s(elev_m, bs = 'cr', k = 10) + #elevation effect
      # smooths for difference between in/out of protected areas
      s(year, by = pa, bs = 'cr', k = 10) +
      s(doy, by = pa, bs = 'cc', k = 10) +
      # tensor interaction smooths
      ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
      ti(y, x, year, bs = c('sos', 'cr'), d = c(2, 1), k = c(500, 5)) +
      ti(y, x, doy, bs = c('sos', 'cc'), d = c(2, 1), k = c(500, 5)),
    family = gaussian(),
    data = slice_sample(d, n = 1e5),
    method = 'fREML',
    discrete = TRUE,
    knots = list(doy = c(0.5, 366.5)),
    nthreads = 50,
    control = gam.control(trace = TRUE))
tictoc::toc()

png('Figures/temp.png', width= 8, height = 8, units = 'in', res = 300)
plot(m, pages = 1, scheme = c(1, 5, rep(1, 5)), scale = 0, too.far = 0.01)
dev.off()

plot(m, select = 1, scheme = 3, n2 = 100, too.far = 0.02)

saveRDS(m, paste0('Models/canada-mean-ndvi-aggr-2-', Sys.Date(), '.rds'))

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

#create table of stats for an ecozone
eco.table <- function(ecozone = 1){
  
  tibble(
    ecozone = ecozone,
    #mean statistics
    mean = mean(VAR$mean[VAR$layer == ecozone]),
    upper.mean.confint = confint.lm(lm(mean ~ 1, VAR[VAR$layer == ecozone,]), level = 0.95)[2],
    lower.mean.confint = confint.lm(lm(mean ~ 1, VAR[VAR$layer == ecozone,]), level = 0.95)[1],
    mean.parks = mean(na.omit(VAR$mean[c(VAR$layer == ecozone & VAR$park == 1)])),
    upper.mean.parks.confint = confint.lm(lm(mean ~ 1, VAR[c(VAR$layer == ecozone & VAR$park == 1),]), level = 0.95)[2],
    lower.mean.parks.confint = confint.lm(lm(mean ~ 1, VAR[c(VAR$layer == ecozone & VAR$park == 1),]), level = 0.95)[1],
    mean.out = mean(na.omit(VAR$mean[c(VAR$layer == ecozone & VAR$park == 0)])),
    upper.mean.out.confint = confint.lm(lm(mean ~ 1, VAR[c(VAR$layer == ecozone & VAR$park == 0),]), level = 0.95)[2],
    lower.mean.out.confint = confint.lm(lm(mean ~ 1, VAR[c(VAR$layer == ecozone & VAR$park == 0),]), level = 0.95)[1],
    #variance statistics
    var = mean(VAR$var[VAR$layer == ecozone]),
    upper.var.confint = confint.lm(lm(var ~ 1, VAR[VAR$layer == ecozone,]), level = 0.95)[2],
    lower.var.confint = confint.lm(lm(var ~ 1, VAR[VAR$layer == ecozone,]), level = 0.95)[1],
    var.parks = mean(na.omit(VAR$var[c(VAR$layer == ecozone & VAR$park == 1)])),
    upper.var.parks.confint = confint.lm(lm(var ~ 1, VAR[c(VAR$layer == ecozone & VAR$park == 1),]), level = 0.95)[2],
    lower.var.parks.confint =  confint.lm(lm(var ~ 1, VAR[c(VAR$layer == ecozone & VAR$park == 1),]), level = 0.95)[1],
    var.out = mean(na.omit(VAR$var[c(VAR$layer == ecozone & VAR$park == 0)])),
    upper.var.out.confint = confint.lm(lm(var ~ 1, VAR[c(VAR$layer == ecozone & VAR$park == 0),]), level = 0.95)[2],
    lower.var.out.confint = confint.lm(lm(var ~ 1, VAR[c(VAR$layer == ecozone & VAR$park == 0),]), level = 0.95)[1],
    
  )
  
}

eco.table(1)

#generate tables for each ecozone & merge them all into 1 table

eco.1 <- eco.table(1)
eco.2 <- eco.table(2)
eco.3 <- eco.table(3)
eco.4 <- eco.table(4)
eco.5 <- eco.table(5)
eco.6 <- eco.table(6)
eco.7 <- eco.table(7)
eco.8 <- eco.table(8)
eco.9 <- eco.table(9)
eco.10 <- eco.table(10)
eco.11 <- eco.table(11)
eco.12 <- eco.table(12)
eco.13 <- eco.table(13)
eco.14 <- eco.table(14)
eco.15 <- eco.table(15)

eco.data <- rbind(eco.1, eco.2, eco.3, eco.4, eco.5, eco.6, eco.7, eco.8, eco.9, eco.10,
                  eco.11, eco.12, eco.13, eco.14, eco.15)

#save

write.csv(eco.data, 'Canada/ecozone.stats.csv')

