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

#import dataset (imports within a few minutes)
d <- readRDS('Data/ndvi-data.rds')

# test spatial smooth only (only using first 100 rasters) ----
if(FALSE) {
  d_0 <- filter(d, date <= as.Date('1981-10-01')) %>%
    mutate(NDVI_scaled = ndvi_to_01(NDVI))
  
  ecodistricts <-
    st_read("Data/ecodistricts/Canada_Ecodistricts.shp") %>%
    st_geometry() %>%
    st_as_sf() %>%
    st_transform(crs = "EPSG:4326") %>%
    mutate(id = 1:n()) %>%
    filter(id %in% unique(d_0$ecodistrict))
  
  plot(ecodistricts)
  plot(rast('Data/ecodistricts/ecodistrict-id.tif'))
  plot(st_geometry(ecodistricts), add = TRUE, lwd = 0.2)
  
  d_0$ecodistrict <- as.factor(d_0$ecodistrict)
  
  # make a list of neighbors
  nb <- poly2nb(ecodistricts, row.names = ecodistricts$id)
  names(nb) <- attr(nb, "region.id")
  
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
d <- d %>%
  filter(prop_water < 0.5) %>% # drop pixels that are mostly water
  select(! c(QA, date, ecodistrict)) %>%
  mutate(across(c(year, doy, pa), .fns = as.integer)) %>%
  na.omit()
gc()

# gaussian model (fits in 46 minutes)
if(file.exists('Models/canada-mean-ndvi-aggr-2-2025-07-08-gaussian.rds')) {
  m_gaus <- readRDS('Models/canada-mean-ndvi-aggr-2-2025-07-08-gaussian.rds')
} else {
  #' `sos` smooths need to be `s(y, x)` to avoid nonsensical plots:
  bam(NDVI ~ s(x, y, bs = 'sos', k = 100),
      family = gaussian(),
      data = slice_sample(d, n = 1e3),
      discrete = TRUE,
      method = 'fREML') %>%
    plot(phi = 0, theta = 180)
  
  bam(NDVI ~ s(y, x, bs = 'sos', k = 100),
      family = gaussian(),
      data = slice_sample(d, n = 1e3),
      discrete = TRUE,
      method = 'fREML') %>%
    plot(phi = 90, theta = 90)
  
  tictoc::tic()
  m_gaus <-
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
      data = d,
      method = 'fREML',
      discrete = TRUE,
      knots = list(doy = c(0.5, 366.5)),
      nthreads = 50,
      control = gam.control(trace = TRUE))
  tictoc::toc()
  
  saveRDS(m_gaus, paste0('Models/canada-mean-ndvi-aggr-2-', Sys.Date(),
                         '-gaussian.rds'))
  
  # gratia::draw() causes R to crash because the model is so large
  png('Figures/full-model-terms-gaussian.png', width = 16, height = 12,
      units = 'in', res = 300)
  plot.gam(m_gaus, rug = FALSE, pages = 1, scale = 0, too.far = 0.05,
           scheme = c(1, 5, rep(1, 5), 3, 5, 5))
  dev.off()
}

# model summary, since it takes some time to run:
# Family: gaussian 
# Link function: identity 
# 
# Formula:
# NDVI ~ s(prop_water, bs = "cr", k = 5) + s(y, x, bs = "sos", 
#     k = 1000) + s(year, bs = "cr", k = 10) + s(doy, bs = "cc", 
#     k = 10) + s(elev_m, bs = "cr", k = 10) + s(year, by = pa, 
#     bs = "cr", k = 10) + s(doy, by = pa, bs = "cc", k = 10) + 
#     ti(year, doy, bs = c("cr", "cc"), k = c(10, 10)) + ti(y, 
#     x, year, bs = c("sos", "cr"), d = c(2, 1), k = c(500, 5)) + 
#     ti(y, x, doy, bs = c("sos", "cc"), d = c(2, 1), k = c(500, 
#         5))
# 
# Parametric coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 6.049e-02  5.049e-05    1198   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                    edf Ref.df        F p-value    
# s(prop_water)    3.999      4   428368  <2e-16 ***
# s(y,x)         998.622    999   242863  <2e-16 ***
# s(year)          9.000      9  1036659  <2e-16 ***
# s(doy)           8.000      8 20108834  <2e-16 ***
# s(elev_m)        9.000      9  1498615  <2e-16 ***
# s(year):pa       8.997      9     4009  <2e-16 ***
# s(doy):pa        8.999      9    34315  <2e-16 ***
# ti(year,doy)    71.994     72   319028  <2e-16 ***
# ti(year,y,x)  1991.747   1996     6940  <2e-16 ***
# ti(doy,y,x)   1496.392   1497    94627  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 4613/4614
# R-sq.(adj) =  0.805   Deviance explained = 80.5%
# fREML = -4.343e+08  Scale est. = 0.012526  n = 563284823

# beta model (fits in 11 hours)
d$NDVI_scaled <- (d$NDVI + 1) / 2

if(file.exists('Models/canada-mean-ndvi-aggr-2-2025-07-08-gaussian.rds')) {
  m_beta <- readRDS('Models/canada-mean-ndvi-aggr-2-2025-07-08-gaussian.rds')
} else {
  tictoc::tic()
  m_beta <-
    bam(
      NDVI_scaled ~
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
      family = betar(link = 'logit'),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      knots = list(doy = c(0.5, 366.5)),
      nthreads = 50,
      control = gam.control(trace = TRUE))
  tictoc::toc()
  
  saveRDS(m_beta, paste0('Models/canada-mean-ndvi-aggr-2-', Sys.Date(),
                         '-beta.rds'))
  
  # gratia::draw() causes R to crash because the model is so large
  png('Figures/full-model-terms-beta.png', width = 16, height = 12,
      units = 'in', res = 300)
  plot.gam(m_beta, rug = FALSE, pages = 1, scale = 0, too.far = 0.05,
           scheme = c(1, 5, rep(1, 5), 3, 5, 5))
  dev.off()
}

# model summary, since it takes some time to run:
# Family: Beta regression(64.537) 
# Link function: logit 
# 
# Formula:
# NDVI_scaled ~ s(prop_water, bs = "cr", k = 5) + s(y, x, bs = "sos", 
#     k = 1000) + s(year, bs = "cr", k = 10) + s(doy, bs = "cc", 
#     k = 10) + s(elev_m, bs = "cr", k = 10) + s(year, by = pa, 
#     bs = "cr", k = 10) + s(doy, by = pa, bs = "cc", k = 10) + 
#     ti(year, doy, bs = c("cr", "cc"), k = c(10, 10)) + ti(y, 
#     x, year, bs = c("sos", "cr"), d = c(2, 1), k = c(500, 5)) + 
#     ti(y, x, doy, bs = c("sos", "cc"), d = c(2, 1), k = c(500, 
#         5))
# 
# Parametric coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.132282   0.000116    1140   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#                    edf Ref.df        F p-value    
# s(prop_water)    3.998      4   377788  <2e-16 ***
# s(y,x)         998.625    999   222171  <2e-16 ***
# s(year)          8.999      9   862456  <2e-16 ***
# s(doy)           7.999      8 19395980  <2e-16 ***
# s(elev_m)        9.000      9  1381785  <2e-16 ***
# s(year):pa       8.995      9     3943  <2e-16 ***
# s(doy):pa        8.998      9    24972  <2e-16 ***
# ti(year,doy)    71.991     72   354161  <2e-16 ***
# ti(year,y,x)  1992.701   1996     8631  <2e-16 ***
# ti(doy,y,x)   1496.417   1497   112841  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 4613/4614
# R-sq.(adj) =  0.809   Deviance explained = 80.4%
# fREML = -2.9496e+08  Scale est. = 1         n = 563284823

#predict the new data from the model--------------------------------
# add residuals to the plot
if(file.exists('Data_annotated/ndvi-data-with-fitted-and-e.rds')) {
  d <- readRDS('Data_annotated/ndvi-data-with-fitted-and-e.rds')
} else {
  d <- mutate(d,
              mu_hat = m_beta$fitted.values, # values on [0, 1] scale
              mu_hat = ndvi_to_11(mu_hat), # map to [-1, 1]
              e = NDVI - mu_hat)
  saveRDS(d, 'Data_annotated/ndvi-data-with-fitted-and-e.rds')
}

# calculate mean + variance spatially from observed values
if(file.exists('Canada/Data_annotated/summarized-spatial-stats.rds')) {
  est <- readRDS('Canada/Data_annotated/summarized-spatial-stats.rds')
} else {
  est <- d %>%
    summarize(mu_hat = mean(mu_hat),
              s2_hat = mean(e^2), #' using `var()` would subtract `mean(e)`
              .by = c(x, y)) %>%
    mutate(cv_hat = sqrt(s2_hat) / mu_hat)
  saveRDS(est, 'Canada/Data_annotated/summarized-spatial-stats.rds')
}

#additional statistics ----------------------------------------------------
#' using `lm()` because the large sample size should be enough to justify
#' the central limit theorem. using `glm()` or `gam()` results in model
#' estimates that depend on the observed values, unlike for Gaussian models

#' *continue editing from here*
#find average variance and confidence intervals
model <- lm(var ~ 1, est)
confint.lm(model, level = 0.95)

#correlation between mean and variance
cor.test(x = est$mean, y = est$var, method = "spearman")

#model variance within and outside of parks

var.park.gamma <- gam(var ~ park, data = est,  family = "Gamma")
plot(var.park.gamma)

#calculate quantiles

quantile(est$mean, probs = 0.7) #0.1127625 
quantile(est$var, probs = 0.3, na.rm = TRUE) #0.002437498 
quantile(est$cv, probs = 0.3) 

est$mean.quant <- with(est, ifelse(mean >= quantile(est$mean, probs = 0.7), 1, 0))
est$var.quant <- with(est, ifelse(var <= quantile(est$var, probs = 0.3, na.rm = TRUE), 1, 0))
est$cv.quant <- with(est, ifelse(cv <= quantile(est$cv, probs = 0.3), 1, 0))

#statistics for quantiles

#average NDVI of top mean quantile
mean(est$mean[est$mean.quant == 1])
model <- lm(mean ~ 1, est[est$mean.quant == 1,])
confint.lm(model, level = 0.95)

#average variance of bottom variance quantile
mean(na.omit(est$var[est$var.quant == 1]))
model <- lm(na.omit(var) ~ 1, na.omit(est[est$var.quant == 1,]))
confint.lm(model, level = 0.95)

#percentage of each quantile found in PAs

est <- na.omit(est)

#amount of high productivity land protected
est %>%
  group_by(park, mean.quant) %>%
  summarise(percent = 100 * n() / nrow(est))
#3.22% of highest productivity land is protected

(3.22*9984670)/100 #sq km of protected land in this quantile
(26.78*9984670)/100 #sq km of remaining land to protect in this quantile

#amount of low variance land protected
est %>%
  group_by(park, var.quant) %>%
  summarise(percent = 100 * n() / nrow(est))
#3.56% of lowest variance land is protected

(3.56*9984670)/100 #sq km of protected land in this quantile
(26.44*9984670)/100 #sq km of remaining land to protect in this quantile

#amount of ideal land protected
est %>%
  group_by(park, cv.quant) %>%
  summarise(percent = 100 * n() / nrow(est))
#3.81% of ideal land is protected

(3.81*9984670)/100 #sq km of protected land in this quantile
(26.19*9984670)/100 #sq km of remaining land to protect in this quantile

#mean and confidence intervals of specific ecozones

#create table of stats for an ecozone (Table S1)
eco.table <- function(ecozone = 1){
  
  tibble(
    ecozone = ecozone,
    #mean statistics
    mean = mean(est$mean[est$layer == ecozone]),
    upper.mean.confint = confint.lm(lm(mean ~ 1, est[est$layer == ecozone,]), level = 0.95)[2],
    lower.mean.confint = confint.lm(lm(mean ~ 1, est[est$layer == ecozone,]), level = 0.95)[1],
    mean.parks = mean(na.omit(est$mean[c(est$layer == ecozone & est$park == 1)])),
    upper.mean.parks.confint = confint.lm(lm(mean ~ 1, est[c(est$layer == ecozone & est$park == 1),]), level = 0.95)[2],
    lower.mean.parks.confint = confint.lm(lm(mean ~ 1, est[c(est$layer == ecozone & est$park == 1),]), level = 0.95)[1],
    mean.out = mean(na.omit(est$mean[c(est$layer == ecozone & est$park == 0)])),
    upper.mean.out.confint = confint.lm(lm(mean ~ 1, est[c(est$layer == ecozone & est$park == 0),]), level = 0.95)[2],
    lower.mean.out.confint = confint.lm(lm(mean ~ 1, est[c(est$layer == ecozone & est$park == 0),]), level = 0.95)[1],
    #variance statistics
    var = mean(est$var[est$layer == ecozone]),
    upper.var.confint = confint.lm(lm(var ~ 1, est[est$layer == ecozone,]), level = 0.95)[2],
    lower.var.confint = confint.lm(lm(var ~ 1, est[est$layer == ecozone,]), level = 0.95)[1],
    var.parks = mean(na.omit(est$var[c(est$layer == ecozone & est$park == 1)])),
    upper.var.parks.confint = confint.lm(lm(var ~ 1, est[c(est$layer == ecozone & est$park == 1),]), level = 0.95)[2],
    lower.var.parks.confint =  confint.lm(lm(var ~ 1, est[c(est$layer == ecozone & est$park == 1),]), level = 0.95)[1],
    var.out = mean(na.omit(est$var[c(est$layer == ecozone & est$park == 0)])),
    upper.var.out.confint = confint.lm(lm(var ~ 1, est[c(est$layer == ecozone & est$park == 0),]), level = 0.95)[2],
    lower.var.out.confint = confint.lm(lm(var ~ 1, est[c(est$layer == ecozone & est$park == 0),]), level = 0.95)[1],
    
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

