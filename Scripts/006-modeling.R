library('sf')        # for spatial data
library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs
library('lubridate') # for working with dates
library('tictoc')    # for estimating time elapsed
library('ggpubr')    # for data wrangling
library('terra')     # for rasters
source('Functions/scale-ndvi.R') # to scale NDVI to [0, 1] and back

#import dataset (imports within a few minutes)
d <- readRDS('Data/ndvi-data.rds')

# test spatial smooth only (only using first 100 rasters) ----
if(FALSE) {
  library('spdep')     # to run neighbour analysis
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

# run full beta model ----------------------------------
# model with sos smooth and k = 1e3 fits in 11 hours
# not using MRF smooth because it results in visibly discrete areas
d <- d %>%
  filter(prop_water < 0.4) %>% #' `prop_water < 0.5` still gave odd values
  select(! c(QA, date, ecodistrict)) %>% # reduce dataset size
  mutate(across(c(year, doy, pa), .fns = as.integer)) %>%
  na.omit() # drop any rows with NAs
gc()

# scale NDVI to fit in [0, 1] interval
d$NDVI_scaled <- ndvi_to_01(d$NDVI)

if(file.exists('Models/canada-mean-ndvi-aggr-2-2025-07-08-beta.rds')) {
  m_beta <- readRDS('Models/canada-mean-ndvi-aggr-2-2025-07-08-beta.rds')
} else {
  tictoc::tic()
  m_beta <-
    bam(
      NDVI_scaled ~
        # global smooths
        s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
        s(y, x, bs = 'sos', k = 1e3) + # smooth of space
        s(year, bs = 'cr', k = 12) + # year effect
        s(doy, bs = 'cc', k = 10) + # seasonal/day of year effect
        s(elev_m, bs = 'cr', k = 10) + # elevation effect
        # smooths for difference between in/out of protected areas
        s(year, by = pa, bs = 'cr', k = 12) +
        s(doy, by = pa, bs = 'cc', k = 10) +
        # tensor interaction smooths
        ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
        ti(y, x, year, bs = c('sos', 'cr'), d = c(2, 1), k = c(500, 6)) +
        ti(y, x, doy, bs = c('sos', 'cc'), d = c(2, 1), k = c(500, 6)),
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

  #' `gratia::draw()` causes R to crash
  png('Figures/full-model-terms-beta.png', width = 16, height = 12,
      units = 'in', res = 300)
  plot.gam(m_beta, rug = FALSE, pages = 1, scale = 0, too.far = 0.05,
           scheme = c(1, 5, rep(1, 5), 3, 5, 5))
  dev.off()
  
  # check spatial smooth
  plot(m_beta, select = 2, n2 = 500, scheme = 5)
  
  summary(m_beta) # output is included below
  
  # fitted and predicted values are different
  head(fitted(m_beta)) #' equivalent to `head(m_beta$fitted.values)`
  head(predict(m_beta, newdata = head(d), type = 'response', discrete = F))
  head(predict(m_beta, newdata = head(d), type = 'response', discrete = T))
}

#' *UPDATE*
# model summary, since it takes a while time to run:
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
              mu_hat = predict(m_beta, newdata = d, type = 'response',
                               discrete = TRUE) %>%
                ndvi_to_11(), # map to [-1, 1]
              e = NDVI - mu_hat)
  saveRDS(d, 'Data_annotated/ndvi-data-with-fitted-and-e.rds')
}

# calculate mean + variance spatially from observed values
if(file.exists('Data_annotated/summarized-spatial-stats.rds')) {
  est <- readRDS('Data_annotated/summarized-spatial-stats.rds')
} else {
  est <- d %>%
    summarize(mu_hat = mean(mu_hat),
              s2_hat = mean(e^2), #' using `var()` would subtract `mean(e)`
              .by = c(x, y)) %>%
    mutate(cv_hat = sqrt(s2_hat) / mu_hat) %>%
    as_tibble()
  saveRDS(est, 'Data_annotated/summarized-spatial-stats.rds')
}

ggplot(est, aes(y, s2_hat)) +
  coord_cartesian(ylim = c(0.1)) +
  geom_hex(bins = 100)

filter(est, y > 70) %>%
  summarize(s2_hat = mean(s2_hat), .by = y) %>%
  filter(s2_hat > 0.0075) %>%
  arrange(y) %>%
  pull(y)

filter(d, y > 70, doy < 30) %>%
  summarize(NDVI = mean(NDVI), .by = y) %>%
  filter(NDVI > 0.0075) %>%
  arrange(y) %>%
  pull(y)

# project the estimates to the canada albers projection
canada_albers <- 'ESRI:102001'

prov <- st_geometry(canadianmaps::PROV) %>%
  st_transform(canada_albers)

if(file.exists('Data_annotated/summarized-spatial-stats-albers.rds')) {
  est_proj <- readRDS('Data_annotated/summarized-spatial-stats-albers.rds')
} else {
  est_proj <- est %>%
    rast() %>%
    `crs<-`('EPSG:4326') %>%
    project(canada_albers) %>%
    as.data.frame(xy = TRUE) %>%
    as_tibble()
  saveRDS(est_proj, 'Data_annotated/summarized-spatial-stats-albers.rds')
}

# check limits
hist(est_proj$mu_hat) # ok
hist(est_proj$s2_hat, breaks = 30) # clamp to (0, 0.06]
abline(v = 0.06, col = 'red')
quantile(est_proj$cv_hat, p = c(0.01, 0.99))
hist(est_proj$cv_hat, breaks = 5e4, xlim = c(-10, 10)) # clamp to [-3, 3]
abline(v = c(-3, 3), col = 'red')

est_proj <- est_proj %>%
  mutate(s2_hat = if_else(s2_hat > 0.06, 0.06, s2_hat),
         cv_hat = if_else(abs(cv_hat) > 3, sign(cv_hat) * 3, cv_hat))

cowplot::plot_grid(
  ggplot() +
    geom_sf(data = prov, fill = 'grey30') +
    geom_raster(aes(x, y, fill = mu_hat), est_proj) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(expression(widehat(E(NDVI))), option = 'B'),
  ggplot() +
    geom_sf(data = prov, fill = 'grey30') +
    geom_raster(aes(x, y, fill = s2_hat), est_proj) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(expression(widehat(Var(NDVI))), option = 'D',
                         limits = c(0, NA)),
  ggplot() +
    geom_sf(data = prov, fill = 'grey30') +
    geom_raster(aes(x, y, fill = cv_hat), est_proj) +
    labs(x = NULL, y = NULL) +
    scale_fill_distiller(expression(widehat(CV(NDVI))), type = 'div',
                         palette = 5),
  nrow = 1)

#additional statistics ----------------------------------------------------
#' using `lm()` because the large sample size should be enough to justify
#' the central limit theorem. using `glm()` or `gam()` results in model
#' estimates that depend on the observed values, unlike for Gaussian models
brms:::inv_logit(0.132282 + c(0, -1, 1) * 1.96 * 0.000116) %>%
  ndvi_to_11() %>%
  round(4)

mutate(est, )
predict.gam(m_beta, newdata = ., )

#' *continue editing from here*
#find average variance and confidence intervals
m_var <- lm(s2_hat ~ 1, est)
coef(m_var)
confint.lm(m_var, level = 0.95)

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

