library('sf')        # for spatial data
library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs
library('lubridate') # for working with dates
library('tictoc')    # for estimating time elapsed
library('ggpubr')    # for data wrangling
library('tidyr')     # for data wrangling
library('terra')     # for rasters
library('purrr')     # for functional programming
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
  hist(fitted(m_sos_0) / fitted(m_sos_0_beta), breaks = 100)
}

# run full beta model ----------------------------------
# model with sos smooth and k = 1e3 fits in about a day
# not using MRF smooth because it results in visibly discrete areas
d <- d %>%
  filter(prop_water < 0.4) %>% #' `prop_water < 0.5` still gave odd values
  select(! c(QA, date, ecodistrict)) %>% # reduce dataset size
  mutate(across(c(year, doy, pa), .fns = as.integer)) %>%
  na.omit() # drop any rows with NAs
gc()

# scale NDVI to fit in [0, 1] interval
d$NDVI_scaled <- ndvi_to_01(d$NDVI)

# plot some rasters to check the data coverage after dropping cloudy pixels
d %>%
  mutate(date = paste(doy, year, sep = '-')) %>%
  filter(date %in% sample(date, 16)) %>%
  ggplot() +
  facet_wrap(~ doy + year) +
  geom_raster(aes(x, y, fill = NDVI_scaled * 1.1 - 0.1)) +
  scale_fill_gradientn(name = 'NDVI', limits = c(-0.1, 1), colours = NDVI_cols)

ggsave(filename = 'temp.png', height = 10, width = 10, bg = 'white')

if(file.exists('Models/canada-mean-ndvi-aggr-2-2025-07-20-beta.rds')) {
  m_beta <- readRDS('Models/canada-mean-ndvi-aggr-2-2025-07-20-beta.rds')
} else {
  # model with Gaussian family fits in 45 minutes
  # model with Beta family fits in 2.7 days 
  tictoc::tic()
  m_beta <-
    bam(
      NDVI_scaled ~
        # global smooths
        s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
        s(y, x, bs = 'sos', k = 500) + # smooth of space
        s(year, bs = 'cr', k = 12) + # year effect
        s(doy, bs = 'cc', k = 10) + # seasonal/day of year effect
        s(elev_m, bs = 'cr', k = 10) + # elevation effect
        # smooths for difference between in/out of protected areas
        s(year, by = pa, bs = 'cr', k = 6) +
        s(doy, by = pa, bs = 'cc', k = 6) +
        # tensor interaction smooths
        ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
        ti(y, x, year, bs = c('sos', 'cr'), d = c(2, 1), k = c(250, 6)) +
        ti(y, x, doy, bs = c('sos', 'cc'), d = c(2, 1), k = c(250, 6)),
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
  pdf('Figures/full-model-terms-beta.pdf', width = 16, height = 12)
  plot.gam(m_beta, rug = FALSE, pages = 1, scale = 0, too.far = 0.05,
           scheme = c(1, 5, rep(1, 5), 3, 5, 5))
  dev.off()
  
  # check spatial smooth alone
  plot.gam(m_beta, select = 2, rug = FALSE, scale = 0, scheme = 5,
           trans = \(x) {
             ndvi_to_11(gratia::inv_link(m_beta)(coef(m_beta)['(Intercept)'] + x))
           })
  
  summary(m_beta) # output is included below
  
  # fitted and predicted values are different
  head(fitted(m_beta)) #' equivalent to `head(m_beta$fitted.values)`
  predict(m_beta, newdata = head(d), type = 'response', discrete = F)
  predict(m_beta, newdata = head(d), type = 'response', discrete = T)
  
  # some small differences between predictions and fitted values
  head(fitted(m_beta)) -
    predict(m_beta, newdata = head(d), type = 'response', discrete = T)
}

# model summary, since it takes a while time to run:
#' `Family: Beta regression(18.374) `
#' `Link function: logit `
#' 
#' `Formula:`
#' `NDVI_scaled ~ s(prop_water, bs = "cr", k = 5) + s(y, x, bs = "sos", `
#' `    k = 500) + s(year, bs = "cr", k = 12) + s(doy, bs = "cc", `
#' `    k = 10) + s(elev_m, bs = "cr", k = 10) + s(year, by = pa, `
#' `    bs = "cr", k = 6) + s(doy, by = pa, bs = "cc", k = 6) + ti(year, `
#' `    doy, bs = c("cr", "cc"), k = c(10, 10)) + ti(y, x, year, `
#' `    bs = c("sos", "cr"), d = c(2, 1), k = c(250, 6)) + ti(y, `
#' `    x, doy, bs = c("sos", "cc"), d = c(2, 1), k = c(250, 6))`
#' 
#' `Parametric coefficients:`
#' `              Estimate Std. Error t value Pr(>|t|)    `
#' `(Intercept) -2.3473613  0.0003099   -7575   <2e-16 ***`
#' `---`
#' `Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1`
#' 
#' `Approximate significance of smooth terms:`
#' `                   edf   Ref.df        F p-value    `
#' `s(prop_water)    3.998    4.000   418014  <2e-16 ***`
#' `s(y,x)         498.945  499.000   460162  <2e-16 ***`
#' `s(year)         10.999   11.000   884807  <2e-16 ***`
#' `s(doy)           8.000    8.000 19045793  <2e-16 ***`
#' `s(elev_m)        9.000    9.000  2150842  <2e-16 ***`
#' `s(year):pa       5.993    6.000     9244  <2e-16 ***`
#' `s(doy):pa        3.893    3.995     2514  <2e-16 ***`
#' `ti(year,doy)    71.994   72.000   338227  <2e-16 ***`
#' `ti(year,y,x)  1244.316 1245.000     9428  <2e-16 ***`
#' `ti(doy,y,x)    995.696  996.000   126284  <2e-16 ***`
#' `---`
#' `Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1`
#' 
#' `Rank: 2855/2856`
#' `R-sq.(adj) =  0.816   Deviance explained = 84.7%`
#' `fREML = -7.0245e+07  Scale est. = 1         n = 554233995`

#predict the new data from the model--------------------------------
# add residuals to the plot
if(file.exists('Data_annotated/ndvi-data-with-fitted-and-e.rds')) {
  d <- readRDS('Data_annotated/ndvi-data-with-fitted-and-e.rds')
} else {
  d <- mutate(d,
              backtrasformed_fitted = ndvi_to_11(fitted(m_beta)),
              e = NDVI - backtrasformed_fitted)
  saveRDS(d, 'Data_annotated/ndvi-data-with-fitted-and-e.rds')
}

# calculate mean + variance spatially from observed values
if(file.exists('Data_annotated/summarized-spatial-stats.rds')) {
  est <- readRDS('Data_annotated/summarized-spatial-stats.rds')
} else {
  SMOOTHS <- gratia::smooths(m_beta)
  EXCLUDE <- SMOOTHS[which(grepl('year', SMOOTHS) | grepl('doy', SMOOTHS))]
  
  est <- d %>%
    summarize(unmodeled_mean = mean(NDVI, na.rm = TRUE), # for fig 1
              mean_e = mean(e, na.rm = TRUE), # for fig S1
              s2_hat = mean(e^2), #' using `var()` would subtract `mean(e)`
              prop_water = prop_water[1], # there's only one unique value
              elev_m = elev_m[1], # there's only one unique value
              pa = pa[1], # there's only one unique value
              year = 0, doy = 0, # values are irrelevant: excluded below
              .by = c(x, y)) %>%
    mutate(
      mu_hat = predict.bam(m_beta, newdata = ., type = 'response',
                           se.fit = FALSE, discrete = FALSE,
                           exclude = EXCLUDE) %>%
                           ndvi_to_11(),
      cv_hat = sqrt(s2_hat) / mu_hat) %>% # for fig 4
    as_tibble()
  saveRDS(est, 'Data_annotated/summarized-spatial-stats.rds')
  
  # test
  ggplot(est) +
    coord_sf(crs = 'EPSG:4326') +
    geom_raster(aes(x, y, fill = unmodeled_mean)) +
    scale_fill_gradientn('Mean NDVI', colours = NDVI_cols)
}

ggplot(est, aes(y, s2_hat)) +
  geom_hex(bins = 100)

# project the estimates to the canada albers projection
canada_albers <- 'ESRI:102001'

prov <- st_transform(st_geometry(canadianmaps::PROV), canada_albers)

if(file.exists('Data_annotated/summarized-spatial-stats-albers.rds')) {
  est_albers <- readRDS('Data_annotated/summarized-spatial-stats-albers.rds')
} else {
  est_albers <- est %>%
    rast() %>%
    `crs<-`('EPSG:4326') %>%
    project(canada_albers)
  
  # rasters of unmodeled mean, mean, variance, mean residuals, and CV
  plot(est_albers)
  
  writeRaster(est_albers$unmodeled_mean, 'Canada/unmodeled-mean.tif')
  writeRaster(est_albers$mu_hat, 'Canada/estimated-mean.tif')
  writeRaster(est_albers$mean_e, 'Canada/mean-residuals.tif')
  writeRaster(est_albers$s2_hat, 'Canada/mean-squared-residuals.tif')
  writeRaster(est_albers$cv_hat, 'Canada/coefficient-of-variation.tif')
  
  # add ecozones
  r_eco <- project(rast('Data/ecodistricts/ecozones.tif'), canada_albers)
  
  est_albers <- est_albers %>%
    as.data.frame(xy = TRUE) %>%
    as_tibble() %>%
    mutate(ecozone = extract(r_eco, select(., x, y))[, 2],
           ecozone = case_when(ecozone == 'Boreal PLain' ~ 'Boreal Plain',
                               ecozone == 'MixedWood Plain' ~ 'Mixedwood Plain',
                               .default = ecozone))
  mean(is.na(est_albers$ecozone)) # some NAs
  
  saveRDS(est_albers, 'Data_annotated/summarized-spatial-stats-albers.rds')
}

# check limits (CVs < 0 should really be tending to Inf)
layout(1:3)
hist(est_albers$mu_hat)
hist(est_albers$s2_hat, breaks = 30)
hist(est_albers$cv_hat, breaks = 5e4, xlim = c(-5, 5))
layout(1)

cowplot::plot_grid(
  ggplot() +
    geom_sf(data = prov, fill = 'grey30') +
    geom_raster(aes(x, y, fill = mu_hat), est_albers) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(expression(widehat(E(NDVI))), option = 'B'),
  ggplot() +
    geom_sf(data = prov, fill = 'grey30') +
    geom_raster(aes(x, y, fill = s2_hat), est_albers) +
    labs(x = NULL, y = NULL) +
    scale_fill_viridis_c(expression(widehat(Var(NDVI))), option = 'D',
                         limits = c(0, 0.05)),
  ggplot() +
    geom_sf(data = prov, fill = 'grey30') +
    geom_raster(aes(x, y, fill = cv_hat), est_albers) +
    labs(x = NULL, y = NULL) +
    scale_fill_distiller(expression(widehat(CV(NDVI))), type = 'div',
                         palette = 5, limits = c(-5, 5)),
  nrow = 1)

# additional statistics ----------------------------------------------------
# mean NDVI and CIs: -0.00400; (-0.00405, -0.00394)
# from the summary:    Intercept + critical value     * standard error
m_beta$family$linkinv(-2.3473613 + c(0, -1, 1) * 1.96 * 0.0003099) %>%
  ndvi_to_11() %>%
  round(5)

#' find average variance and confidence intervals
#' using `lm()` because the large sample size should be enough to justify
#' the central limit theorem. using `glm()` or `gam()` results in model
#' estimates that depend on the observed values, unlike for Gaussian models
m_var <- lm(s2_hat ~ 1, est_albers) # using equal-area projected data
coef(m_var)
confint.lm(m_var, level = 0.95)

# correlation between mean and variance
# p-value is < 2e-16 even for the first 1e3 values only
cor.test(x = est_albers$mu_hat, y = est_albers$s2_hat, method = "spearman")

#model variance within and outside of parks

var.pa.gamma <- gam(s2_hat ~ pa, data = est_albers, family = Gamma(link = 'log'))
summary(var.pa.gamma)

#calculate quantiles
quantile(est_albers$mu_hat, probs = 0.7)
quantile(est_albers$s2_hat, probs = 0.3, na.rm = TRUE)
quantile(est_albers$cv_hat, probs = 0.3) 

est_albers <- mutate(est_albers,
                     mu_hat_q = mu_hat >= quantile(mu_hat, probs = 0.7),
                     s2_hat_q = s2_hat <= quantile(s2_hat, probs = 0.3),
                     cv_hat_q = cv_hat <= quantile(cv_hat, probs = 0.3))

#statistics for quantiles ----
#average NDVI of top mean quantile
m_mu_q <- lm(mu_hat ~ 1, filter(est_albers, mu_hat_q))
coef(m_mu_q)['(Intercept)']
confint.lm(m_mu_q, level = 0.95)

#average variance of bottom variance quantile
m_s2_q <- lm(s2_hat ~ 1, filter(est_albers, s2_hat_q))
coef(m_s2_q)['(Intercept)']
confint.lm(m_s2_q, level = 0.95)

# find how much of each "good" quantile is protected ----
# i.e. what percentage of each good quantile is in a PAs? 
# total area of canada
prov_area_km2 <- prov %>%
  st_area() %>%
  sum() %>%
  `units<-`('km^2') %>%
  as.numeric()

# area corresponding to each "best" quantile (equivalent across all 3)
q_area_km2 <- prov_area * mean(est_albers$mu_hat_q)

est_albers %>%
  select(pa, mu_hat_q, s2_hat_q, cv_hat_q) %>%
  pivot_longer(! pa, names_to = 'parameter', values_to = 'in_q') %>%
  filter(in_q) %>% # filter to inside best quantiles only
  select(! in_q) %>% # no longer necessary
  # find approximate proportion of quantile that is in a PA
  summarize(proportion = mean(pa, na.rm = TRUE), .by = parameter) %>%
  mutate(percentage = proportion * 100,
         area_km = proportion * q_area_km2)

# mean and confidence intervals of specific ecozones (Table S1) ----
options(scipen = 10) # to prevent scientific notation in the table

# create the table of stats for an ecozone
ecozone_summary <-
  expand_grid(Ecozone = as.character(t(unique(r_eco))),
              Subset = c('Overall', 'In protected areas', 'Outside protected areas')) %>%
  mutate(output = map2(Ecozone, Subset, function(.eco, .pa) {
    # subset the data
    .d <- filter(est_albers, ecozone == .eco)
    if(.pa == 'Inside protected areas') .d <- filter(.d, pa == 1)
    if(.pa == 'Outside protected areas') .d <- filter(.d, pa == 0)
    
    # fit the models
    m_mu <- lm(mu_hat ~ 1, .d)
    m_s2 <- lm(s2_hat ~ 1, .d)
    
    make_ci <- function(x) {
      paste0('(', round(x[1], 4), ', ', round(x[2], 4), ')')
    }
    
    # extract the estimates and CIs of mean and variance
    tibble(mu_hat_est = round(coef(m_mu)['(Intercept)'], 4),
           mu_hat_95_ci = make_ci(confint.lm(m_mu, level = 0.95)),
           s2_hat_est = coef(m_s2)['(Intercept)'],
           s2_hat_95_ci = make_ci(confint.lm(m_s2, level = 0.95))) %>%
      return()
  })) %>%
  unnest(output)

print(ecozone_summary, n = 15 * 3)

readr::write_csv(ecozone_summary, 'Canada/ecozone-summary-statistics.csv')
