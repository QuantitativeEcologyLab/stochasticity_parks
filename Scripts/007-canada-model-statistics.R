library('sf')      # for simple features
library('dplyr')   # for data wrangling
library('mgcv')    # for GAMs
library('tidyr')   # for data wrangling
library('purrr')   # for functional programming
library('ggplot2') # for fancy plots

canada_albers <- 'ESRI:102001'
prov <- st_transform(st_geometry(canadianmaps::PROV), canada_albers)
est_albers <- readRDS('Data_annotated/summarized-spatial-stats-albers.rds')

#' find average variance and confidence intervals
#' using `lm()` because the large sample size should be enough to justify
#' the central limit theorem.
m_var <- lm(s2_hat ~ 1, est_albers) # using equal-area projected data
coef(m_var)
confint.lm(m_var, level = 0.95)

# correlation between mean and variance
# p-value is < 2e-16 even for the first 1e3 values only
cor.test(x = est_albers$mu_hat, y = est_albers$s2_hat, method = "spearman")

# model variance within and outside of parks
m_var_pa <- lm(s2_hat ~ pa, data = mutate(est_albers, pa = factor(pa)))
summary(m_var_pa)

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
q_area_km2 <- prov_area_km2 * mean(est_albers$mu_hat_q)

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
  expand_grid(Ecozone = unique(est_albers$ecozone),
              Subset = c('Overall', 'Inside protected areas', 'Outside protected areas')) %>%
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

readr::write_csv(ecozone_summary, 'Outputs/ecozone-summary-statistics.csv')
