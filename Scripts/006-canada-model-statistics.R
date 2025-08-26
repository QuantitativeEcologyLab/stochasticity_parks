library('sf')      # for simple features
library('dplyr')   # for data wrangling
library('mgcv')    # for GAMs
library('tidyr')   # for data wrangling
library('purrr')   # for functional programming
library('ggplot2') # for fancy plots
library('terra')   # for rasters

options(digits = 4) # change rounding to be 4 digits

canada_albers <- 'ESRI:102001'
prov <- rgeoboundaries::geoboundaries("Canada") %>%
  st_transform(canada_albers)
est_albers <-
  readRDS('Data_annotated/summarized-spatial-stats-albers.rds') %>%
  mutate(cv_hat = if_else(cv_hat < 0, Inf, cv_hat),
         richness = extract(rast('Data/species_richness.tif'),
                            tibble(x, y))[, 2],
         extremes = extract(rast('Outputs/n-extreme-months-1981-2024.tif'),
                            tibble(x, y))[, 2])

# correlation between mean and variance
# p-value is < 2e-16 even for the first 1e3 values only
cor.test(x = est_albers$mu_hat, y = est_albers$s2_hat, method = "spearman")

#calculate quantiles
quantile(est_albers$mu_hat, probs = 0.7)
quantile(est_albers$s2_hat, probs = 0.3, na.rm = TRUE)
quantile(est_albers$cv_hat, probs = 0.3)

# find cells within the quantiles
est_albers <- est_albers %>%
  mutate(mu_hat_q = mu_hat >= quantile(mu_hat, probs = 0.7),
         s2_hat_q = s2_hat <= quantile(s2_hat, probs = 0.3),
         cv_hat_q = cv_hat <= quantile(cv_hat, probs = 0.3),
         rich_q = richness >= quantile(richness, probs = 0.7, na.rm = TRUE),
         extr_q = extremes <= quantile(extremes, probs = 0.3, na.rm = TRUE))

#statistics for quantiles ----
#average NDVI of top mean quantile
m_mu_q <- lm(mu_hat ~ 1, filter(est_albers, mu_hat_q))
coef(m_mu_q)['(Intercept)']
confint.lm(m_mu_q, level = 0.95)

#average variance of bottom variance quantile
m_s2_q <- lm(s2_hat ~ 1, filter(est_albers, s2_hat_q))
coef(m_s2_q)['(Intercept)']
confint.lm(m_s2_q, level = 0.95)

# average cv of bottom cv quantile
m_cv_q <- lm(cv_hat ~ 1, filter(est_albers, cv_hat_q))
coef(m_cv_q)['(Intercept)']
confint.lm(m_cv_q, level = 0.95)

# find how much of each "good" quantile is protected ----
# i.e. what percentage of each good quantile is in a PAs? 
# total area of canada
prov_area_km2 <- prov %>%
  st_area() %>%
  sum() %>%
  `units<-`('km^2') %>%
  as.numeric()

# 30% of the area of canada
q_area_km2 <- prov_area_km2 * mean(est_albers$mu_hat_q)

# area corresponding to each "best" quantile (equivalent across all 3)
est_albers %>%
  select(pa, mu_hat_q, s2_hat_q, cv_hat_q, rich_q, extr_q) %>%
  pivot_longer(! pa, names_to = 'parameter', values_to = 'in_q') %>%
  # find approximate proportion of quantile that is in a PA
  summarize(proportion = mean(pa, na.rm = TRUE),
            .by = c(parameter, in_q)) %>%
  mutate(percentage = proportion * 100,
         ref_area = if_else(in_q, q_area_km2, prov_area_km2 - q_area_km2),
         protected_km2 = proportion * ref_area,
         unprotected_km2 = ref_area - protected_km2,
         parameter = factor(parameter, levels = c('mu_hat_q', 's2_hat_q',
                                                  'cv_hat_q', 'rich_q',
                                                  'extr_q'))) %>%
  arrange(parameter, in_q) %>%
  filter(! is.na(in_q)) # missing some richness and extreme values
