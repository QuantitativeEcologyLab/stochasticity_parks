#This script explores the relationship between environmental stochasticity
#and species richness across Canada. It is based on the outputs of the 
#country-wide modelling in the script entitled 005-modelling.R


#-------------------------------------------------------------
# Workspace and data preparation
#-------------------------------------------------------------


#Load in the necessary packages
library(terra)
library(mgcv)
library(gratia)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scico)


#Import the rasters describing the mean and variance in NDVI
mean_ndvi <- rast("Outputs/estimated-mean.tif")
var_ndvi <- rast("Outputs/estimated-variance.tif")

#Import the 2024 species richness data obtained from the IUCN
#https://www.iucnredlist.org/resources/other-spatial-downloads
richness <- rast("Data/species_richness.tif")

#Convert the species richness raster to a data frame
richness_df <- terra::as.data.frame(richness, xy = T, row.names = F)
row.names(richness_df) <- NULL
names(richness_df)[3] <- "richness"
#get mean ndvi at coordinates where there is species richness
richness_df$mu_hat <- terra::extract(mean_ndvi, richness_df[,c("x","y")])[,"mu_hat"]
#get variance in ndvi at coordinates where there is species richness
richness_df$var_hat <- terra::extract(var_ndvi, richness_df[,c("x","y")])[,"s2_hat"]

#Drop extreme values of var_hat
hist(richness_df$var_hat)
mean(richness_df$var_hat > 0.04)
sum(richness_df$var_hat > 0.04)
richness_df <- richness_df[richness_df$var_hat < 0.04,]


#-------------------------------------------------------------
# Modelling
#-------------------------------------------------------------

#Fit a gam predicting richness from mean NDVI alone
fit_mu <-
  bam(richness ~ s(mu_hat, k = 10), # effect of mean ndvi
      family = Gamma(link = "log"),
      data = richness_df,
      discrete = TRUE,
      method = 'fREML')
draw(fit_mu, n = 250, rug = FALSE) # plot smooths
summary(fit_mu)

#Fit a gam predicting richness from both mean and variance in NDVI
fit_var <-
  bam(richness ~
        s(mu_hat, k = 10) + # effect of mean ndvi
        s(sqrt(var_hat), k = 10) + # effect of the variance in ndvi
        ti(mu_hat, sqrt(var_hat), k = 3), # interactive effect of mean and variance in ndvi
      family = Gamma(link = "log"),
      data = richness_df,
      discrete = TRUE,
      method = 'fREML')

summary(fit_var)

draw(fit_var, n = 250, rug = FALSE, nrow = 1) # plot smooths

#Does including the variance improve the fit to the data?
#'Note: the AIC values are corrected for use of REML see `?mgcv::AIC.gam`
AIC(fit_mu, fit_var) %>%
  mutate(delta_aic = AIC - AIC[2])

#likelihood ratio test
anova(fit_mu, fit_var, test = "Chisq")


#-------------------------------------------------------------
# Figure generation
#-------------------------------------------------------------

# look at the means of each predictor relative to the full distribution
layout(t(1:2))
hist(richness_df$mu_hat)
abline(v = mean(richness_df$mu_hat), col = 'red')
hist(richness_df$var_hat)
abline(v = mean(richness_df$var_hat), col = 'red')
layout(1)

#----------
#Panel A effect of mean NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(min(richness_df$mu_hat),
                                  max(richness_df$mu_hat),
                                  length.out = 500),
                     var_hat = mean(richness_df$var_hat) # 0.01129803
)

#Predict from the fitted model excluding the terms related to the variance and spatial structure
mu_preds <- predict(fit_var, newdata = mu_hat, type = "link", se = TRUE)

#Convert structure for plotting
mu_hat$richness_hat <- exp(mu_preds$fit)
mu_hat$richness_min <- exp(mu_preds$fit - 1.96*mu_preds$se.fit)
mu_hat$richness_max <- exp(mu_preds$fit + 1.96*mu_preds$se.fit)

# center the function at 1
mu_hat$richness_min <- mu_hat$richness_min / mean(mu_hat$richness_hat)
mu_hat$richness_max <- mu_hat$richness_max / mean(mu_hat$richness_hat)
mu_hat$richness_hat <- mu_hat$richness_hat / mean(mu_hat$richness_hat)

#Generate panel A
A <- 
  ggplot()+
  geom_hline(aes(yintercept = 1), col = "grey80", linetype = "dashed") +
  geom_ribbon(data = mu_hat,
              aes(mu_hat,
                  ymin = richness_min,
                  ymax = richness_max),
              fill = "darkgreen",
              alpha = 0.3) +
  geom_line(data = mu_hat,
            aes(mu_hat,
                richness_hat),
            colour = "darkgreen",
            linewidth = 1) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        axis.text.x  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  ylab("Multiplicative change in species richness") +
  xlab("Mean NDVI") +
  scale_x_continuous(expand = c(0,0.005))

#----------
#Panel B effect of mean NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
var_hat <- data.frame(var_hat = seq(min(richness_df$var_hat), 0.04,
                                    length.out = 1000),
                      mu_hat = mean(richness_df$mu_hat) # 0.2066355
)

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_preds <- predict(fit_var, newdata = var_hat, type = "link", se = TRUE)

#Convert structure for plotting
var_hat$richness_hat <- exp(var_preds$fit)
var_hat$richness_min <- exp(var_preds$fit - 1.96*var_preds$se.fit)
var_hat$richness_max <- exp(var_preds$fit + 1.96*var_preds$se.fit)

# center the function at 1
var_hat$richness_min <- var_hat$richness_min / mean(var_hat$richness_hat)
var_hat$richness_max <- var_hat$richness_max / mean(var_hat$richness_hat)
var_hat$richness_hat <- var_hat$richness_hat / mean(var_hat$richness_hat)

B <- 
  ggplot()+
  geom_hline(aes(yintercept = 1), col = "grey80", linetype = "dashed") +
  geom_ribbon(data = var_hat,
              aes(var_hat,
                  ymin = richness_min,
                  ymax = richness_max),
              fill = "dodgerblue3",
              alpha = 0.3) +
  geom_line(data = var_hat,
            aes(var_hat,
                richness_hat),
            colour = "dodgerblue3",
            linewidth = 1) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        axis.text.x  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  ylab("Multiplicative change in species richness") +
  xlab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0.0004)) +
  annotation_logticks(sides= "b",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)

#----------
#Panel C interactive effect of mean and variance in NDVI on richness
#----------
var_mu <- as.data.frame(expand_grid(mu_hat = seq(min(richness_df$mu_hat),
                                                 max(richness_df$mu_hat), length.out = 400),
                                    var_hat = seq(min(richness_df$var_hat),
                                                  max(richness_df$var_hat), length.out = 400)))

var_mu$too_far <- too_far(var_mu$mu_hat, var_mu$var_hat,
                          richness_df$mu_hat, richness_df$var_hat,
                          dist = 0.05)
var_mu <- var_mu[! var_mu$too_far, ]

#Predict from the fitted model excluding the spatial structure, and marginal terms
var_mu$richness_hat <- predict(fit_var, newdata = var_mu, type = "response")

# re-center the function near 1
hist(var_mu$richness_hat)
var_mu$richness_hat <- var_mu$richness_hat / 150
var_mu$z <- log2(var_mu$richness_hat)

fill_mean_lab <- bquote(atop(bold('Multiplicative change in richness')))


C <-
  ggplot(var_mu[!var_mu$too_far,]) +
  geom_raster(aes(mu_hat, var_hat, fill = z)) +
  geom_contour(aes(mu_hat, var_hat, z = z), color = 'black') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_scico(fill_mean_lab, palette = "vik", direction = -1,
                   limits = c(-1, 1) * max(abs(log2(var_mu$richness))),
                   labels = function(.x) round(2^.x, 2)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size=10, family = "sans", face = "bold"),
    axis.title.x = element_text(size=10, family = "sans", face = "bold"),
    axis.text.y = element_text(size=8, family = "sans"),
    axis.text.x  = element_text(size=8, family = "sans"),
    plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(size=8, family = "sans", face = "bold"),
    panel.background = element_rect(fill = "grey40"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  xlab("Mean NDVI") +
  ylab("Variance in NDVI") +
  guides(fill = guide_colorbar(ticks.colour = NA, barwidth = 10,
                               barheight = 1, direction = "horizontal"))

# final figure ----
cowplot::plot_grid(A, B, C, labels = 'AUTO', nrow = 1)

ggsave('Figures/figure-4-richness_regression.png',
       width = 6.86, height = (6.86 / 3), units = 'in',
       dpi = 600, bg = 'white', scale = 2)
