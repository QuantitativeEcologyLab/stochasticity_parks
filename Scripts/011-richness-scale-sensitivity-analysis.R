
#Load in the necessary packages
library(terra)
library(sf)
library(mgcv)
library(gratia)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scico)
library(khroma)
library(cowplot)

#-------------------------------------------------------------
# Data import and pre-processing
#-------------------------------------------------------------

#Define project CRS
canada_albers <- 'ESRI:102001' # albers conic projection for canada

#canada shapefile
canada <- rgeoboundaries::geoboundaries("Canada") %>%
  st_transform(canada_albers)

#Import the rasters describing the mean and variance in NDVI
#and then systematically coarsen the rasters
mean_ndvi <- rast("Outputs/estimated-mean.tif")
mean_ndvi_23km <- aggregate(mean_ndvi, fact = 4, na.rm=TRUE)
mean_ndvi_50km <- aggregate(mean_ndvi, fact = 9, na.rm=TRUE)
mean_ndvi_75km <- aggregate(mean_ndvi, fact = 13, na.rm=TRUE)
var_ndvi <- rast("Outputs/estimated-variance.tif")
var_ndvi_23km <- aggregate(var_ndvi, fact = 4, na.rm=TRUE)
var_ndvi_50km <- aggregate(var_ndvi, fact = 9, na.rm=TRUE)
var_ndvi_75km <- aggregate(var_ndvi, fact = 13, na.rm=TRUE)

#Import the species richness data obtained from the IUCN
#https://www.iucnredlist.org/resources/other-spatial-downloads
#and then systematically coarsen the raster
richness <- rast("Data/species_richness.tif")
richness_23km <- aggregate(richness, fact = 4, na.rm=TRUE)
richness_50km <- aggregate(richness, fact = 9, na.rm=TRUE)
richness_75km <- aggregate(richness, fact = 13, na.rm=TRUE)



#Get the richness data frames at each resolution
#while also dropping the extreme values of the variance
richness_df <- as.data.frame(richness, xy = T); names(richness_df)[3] <- "richness"
richness_df$mu_hat <- terra::extract(mean_ndvi, richness_df[,c("x","y")])[,"mu_hat"]
richness_df$var_hat <- terra::extract(var_ndvi, richness_df[,c("x","y")])[,"s2_hat"]
richness_df <- richness_df[richness_df$var_hat < 0.04,]

#25km
richness_df_25km <- as.data.frame(richness_23km, xy = T); names(richness_df_25km)[3] <- "richness"
richness_df_25km$mu_hat <- terra::extract(mean_ndvi_23km, richness_df_25km[,c("x","y")])[,"mu_hat"]
richness_df_25km$var_hat <- terra::extract(var_ndvi_23km, richness_df_25km[,c("x","y")])[,"s2_hat"]
richness_df_25km <- richness_df_25km[richness_df_25km$var_hat < 0.04,]

#50km
richness_df_50km <- as.data.frame(richness_50km, xy = T); names(richness_df_50km)[3] <- "richness"
richness_df_50km$mu_hat <- terra::extract(mean_ndvi_50km, richness_df_50km[,c("x","y")])[,"mu_hat"]
richness_df_50km$var_hat <- terra::extract(var_ndvi_50km, richness_df_50km[,c("x","y")])[,"s2_hat"]
richness_df_50km <- richness_df_50km[richness_df_50km$var_hat < 0.04,]

#75km
richness_df_75km <- as.data.frame(richness_75km, xy = T); names(richness_df_75km)[3] <- "richness"
richness_df_75km$mu_hat <- terra::extract(mean_ndvi_75km, richness_df_75km[,c("x","y")])[,"mu_hat"]
richness_df_75km$var_hat <- terra::extract(var_ndvi_75km, richness_df_75km[,c("x","y")])[,"s2_hat"]
richness_df_75km <- richness_df_75km[richness_df_75km$var_hat < 0.04,]


#-------------------------------------------------------------
# Modelling
#-------------------------------------------------------------

#Fit the gams predicting richness at the different spatial scales

fit <-
  bam(richness ~
        s(mu_hat, k = 10) + # effect of mean ndvi
        s(sqrt(var_hat), k = 10) + # effect of the variance in ndvi
        ti(mu_hat, sqrt(var_hat), k = 3), # interactive effect of mean and variance in ndvi
      family = tw(link = "log"),
      data = richness_df,
      discrete = TRUE,
      method = 'fREML')


#25km
fit_25 <-
  bam(richness ~
        s(mu_hat, k = 10) + # effect of mean ndvi
        s(sqrt(var_hat), k = 10) + # effect of the variance in ndvi
        ti(mu_hat, sqrt(var_hat), k = 3), # interactive effect of mean and variance in ndvi
      family = tw(link = "log"),
      data = richness_df_25km,
      discrete = TRUE,
      method = 'fREML')


#50km
fit_50 <-
  bam(richness ~
        s(mu_hat, k = 10) + # effect of mean ndvi
        s(sqrt(var_hat), k = 10) + # effect of the variance in ndvi
        ti(mu_hat, sqrt(var_hat), k = 3), # interactive effect of mean and variance in ndvi
      family = tw(link = "log"),
      data = richness_df_50km,
      discrete = TRUE,
      method = 'fREML')


#75km
fit_75 <-
  bam(richness ~
        s(mu_hat, k = 10) + # effect of mean ndvi
        s(sqrt(var_hat), k = 10) + # effect of the variance in ndvi
        ti(mu_hat, sqrt(var_hat), k = 3), # interactive effect of mean and variance in ndvi
      family = tw(link = "log"),
      data = richness_df_75km,
      discrete = TRUE,
      method = 'fREML')


#-------------------------------------------------------------
# Figure generation
#-------------------------------------------------------------

#------------------------------------
# Full resolution
#------------------------------------

map_5km <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = richness), richness_df) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Species richness - ~5km x 5km", palette = 'batlow') +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 12, barheight = 0.5,
                               direction = "horizontal"))



#----------
#Panel of mean NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(min(richness_df$mu_hat),
                                  max(richness_df$mu_hat),
                                  length.out = 500),
                     var_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the variance and spatial structure
mu_preds <- predict(fit, newdata = mu_hat, type = "link", se = TRUE,
                    terms = 's(mu_hat)')

#Convert structure for plotting
mu_hat$richness_hat <- exp(mu_preds$fit)
mu_hat$richness_min <- exp(mu_preds$fit - 1.96*mu_preds$se.fit)
mu_hat$richness_max <- exp(mu_preds$fit + 1.96*mu_preds$se.fit)

#Generate panel
mu_5km <- 
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
  ylab("Multiplicative change in richness") +
  xlab("Mean NDVI") +
  scale_x_continuous(expand = c(0,0.005))

#----------
#Panel of variance in NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
var_hat <- data.frame(var_hat = seq(min(richness_df$var_hat), 0.04,
                                    length.out = 1000),
                      mu_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_preds <- predict(fit, newdata = var_hat, type = "link", se = TRUE,
                     terms = 's(sqrt(var_hat))')

#Convert structure for plotting
var_hat$richness_hat <- exp(var_preds$fit)
var_hat$richness_min <- exp(var_preds$fit - 1.96*var_preds$se.fit)
var_hat$richness_max <- exp(var_preds$fit + 1.96*var_preds$se.fit)

var_5km <- 
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
  ylab("Multiplicative change in richness") +
  xlab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0.0004)) +
  annotation_logticks(sides= "b",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)

#----------
#Panel of interactive effect of mean and variance in NDVI on richness
#----------
var_mu <- expand_grid(mu_hat = seq(min(richness_df$mu_hat),
                                   max(richness_df$mu_hat),
                                   length.out = 400),
                      var_hat = seq(min(richness_df$var_hat),
                                    max(richness_df$var_hat),
                                    length.out = 400))

var_mu$too_far <- too_far(var_mu$mu_hat, var_mu$var_hat,
                          richness_df$mu_hat, richness_df$var_hat,
                          dist = 0.02)
var_mu <- var_mu[! var_mu$too_far, ]

#Predict from the fitted model excluding the spatial structure, and marginal terms
var_mu$richness_hat <- predict(fit, newdata = var_mu, type = "response")
var_mu$richness_hat <- ifelse(var_mu$richness_hat > 300, 300, var_mu$richness_hat)

var_mu_5km <-
  ggplot(var_mu[!var_mu$too_far,]) +
  geom_raster(aes(mu_hat, var_hat, fill = richness_hat)) +
  geom_contour(aes(mu_hat, var_hat, z = richness_hat), color = 'black',
               bins = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_batlow(name = 'Estimated richness', reverse = TRUE) +
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




#------------------------------------
# 25km grid
#------------------------------------

map_25km <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = richness), richness_df_23km) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Species richness - ~25km x 25km", palette = 'batlow') +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 12, barheight = 0.5,
                               direction = "horizontal"))

#----------
#Panel of mean NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(min(richness_df$mu_hat),
                                  max(richness_df$mu_hat),
                                  length.out = 500),
                     var_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the variance and spatial structure
mu_preds <- predict(fit_25, newdata = mu_hat, type = "link", se = TRUE,
                    terms = 's(mu_hat)')

#Convert structure for plotting
mu_hat$richness_hat <- exp(mu_preds$fit)
mu_hat$richness_min <- exp(mu_preds$fit - 1.96*mu_preds$se.fit)
mu_hat$richness_max <- exp(mu_preds$fit + 1.96*mu_preds$se.fit)

#Generate panel
mu_25km <- 
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
  ylab("Multiplicative change in richness") +
  xlab("Mean NDVI") +
  scale_x_continuous(expand = c(0,0.005))

#----------
#Panel of variance in NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
var_hat <- data.frame(var_hat = seq(min(richness_df$var_hat), 0.04,
                                    length.out = 1000),
                      mu_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_preds <- predict(fit_25, newdata = var_hat, type = "link", se = TRUE,
                     terms = 's(sqrt(var_hat))')

#Convert structure for plotting
var_hat$richness_hat <- exp(var_preds$fit)
var_hat$richness_min <- exp(var_preds$fit - 1.96*var_preds$se.fit)
var_hat$richness_max <- exp(var_preds$fit + 1.96*var_preds$se.fit)

var_25km <- 
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
  ylab("Multiplicative change in richness") +
  xlab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0.0004)) +
  annotation_logticks(sides= "b",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)

#----------
#Panel of interactive effect of mean and variance in NDVI on richness
#----------
var_mu <- expand_grid(mu_hat = seq(min(richness_df$mu_hat),
                                   max(richness_df$mu_hat),
                                   length.out = 400),
                      var_hat = seq(min(richness_df$var_hat),
                                    max(richness_df$var_hat),
                                    length.out = 400))

var_mu$too_far <- too_far(var_mu$mu_hat, var_mu$var_hat,
                          richness_df$mu_hat, richness_df$var_hat,
                          dist = 0.02)
var_mu <- var_mu[! var_mu$too_far, ]

#Predict from the fitted model excluding the spatial structure, and marginal terms
var_mu$richness_hat <- predict(fit_25, newdata = var_mu, type = "response")
var_mu$richness_hat <- ifelse(var_mu$richness_hat > 300, 300, var_mu$richness_hat)

var_mu_25km <-
  ggplot(var_mu[!var_mu$too_far,]) +
  geom_raster(aes(mu_hat, var_hat, fill = richness_hat)) +
  geom_contour(aes(mu_hat, var_hat, z = richness_hat), color = 'black',
               bins = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_batlow(name = 'Estimated richness', reverse = TRUE) +
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




#------------------------------------
# 50km grid
#------------------------------------

map_50km <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = richness), richness_df_50km) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Species richness - ~50km x 50km", palette = 'batlow') +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 12, barheight = 0.5,
                               direction = "horizontal"))

#----------
#Panel of mean NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(min(richness_df$mu_hat),
                                  max(richness_df$mu_hat),
                                  length.out = 500),
                     var_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the variance and spatial structure
mu_preds <- predict(fit_50, newdata = mu_hat, type = "link", se = TRUE,
                    terms = 's(mu_hat)')

#Convert structure for plotting
mu_hat$richness_hat <- exp(mu_preds$fit)
mu_hat$richness_min <- exp(mu_preds$fit - 1.96*mu_preds$se.fit)
mu_hat$richness_max <- exp(mu_preds$fit + 1.96*mu_preds$se.fit)

#Generate panel
mu_50km <- 
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
  ylab("Multiplicative change in richness") +
  xlab("Mean NDVI") +
  scale_x_continuous(expand = c(0,0.005))

#----------
#Panel of variance in NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
var_hat <- data.frame(var_hat = seq(min(richness_df$var_hat), 0.04,
                                    length.out = 1000),
                      mu_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_preds <- predict(fit_50, newdata = var_hat, type = "link", se = TRUE,
                     terms = 's(sqrt(var_hat))')

#Convert structure for plotting
var_hat$richness_hat <- exp(var_preds$fit)
var_hat$richness_min <- exp(var_preds$fit - 1.96*var_preds$se.fit)
var_hat$richness_max <- exp(var_preds$fit + 1.96*var_preds$se.fit)

var_50km <- 
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
  ylab("Multiplicative change in richness") +
  xlab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0.0004)) +
  annotation_logticks(sides= "b",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)

#----------
#Panel of interactive effect of mean and variance in NDVI on richness
#----------
var_mu <- expand_grid(mu_hat = seq(min(richness_df$mu_hat),
                                   max(richness_df$mu_hat),
                                   length.out = 400),
                      var_hat = seq(min(richness_df$var_hat),
                                    max(richness_df$var_hat),
                                    length.out = 400))

var_mu$too_far <- too_far(var_mu$mu_hat, var_mu$var_hat,
                          richness_df$mu_hat, richness_df$var_hat,
                          dist = 0.02)
var_mu <- var_mu[! var_mu$too_far, ]

#Predict from the fitted model excluding the spatial structure, and marginal terms
var_mu$richness_hat <- predict(fit_50, newdata = var_mu, type = "response")
var_mu$richness_hat <- ifelse(var_mu$richness_hat > 300, 300, var_mu$richness_hat)

var_mu_50km <-
  ggplot(var_mu[!var_mu$too_far,]) +
  geom_raster(aes(mu_hat, var_hat, fill = richness_hat)) +
  geom_contour(aes(mu_hat, var_hat, z = richness_hat), color = 'black',
               bins = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_batlow(name = 'Estimated richness', reverse = TRUE) +
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




#------------------------------------
# 75km grid
#------------------------------------

map_75km <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = richness), richness_df_75km) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Species richness - ~75km x 75km", palette = 'batlow') +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 12, barheight = 0.5,
                               direction = "horizontal"))

#----------
#Panel of mean NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(min(richness_df$mu_hat),
                                  max(richness_df$mu_hat),
                                  length.out = 500),
                     var_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the variance and spatial structure
mu_preds <- predict(fit_75, newdata = mu_hat, type = "link", se = TRUE,
                    terms = 's(mu_hat)')

#Convert structure for plotting
mu_hat$richness_hat <- exp(mu_preds$fit)
mu_hat$richness_min <- exp(mu_preds$fit - 1.96*mu_preds$se.fit)
mu_hat$richness_max <- exp(mu_preds$fit + 1.96*mu_preds$se.fit)

#Generate panel
mu_75km <- 
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
  ylab("Multiplicative change in richness") +
  xlab("Mean NDVI") +
  scale_x_continuous(expand = c(0,0.005))

#----------
#Panel of variance in NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
var_hat <- data.frame(var_hat = seq(min(richness_df$var_hat), 0.04,
                                    length.out = 1000),
                      mu_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_preds <- predict(fit_75, newdata = var_hat, type = "link", se = TRUE,
                     terms = 's(sqrt(var_hat))')

#Convert structure for plotting
var_hat$richness_hat <- exp(var_preds$fit)
var_hat$richness_min <- exp(var_preds$fit - 1.96*var_preds$se.fit)
var_hat$richness_max <- exp(var_preds$fit + 1.96*var_preds$se.fit)

var_75km <- 
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
  ylab("Multiplicative change in richness") +
  xlab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0.0004)) +
  annotation_logticks(sides= "b",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)

#----------
#Panel of interactive effect of mean and variance in NDVI on richness
#----------
var_mu <- expand_grid(mu_hat = seq(min(richness_df$mu_hat),
                                   max(richness_df$mu_hat),
                                   length.out = 400),
                      var_hat = seq(min(richness_df$var_hat),
                                    max(richness_df$var_hat),
                                    length.out = 400))

var_mu$too_far <- too_far(var_mu$mu_hat, var_mu$var_hat,
                          richness_df$mu_hat, richness_df$var_hat,
                          dist = 0.02)
var_mu <- var_mu[! var_mu$too_far, ]

#Predict from the fitted model excluding the spatial structure, and marginal terms
var_mu$richness_hat <- predict(fit_75, newdata = var_mu, type = "response")
var_mu$richness_hat <- ifelse(var_mu$richness_hat > 300, 300, var_mu$richness_hat)

var_mu_75km <-
  ggplot(var_mu[!var_mu$too_far,]) +
  geom_raster(aes(mu_hat, var_hat, fill = richness_hat)) +
  geom_contour(aes(mu_hat, var_hat, z = richness_hat), color = 'black',
               bins = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_batlow(name = 'Estimated richness', reverse = TRUE) +
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
s_9 <- cowplot::plot_grid(map_5km, map_25km, map_50km, map_75km,
                          mu_5km, mu_25km, mu_50km, mu_75km,
                          var_5km, var_25km, var_50km, var_75km,
                          var_mu_5km, var_mu_25km, var_mu_50km, var_mu_75km,
                          labels = 'AUTO', nrow = 4)



ggsave('Figures/figure-s9-richness-scale-sensitivity.png', s_9,
       width = 6.86, height = 6.86/1, units = 'in',
       dpi = 600, bg = 'white', scale = 2)


