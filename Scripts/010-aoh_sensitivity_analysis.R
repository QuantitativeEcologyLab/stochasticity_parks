
#Load in the necessary packages
library(terra)
library(mgcv)
library(gratia)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scico)
library(khroma)
library(cowplot)
library(sf)

#-------------------------------------------------------------
# Data pre-processing
#-------------------------------------------------------------

#Define project CRS
canada_albers <- 'ESRI:102001' # albers conic projection for canada

#canada shapefile
canada <- rgeoboundaries::geoboundaries("Canada") %>%
  st_transform(canada_albers)

#Import the rasters describing the mean and variance in NDVI
mean_ndvi <- rast("Outputs/estimated-mean.tif")
var_ndvi <- rast("Outputs/estimated-variance.tif")
ecozone <- rast("Outputs/estimated-variance.tif")

#Import the 2024 species richness data obtained from the IUCN
#https://www.iucnredlist.org/resources/other-spatial-downloads
richness <- rast("Data/species_richness.tif")

#Import mammal AOH
mammals_AOH <- rast("Data/Mammals_AoH_2021_Richness.tif")
mammals_AOH <- project(mammals_AOH, richness)
mammals_AOH <- crop(mammals_AOH, richness)
mammals_AOH <- mask(mammals_AOH, richness)

#Import mammal richness
mammals_SR <- rast("Data/Mammals_SR_2025.tif")
mammals_SR <- project(mammals_SR, richness)
mammals_SR <- crop(mammals_SR, richness)
mammals_SR <- mask(mammals_SR, richness)
mammals_SR <- mask(mammals_SR, mammals_AOH)


#Get both richness datasets in a single data frame
richness_df <- as.data.frame(mammals_AOH, xy = T); names(richness_df)[3] <- "richness_aoh"
richness_df$richness <- terra::extract(mammals_SR, richness_df[,c("x","y")])[,2]


#get mean ndvi at coordinates where there is species richness
richness_df$mu_hat <- terra::extract(mean_ndvi, richness_df[,c("x","y")])[,"mu_hat"]
#get variance in ndvi at coordinates where there is species richness
richness_df$var_hat <- terra::extract(var_ndvi, richness_df[,c("x","y")])[,"s2_hat"]



# add ecozones
r_eco <- project(rast('Data/ecodistricts/ecozones.tif'), canada_albers)

richness_df <- 
  richness_df %>%
  mutate(ecozone = terra::extract(r_eco, select(., x, y))[, 2],
         ecozone = case_when(ecozone == 'Boreal PLain' ~ 'Boreal Plain',
                             ecozone == 'MixedWood Plain' ~ 'Mixedwood Plain',
                             .default = ecozone))




#-------------------------------------------------------------
# Figure S6 - Differences in the richness data
#-------------------------------------------------------------
# shapefile of ecozones
ecozones <-
  read_sf("Data/ecodistricts/Canada_Ecodistricts.shp") %>%
  st_transform(canada_albers) %>%
  summarize(geometry = st_union(geometry), .by = zone_name)

richness_df$rich_diff <-  richness_df$richness - richness_df$richness_aoh

A <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = richness), richness_df) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Mammalian species richness - IUCN range maps", palette = 'batlow') +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 18, barheight = 0.5,
                               direction = "horizontal"))

B <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = richness_aoh), richness_df) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Mammalian species richness - AOH range maps", palette = 'batlow') +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 18, barheight = 0.5,
                               direction = "horizontal"))

#spatial mean of residuals
C <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = rich_diff), richness_df) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Richness differences (IUCN - AOH)", palette = 'davos', direction = -1) +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 18, barheight = 0.5,
                               direction = "horizontal"))



boxdiff <-
  richness_df %>%
  filter(! is.na(ecozone)) %>%
  mutate(median_diff = median(rich_diff), .by = ecozone) %>%
  arrange(median_diff) %>%
  mutate(ecozone = factor(ecozone, levels = unique(ecozone))) %>%
  
  ggplot(aes(x = rich_diff, y = ecozone, group = ecozone)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.2, outlier.alpha = 0.1, fill = "#839E9C") +
  xlab("Richness differences (IUCN - AOH)") +
  ylab("Ecozone") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12, family = "sans", face = "bold"),
        axis.text = element_text(size = 10, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.7, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
        legend.position = 'none', legend.position.inside = c(0.9, 0.2))




plot_grid(
  plot_grid(A, B, C, ncol = 1),
  plot_grid(boxdiff, ncol = 1),
  ncol = 2,
  labels = 'AUTO'
)


ggsave('Figures/figure-s6-mammal-richness-differences.png',
       width = 6.86, height = 6.86/1.2, units = 'in',
       dpi = 600, bg = 'white', scale = 2)



#-------------------------------------------------------------
# Modelling
#-------------------------------------------------------------

#Fit a gam predicting standard mammal richness data
fit_var <-
  bam(richness ~
        s(mu_hat, k = 10) + # effect of mean ndvi
        s(sqrt(var_hat), k = 10) + # effect of the variance in ndvi
        ti(mu_hat, sqrt(var_hat), k = 3), # interactive effect of mean and variance in ndvi
      family = tw(link = "log"),
      data = richness_df,
      discrete = TRUE,
      method = 'fREML')


draw(fit_var, n = 250, rug = FALSE, nrow = 1) # plot smooths



#Fit a gam predicting AOH mammal richness data
fit_var_AOH <-
  bam(richness_aoh ~
        s(mu_hat, k = 10) + # effect of mean ndvi
        s(sqrt(var_hat), k = 10) + # effect of the variance in ndvi
        ti(mu_hat, sqrt(var_hat), k = 3), # interactive effect of mean and variance in ndvi
      family = tw(link = "log"),
      data = richness_df,
      discrete = TRUE,
      method = 'fREML')


draw(fit_var_AOH, n = 250, rug = FALSE, nrow = 1) # plot smooths




#-------------------------------------------------------------
# Figure generation
#-------------------------------------------------------------

#----------
#Panel A effect of mean NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(min(richness_df$mu_hat),
                                  max(richness_df$mu_hat),
                                  length.out = 500),
                     var_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the variance and spatial structure
mu_preds <- predict(fit_var, newdata = mu_hat, type = "link", se = TRUE,
                    terms = 's(mu_hat)')

#Convert structure for plotting
mu_hat$richness_hat <- exp(mu_preds$fit)
mu_hat$richness_min <- exp(mu_preds$fit - 1.96*mu_preds$se.fit)
mu_hat$richness_max <- exp(mu_preds$fit + 1.96*mu_preds$se.fit)

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
                      mu_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_preds <- predict(fit_var, newdata = var_hat, type = "link", se = TRUE,
                     terms = 's(sqrt(var_hat))')

#Convert structure for plotting
var_hat$richness_hat <- exp(var_preds$fit)
var_hat$richness_min <- exp(var_preds$fit - 1.96*var_preds$se.fit)
var_hat$richness_max <- exp(var_preds$fit + 1.96*var_preds$se.fit)

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
var_mu$richness_hat <- predict(fit_var, newdata = var_mu, type = "response")
var_mu$richness_hat <- ifelse(var_mu$richness_hat > 300, 300, var_mu$richness_hat)

C <-
  ggplot(var_mu[!var_mu$too_far,]) +
  geom_raster(aes(mu_hat, var_hat, fill = richness_hat)) +
  geom_contour(aes(mu_hat, var_hat, z = richness_hat), color = 'black',
               bins = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_batlow(name = 'Estimated species richness', reverse = TRUE) +
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






#----------
#Panel D effect of mean NDVI on AOH richness
#----------

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(min(richness_df$mu_hat),
                                  max(richness_df$mu_hat),
                                  length.out = 500),
                     var_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the variance and spatial structure
mu_preds <- predict(fit_var_AOH, newdata = mu_hat, type = "link", se = TRUE,
                    terms = 's(mu_hat)')

#Convert structure for plotting
mu_hat$richness_hat <- exp(mu_preds$fit)
mu_hat$richness_min <- exp(mu_preds$fit - 1.96*mu_preds$se.fit)
mu_hat$richness_max <- exp(mu_preds$fit + 1.96*mu_preds$se.fit)

#Generate panel D
D <- 
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
#Panel E effect of mean NDVI on AOH richness
#----------

#Build a dataframe for generating the predictions for the figure 
var_hat <- data.frame(var_hat = seq(min(richness_df$var_hat), 0.04,
                                    length.out = 1000),
                      mu_hat = 0) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_preds <- predict(fit_var_AOH, newdata = var_hat, type = "link", se = TRUE,
                     terms = 's(sqrt(var_hat))')

#Convert structure for plotting
var_hat$richness_hat <- exp(var_preds$fit)
var_hat$richness_min <- exp(var_preds$fit - 1.96*var_preds$se.fit)
var_hat$richness_max <- exp(var_preds$fit + 1.96*var_preds$se.fit)

E <- 
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
#Panel F interactive effect of mean and variance in NDVI on AOH richness
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
var_mu$richness_hat <- predict(fit_var_AOH, newdata = var_mu, type = "response")
var_mu$richness_hat <- ifelse(var_mu$richness_hat > 300, 300, var_mu$richness_hat)

F <-
  ggplot(var_mu[!var_mu$too_far,]) +
  geom_raster(aes(mu_hat, var_hat, fill = richness_hat)) +
  geom_contour(aes(mu_hat, var_hat, z = richness_hat), color = 'black',
               bins = 10) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_batlow(name = 'Estimated species richness', reverse = TRUE) +
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
s_7 <- cowplot::plot_grid(A, B, C,
                          D, E, F, labels = 'AUTO', nrow = 2)

ggsave('Figures/figure-s7-richness-regression-sensitivity-mammals.png', s_7,
       width = 6.86, height = 6.86/1.5, units = 'in',
       dpi = 600, bg = 'white', scale = 2)









#plot richness model differences (Fig S8 in appendix) --------------------------------

#Add model differences into the richness dataframe
richness_df$richness_hat <- predict(fit_var, newdata = richness_df, type = "response")
richness_df$richness_hat_AOH <- predict(fit_var_AOH, newdata = richness_df, type = "response")
richness_df$diff <- richness_df$richness_hat - richness_df$richness_hat_AOH
richness_df$diff[richness_df$diff > quantile(richness_df$diff, .99)] <- quantile(richness_df$diff, .99)



A <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = richness_hat), richness_df) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Estimated mammalian species richness - IUCN range maps", palette = 'batlow') +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 18, barheight = 0.5,
                               direction = "horizontal"))

B <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = richness_hat_AOH), richness_df) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Estimated mammalian species richness - AOH range maps", palette = 'batlow') +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 18, barheight = 0.5,
                               direction = "horizontal"))

#spatial mean of residuals
C <-
  ggplot() +
  geom_sf(data = canada, fill = 'grey', color = "transparent") + 
  geom_raster(aes(x, y, fill = diff), richness_df) +
  geom_sf(data = canada, fill = NA, color = "black", lwd = 0.2) + 
  geom_sf(data = ecozones, fill = NA, color = "black", lwd = 0.2) + 
  scale_fill_scico(name = "Difference between estimates (IUCN - AOH)", palette = 'davos', direction = -1) +
  theme_void() +
  theme(plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 18, barheight = 0.5,
                               direction = "horizontal"))



boxdiff <-
  richness_df %>%
  filter(! is.na(ecozone)) %>%
  mutate(median_diff = median(diff), .by = ecozone) %>%
  arrange(median_diff) %>%
  mutate(ecozone = factor(ecozone, levels = unique(ecozone))) %>%
  
  ggplot(aes(y = diff, x = ecozone, group = ecozone)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.2, outlier.alpha = 0.1, fill = "#839E9C") +
  ylab("Difference between estimates (IUCN - AOH)") +
  xlab("Ecozone") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.7, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
        legend.position = 'none', legend.position.inside = c(0.9, 0.2))




mu_diff_fit <- bam(diff ~
                     s(mu_hat, k = 6),
                   family = Gamma(link = "log"),
                   data = richness_df,
                   discrete = TRUE,
                   method = 'fREML')

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(0, 0.5,
                                  length.out = 1000)) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the mean and spatial structure
mu_hat$richness_hat <- predict(mu_diff_fit, newdata = mu_hat, type = "response")



mean_ndvi_diff <-
  ggplot(richness_df) +
  geom_hex(aes(x = mu_hat, y = diff)) +
  # geom_line(data = mu_hat,
  #           aes(mu_hat,
  #               richness_hat),
  #           colour = "black",
  #           linewidth = 1) +
  xlab("Mean NDVI")  +
  ylab("Difference between estimates (IUCN - AOH)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(face = 'bold'),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "top") +
  scale_fill_lapaz(name = 'Count', reverse = TRUE,
                   labels = function(x) {
                     .e <- floor(log10(x))
                     paste0(round(x / 10^.e), 'e', .e)
                   }) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 8, barheight = 0.5,
                               direction = "horizontal"))





var_diff_fit <- bam(diff ~
                      s(var_hat, k = 6),
                    family = Gamma(link = "log"),
                    data = richness_df,
                    discrete = TRUE,
                    method = 'fREML')

#Build a dataframe for generating the predictions for the figure 
var_hat <- data.frame(var_hat = seq(0.001, 0.025,
                                    length.out = 1000)) # excluded, so the value doesn't matter

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_hat$richness_hat <- predict(var_diff_fit, newdata = var_hat, type = "response")


var_ndvi_diff <-
  ggplot() +
  geom_hex(data = richness_df, aes(x = var_hat, y = diff)) +
  # geom_line(data = var_hat,
  #           aes(var_hat,
  #               richness_hat),
  #           colour = "black",
  #           linewidth = 1) +
  
  xlab("Variance in NDVI")  +
  ylab("Difference between estimates (IUCN - AOH)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9, family = "sans", face = "bold"),
        axis.text = element_text(size = 8, family = "sans"),
        legend.title = element_text(face = 'bold'),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "top") +
  scale_fill_lapaz(name = 'Count', reverse = TRUE,
                   labels = function(x) {
                     .e <- floor(log10(x))
                     paste0(round(x / 10^.e), 'e', .e)
                   }) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA,
                               barwidth = 8, barheight = 0.5,
                               direction = "horizontal")) +
  coord_cartesian(xlim = c(0,0.04))



s_8 <- plot_grid(
  plot_grid(A, B, C, nrow = 1),
  plot_grid(boxdiff, nrow = 1),
  plot_grid(mean_ndvi_diff, var_ndvi_diff, nrow = 1),
  nrow = 3,
  labels = 'AUTO'
)


ggsave('Figures/figure-s8-richness-regression-sensitivity-biases.png', s_8,
       width = 6.86, height = 6.86/.75, units = 'in',
       dpi = 600, bg = 'white', scale = 2)
