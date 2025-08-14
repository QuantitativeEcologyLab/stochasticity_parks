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
richnes_df <- terra::as.data.frame(richness, xy = T, row.names = F)
row.names(richnes_df) <- NULL
names(richnes_df)[3] <- "richness"
#get mean ndvi at coordinates where there is species richness
richnes_df$mu_hat <- terra::extract(mean_ndvi, richnes_df[,c("x","y")])[,"mu_hat"]
#get variance in ndvi at coordinates where there is species richness
richnes_df$var_hat <- terra::extract(var_ndvi, richnes_df[,c("x","y")])[,"s2_hat"]

#Drop extreme values of var_hat
#richnes_df <- richnes_df[richnes_df$var_hat < quantile(richnes_df$var_hat, .99),]


#-------------------------------------------------------------
# Modelling
#-------------------------------------------------------------

#Fit a gam predicting richness from mean NDVI alone
fit_mu <-
  bam(richness ~ # not scaling because range is in (0, 1)
        s(mu_hat, k = 5) + # effect of mean ndvi
        s(x,y, k = 80), # x,y smooth to deal with the spatial autocorrelation
      
      family = Gamma(link = "log"),
      data = richnes_df,
      discrete = TRUE,
      method = 'fREML')


#Fit a gam predicting richness from both mean and variance in NDVI
fit_var <-
  bam(richness ~ # not scaling because range is in (0, 1)
        s(mu_hat, k = 5) + # effect of mean ndvi
        s(log(var_hat), k = 5) + # effect of the variance in ndvi
        ti(mu_hat, log(var_hat), k = 6) + # interactive effect of mean and variance in ndvi
        s(x,y, k = 80), # x,y smooth to deal with the spatial autocorrelation
      
      family = Gamma(link = "log"),
      data = richnes_df,
      discrete = TRUE,
      method = 'fREML')


saveRDS(fit_var, paste0('Outputs/species_richness_gam.rds'))

summary(fit_var)

gratia::draw(fit_var, n = 250, rug = FALSE, select = "s(mu_hat)") # plot smooths


#Does including the variance improve the fit to the data?
#Note: fit with REML so AIC should be interpreted with some caution,
#but results are fairly clear
AIC(fit_mu, fit_var)

#likelihood ratio test
anova(fit_mu, fit_var, test = "Chisq")


#-------------------------------------------------------------
# Figure generation
#-------------------------------------------------------------

#----------
#Panel A effect of mean NDVI on richness
#----------

#Build a dataframe for generating the predictions for the figure 
mu_hat <- data.frame(mu_hat = seq(min(richnes_df$mu_hat),
                                  max(richnes_df$mu_hat),
                                  length.out = 500),
                     var_hat = mean(richnes_df$var_hat),
                     x = mean(richnes_df$x),
                     y = mean(richnes_df$y)
)

#Predict from the fitted model excluding the terms related to the variance and spatial structure
mu_preds <- predict(fit_var,
                    newdata = mu_hat,
                    terms = s(mu_hat),
                    #exclude = c("s(log(var_hat))","ti(mu_hat,log(var_hat))","s(x,y)", "(Intercept)"),
                    type = "link",
                    se = TRUE)

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
  ggtitle("A") +
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
  scale_x_continuous(expand = c(0,0.005)) +
  coord_cartesian(ylim = c(-1,20))



#----------
#Panel B effect of mean NDVI on richness
#----------


#Build a dataframe for generating the predictions for the figure 
var_hat <- data.frame(var_hat = seq(min(richnes_df$var_hat),
                                    0.02,
                                    length.out = 1000),
                      mu_hat = mean(richnes_df$mu_hat),
                      x = mean(richnes_df$x),
                      y = mean(richnes_df$y)
)

#Predict from the fitted model excluding the terms related to the mean and spatial structure
var_preds <- predict(fit_var, newdata = var_hat,
                     terms = s(log(var_hat)),
                     #exclude = c("s(mu_hat)", "ti(mu_hat,log(var_hat))", "s(x,y)"),
                     type = "terms",
                     se = TRUE)

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
  ggtitle("B") +
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
  # scale_x_log10(expand = c(0,0.01),
  #               breaks = c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03),
  #               labels = c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03)) +
  annotation_logticks(sides= "b",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2) +
  coord_cartesian(ylim = c(-1,20))




#----------
#Panel C interactive effect of mean and variance in NDVI on richness
#----------


var_mu <- as.data.frame(expand_grid(mu_hat = seq(min(richnes_df$mu_hat),
                                                 max(richnes_df$mu_hat), length.out = 100),
                                    var_hat = exp(seq(log(min(richnes_df$var_hat)),
                                                      log(max(richnes_df$var_hat)), length.out = 100))))
var_mu$x <- mean(richnes_df$x)
var_mu$y <- mean(richnes_df$y)


var_mu$too_far <- too_far(var_mu$mu_hat, var_mu$var_hat, richnes_df$mu_hat, richnes_df$var_hat, dist = 0.1)


#Predict from the fitted model excluding the spatial structure, and marginal terms
var_preds <- predict(fit, newdata = var_mu,
                     terms = c(ti(mu_hat,log(var_hat))),
                     #exclude = c("s(x,y)", '(Intercept)'),
                     type = "terms",
                     se = TRUE)

var_mu$richness_hat <- exp(var_preds$fit)

#Clamp multiplicative effect at 20 to match axes of A and B
var_mu[var_mu$richness_hat > 20,"richness_hat"] <- 20


fill_mean_lab <- bquote(atop(bold('Multiplicative change in richness')))


C <-
  ggplot(var_mu[!var_mu$too_far,]) +
  geom_raster(aes(mu_hat, var_hat, fill = log(richness_hat))) +
  geom_contour(aes(mu_hat, var_hat, z = log(richness_hat)), color = 'black') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0,0),
                breaks = c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03),
                labels = c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03)) +
  annotation_logticks(sides= "b",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2) +
  scale_fill_scico(fill_mean_lab,
                   palette = "vik",
                   direction = -1,
                   midpoint = 0,
                   # breaks = c(0, 2, 4, 6, 8),
                   # labels = c(0, 2, 4, 6, 8),
                   breaks = log(c(0.001, 0.01, 0.1, 1, 10)),
                   labels = c(0.001, 0.01, 0.1, 1, 10),
                   #limits = c(-2, 16)
  ) +
  ggtitle("C") +
  theme_bw() +
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
    # legend.title = element_text(size=5, family = "sans", face = "bold"),
    legend.title = element_text(size=8, family = "sans", face = "bold"),
    # legend.background = element_rect(fill = "transparent"),
    # legend.key.size = unit(0.3, 'cm'),
    # legend.spacing.y = unit(0.2, 'cm'),
    panel.background = element_rect(fill = "grey40"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  xlab("Mean NDVI") +
  ylab("Variance in NDVI") +
  guides(fill = guide_colorbar(ticks.colour = NA, barwidth = 10,
                               barheight = 1, direction = "horizontal"))


# final figure ----
# mean
cowplot::plot_grid(A, B, C,
                   nrow = 1)

ggsave('Figures/richness_regression.png',
       width = 6.86, height = (6.86 / 3), units = 'in',
       dpi = 600, bg = 'white', scale = 2)


