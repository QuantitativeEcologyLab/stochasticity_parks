library('dplyr')
library('ggplot2') #for fancy figures
library('ggpubr') #arrange multiple figures
library('lubridate') #for decimal_date function
library('mgcv') #for predict.gam
library('sf') #for importing ecozone data
library('terra') #for basemaps
library('rgeoboundaries') #for canada shapefile
library('basemaps') #for figure basemaps
library('tidyterra')
library('ggspatial')
library('raster') #create rasters to work with, cv function
library('ncdf4') #to write ncdf files
library('fasterize')

setwd("/run/user/1472454/gvfs/smb-share:server=files.ok.ubc.ca,share=rsmarcus/NDVI")

#import data -----------------------------------------------

#complete dataset with predicted values for each observation
r <- readRDS("Canada/Data_annotated/DATA_preds_full.rds")

#canada shapefile
canada <- geoboundaries("Canada")
canada <- st_transform(canada, "EPSG:4326")

#ecodistricts
ecodistricts <- sf::st_read("Canada/ecozones/ecodistrict_shp/Ecodistricts/ecodistricts.shp")
ecodistricts <- st_transform(ecodistricts, crs = "EPSG:4326")

#ecozones
ecozones <- sf::st_read("Canada/ecozones/ecozone_shp/Ecozones/ecozones.shp")
ecozones <- st_transform(ecozones, crs = "EPSG:4326")

#protected area data
parks <- readRDS("Canada/PAs/parks_fixed.rds")

#create color palettes
colors <- c('lightcyan', 'lightskyblue2', 'steelblue1', 'dodgerblue3', 'royalblue4', 'midnightblue')
colors.cv <- rev(c("#22052d","#3e2248","#5b3f64","#775c7f","#93799a","#b096b6","#ccb3d1"))

NDVI_cols <- rev(c("#0f2902", "#1d3900","#193401","#274009","#2e4511",
                   "#3d4f21", "#485921","#536321","#69761f","#868924",
                   "#8d8e37","#aaa263","#b5a975","#c2b58c","#c7b995",
                   "#cdbf9f","#e3d6c6","#e7dbce"))

#some data wrangling -------------------------------------------------------------
#calculate mean +  variance spatially 

VAR <- r %>%
  group_by(x, y, park) %>%
  summarize(mean = mean(preds),
            var = var(res),
            mean_res = mean(res),
            cv = cv(preds, aszero = T),
            cv.2 = sd(preds)/mean(preds),
            mean.2 = mean(NDVI))

VAR <- readRDS('Canada/VAR.rds')

VAR <- as.data.frame(VAR)

#create rasters for this object to plot easier 
rast <- rast(VAR, type = 'xyz')

mean.rast <- rast$mean
terra::writeCDF(mean.rast, 'Canada/mean_predNDVI_raster.nc', overwrite = TRUE)

var.rast <- rast$var
terra::writeCDF(var.rast, 'Canada/var_residuals_raster.nc', overwrite = TRUE)

meanres.rast <- rast$mean_res
terra::writeCDF(meanres.rast, 'Canada/mean_residuals_raster.nc', overwrite = TRUE)

cv.rast <- rast$cv
terra::writeCDF(cv.rast, 'Canada/coeff_variation_raster.nc', overwrite = TRUE)

meanNDVI.rast <- rast$mean.2
terra::writeCDF(meanNDVI.rast, 'Canada/mean_rawNDVI_raster.nc', overwrite = TRUE)

RES <- r %>%
  group_by(year, doy) %>%
  summarize(mean = mean(res))

saveRDS(RES, "Canada/results/RES.rds")

ELEV <- r %>%
  group_by(elevation) %>%
  summarize(mean = mean(res))

saveRDS(ELEV, "Canada/results/ELEV.rds")

MEAN.DOY <- r %>%
  group_by(doy, park) %>%
  summarize(mean = mean(preds))

saveRDS(MEAN.DOY, "Canada/results/MEAN.DOY.rds")

MEAN.YEAR <- r %>%
  group_by(year, park) %>%
  summarize(mean = mean(preds))

saveRDS(MEAN.YEAR, "Canada/results/MEAN.YEAR.rds")

VAR.DOY <- r %>%
  group_by(doy, park) %>%
  summarize(var = var(res))

saveRDS(VAR.DOY, "Canada/results/VAR.DOY.rds")

VAR.YEAR <- r %>%
  group_by(year, park) %>%
  summarize(var = var(res))

saveRDS(VAR.YEAR, "Canada/results/VAR.YEAR.rds")

CV.YEAR <- r %>%
  group_by(year, park) %>%
  summarize(cv = cv(preds, aszero = T))

saveRDS(CV.YEAR, "Canada/results/CV.YEAR.rds")

CV.DOY <- r %>%
  group_by(doy, park) %>%
  summarize(cv = cv(preds, aszero = T))

saveRDS(CV.DOY, "Canada/results/CV.DOY.rds")

ECOZ <- r %>%
  group_by(ecozone) %>%
  summarize(mean = mean(NDVI),
            var = var(preds),
            cv = cv(preds, aszero = T))

#plot data (figure 1) ------------------------------------------------------------

parks.plot <- ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = mean.2)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  geom_sf(data = parks, fill = NA, color = "black", linewidth = 0.2) +
  scale_fill_gradientn('Mean NDVI', colours = NDVI_cols) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, -0.5, -0.2, -0.5), "cm"), #top, right, bottom, left
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5, face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 10,
                               barheight = 0.2, title="Mean NDVI",  direction = "horizontal"))

ecoz.plot <- ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = mean.2)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  geom_sf(data = ecozones, fill = NA, color = "black") +
  scale_fill_gradientn('Mean NDVI', colours = NDVI_cols) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, -0.5, -0.2, -0.5), "cm"), #top, right, bottom, left
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) 

ecod.plot <- ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = mean.2)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  geom_sf(data = ecodistricts, fill = NA, color = "black") +
  scale_fill_gradientn('Mean NDVI', colours = NDVI_cols) +
  annotation_scale(location = "bl",
                   width_hint = 0.2,
                   pad_x = unit(0.6, "in"),
                   pad_y = unit(0.35, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         height = unit(0.9, "cm"),
                         width = unit(0.75, "cm"),
                         pad_y = unit(0.35, "in"),
                         pad_x = unit(0.2, "in"),
                         style = north_arrow_orienteering) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, -0.5, -0.2, -0.5), "cm"), #top, right, bottom, left
        legend.position = "none") 


ggarrange(parks.plot, ecoz.plot, ecod.plot, nrow = 3, common.legend = T, labels = "auto")

ggsave("figure1.png",
       units = "in",
       width = 3.23,
       height = 9.14,
       bg = "transparent",
       dpi = 600)

#plot mean trends (figure 2) -----------------------------------------------------

#spatial mean trends
mean <- 
  ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = mean)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn('Mean', colours = NDVI_cols) + 
  annotation_scale(location = "bl",
                   width_hint = 0.2,
                   pad_x = unit(1.2, "in"),
                   pad_y = unit(0.35, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         height = unit(0.9, "cm"),
                         width = unit(0.75, "cm"),
                         pad_y = unit(0.35, "in"),
                         pad_x = unit(0.4, "in"),
                         style = north_arrow_orienteering) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 20,
                               barheight = 0.5, title = "NDVI",  direction = "horizontal"))

#plot smooths

#separate into park and non-park data
doyparkyes <- MEAN.DOY[MEAN.DOY$park == 1,]
doyparkno <- MEAN.DOY[MEAN.DOY$park == 0,]

yearparkyes <- MEAN.YEAR[MEAN.YEAR$park == 1,]
yearparkno <- MEAN.YEAR[MEAN.YEAR$park == 0,]

#mean trends by year
parkymean <- 
  ggplot() +
    geom_point(yearparkyes, mapping = aes(year, mean), size = 0.3, color = "darkgreen") +
    geom_point(yearparkno, mapping = aes(year, mean), size = 0.3, color = "#A7C957") +
  geom_smooth(yearparkyes, mapping = aes(year, mean, colour = "Within Parks"), span = 0.25, se = FALSE) +
  geom_smooth(yearparkno, mapping = aes(year, mean, colour = "Outside Parks"),  span = 0.25, se = FALSE) +
  scale_colour_manual(name = "", values=c("#A7C957", "darkgreen")) +
  xlab("Year") +
  ylab("Mean NDVI") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.25, 0.9))

#mean trends by day of year
parkdoymean <- 
  ggplot() +
  geom_point(doyparkyes, mapping = aes(doy, mean), size = 0.3, color = "darkgreen") +
  geom_point(doyparkno, mapping = aes(doy, mean), size = 0.3, color = "#A7C957") +
  geom_smooth(doyparkyes, mapping = aes(doy, mean, colour = "Within Parks"), span = 0.1, se = FALSE) +
  geom_smooth(doyparkno, mapping = aes(doy, mean, colour = "Outside Parks"),  span = 0.1, se = FALSE) +
  scale_colour_manual(name = "", values=c("#A7C957", "darkgreen")) +
  xlab("Day of Year") +
  ylab("Mean NDVI") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

#plot mean trends by ecozone

#annotate ecozones into dataset
coords <- VAR[,c(1:2)]
ecoz <- unique(ecozones$ECOZONE)
ecozones$ECOZONE <- as.numeric(factor(ecozones$ECOZONE))
ecozones.rast <- rast(fasterize(ecozones, raster::raster(ecozones, resolution = 0.05), field = "ECOZONE"))
eco <- terra::extract(ecozones.rast, coords)
VAR <- cbind(VAR, eco)

VAR$layer <- as.character(VAR$layer)
names(VAR)[11] <- "ecozone"

saveRDS(VAR, "Canada/results/VAR.rds")

#boxplot to plot mean in different ecozones
boxmean <- 
  ggplot(edata, aes(x = layer, y = mean, fill = park)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.2) +
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values=c("#A7C957", "darkgreen")) +
  scale_x_discrete(labels=c("2" = "Northern Arctic", "1" = "Arctic Cordillera", "3" = "Southern Arctic", "9" = "Boreal Plain",
                            "6" = "Boreal Shield", "15" = "Hudson Plain", "14" = "Montane Cordillera", "8" = "Mixedwood Plain",
                            "13" = "Pacific Maritime", "10" = "Prairies", "11" = "Taiga Cordillera", "4" = "Taiga Plain",
                            "5" = "Taiga Shield", "7" = "Atlantic Maritime", "12" = "Boreal Cordillera")) +
  xlab("Ecozone") +
  ylab("Mean NDVI") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1, size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.7, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.1, 0.9))

bc <- ggarrange(parkymean, parkdoymean, nrow = 2, labels = c("b", "c"))
abc <- ggarrange(mean, bc, ncol = 2, labels = c("a", ""), widths = c(0.75, 0.5))
ggarrange(abc, boxmean, nrow = 2, labels = c("", "d"), heights = c(0.75, 0.5))

ggsave('figure2.png',
       units = "in",
       width = 6.86, 
       height = 8.52,
       bg = "transparent",
       dpi = 600)

#plot variance trends (figure 3) -------------------------------------------------

#spatial variance trends
var <- ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = var)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn('Variance in NDVI', colours = colors, limits = c(0, 0.05)) +
  annotation_scale(location = "bl",
                   width_hint = 0.2,
                   pad_x = unit(1.2, "in"),
                   pad_y = unit(0.35, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         height = unit(0.9, "cm"),
                         width = unit(0.75, "cm"),
                         pad_y = unit(0.35, "in"),
                         pad_x = unit(0.4, "in"),
                         style = north_arrow_orienteering) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 18,
                               barheight = 0.5, title = "Variance in NDVI",  direction = "horizontal"))

#plot smooths

#separate data by within/outside parks
doyparkyes <- VAR.DOY[VAR.DOY$park == 1,]
doyparkno <- VAR.DOY[VAR.DOY$park == 0,]

yearparkyes <- VAR.YEAR[VAR.YEAR$park == 1,]
yearparkno <- VAR.YEAR[VAR.YEAR$park == 0,]

parkyvar <- 
  ggplot() +
  geom_point(yearparkyes, mapping = aes(year, var), size = 0.3, color = "dodgerblue3") +
  geom_point(yearparkno, mapping = aes(year, var), size = 0.3, color = "lightskyblue2") +
  geom_smooth(yearparkyes, mapping = aes(year, var, colour = "Within Parks"), span = 0.25, se = FALSE) +
  geom_smooth(yearparkno, mapping = aes(year, var, colour = "Outside Parks"),  span = 0.25, se = FALSE) +
  scale_colour_manual(name = "", values=c('lightskyblue2', 'dodgerblue3')) +
  xlab("Year") +
  ylab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.25, 0.9))

#residual trends by day of year
parkdoyvar <- 
  ggplot() +
  geom_point(doyparkyes, mapping = aes(doy, var), size = 0.3, color = "dodgerblue3") +
  geom_point(doyparkno, mapping = aes(doy, var), size = 0.3, color = "lightskyblue2") +
  geom_smooth(doyparkyes, mapping = aes(doy, var, colour = "Within Parks"), span = 0.15, se = FALSE) +
  geom_smooth(doyparkno, mapping = aes(doy, var, colour = "Outside Parks"),  span = 0.15, se = FALSE) +
  scale_colour_manual(name = "", values=c('lightskyblue2', 'dodgerblue3')) +
  xlab("Day of Year") +
  ylab("Variance in NDVI") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

#plot mean trends by ecozone

#boxplot to plot mean in different ecozones
boxvar <- 
  ggplot(edata, aes(x = layer, y = var, fill = park)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.2) +
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values=c('lightskyblue2', 'dodgerblue3')) +
  scale_x_discrete(labels=c("2" = "Northern Arctic", "1" = "Arctic Cordillera", "3" = "Southern Arctic", "9" = "Boreal Plain",
                            "6" = "Boreal Shield", "15" = "Hudson Plain", "14" = "Montane Cordillera", "8" = "Mixedwood Plain",
                            "13" = "Pacific Maritime", "10" = "Prairies", "11" = "Taiga Cordillera", "4" = "Taiga Plain",
                            "5" = "Taiga Shield", "7" = "Atlantic Maritime", "12" = "Boreal Cordillera")) +
  xlab("Ecozone") +
  ylab("Variance in NDVI") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1, size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.7, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.1, 0.9))

bc <- ggarrange(parkyvar, parkdoyvar, nrow = 2, labels = c("b", "c"))
abc <- ggarrange(var, bc, ncol = 2, labels = c("a", ""), widths = c(0.75, 0.5))
ggarrange(abc, boxvar, nrow = 2, labels = c("", "d"), heights = c(0.75, 0.5))

ggsave('figure3.png',
       units = "in",
       width = 6.86, 
       height = 8.52,
       bg = "transparent",
       dpi = 600)

#plot quantiles (figure 4) ---------------------------------------------------------

#calculate quantiles
quantile(VAR$mean, probs = 0.7) #0.1127625 
quantile(VAR$var, probs = 0.3, na.rm = TRUE) #0.002437498 
quantile(VAR$cv, probs = 0.3) 

VAR$mean.quant <- with(VAR, ifelse(mean >= quantile(VAR$mean, probs = 0.7), 1, 0))
VAR$var.quant <- with(VAR, ifelse(var <= quantile(VAR$var, probs = 0.3, na.rm = TRUE), 1, 0))
VAR$cv.quant <- with(VAR, ifelse(cv <= quantile(VAR$cv, probs = 0.3), 1, 0))

VAR$mean.quant <- as.factor(VAR$mean.quant)
VAR$var.quant <- as.factor(VAR$var.quant)
VAR$cv.quant <- as.factor(VAR$cv.quant)

#save quantiles as rasters to speed up subsequent analyses
rast <- rast(VAR, type = 'xyz')

mean.quant.rast <- rast$mean.quant
terra::writeCDF(mean.quant.rast, 'Canada/mean_predNDVI_quantile_raster.nc', overwrite = TRUE)

var.quant.rast <- rast$var.quant
terra::writeCDF(var.quant.rast, 'Canada/mean_predNDVI_quantile_raster.nc', overwrite = TRUE)

cv.quant.rast <- rast$cv.quant
terra::writeCDF(cv.quant.rast, 'Canada/mean_predNDVI_quantile_raster.nc', overwrite = TRUE)

mean.quant <- 
  ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = mean.quant, alpha = mean.quant)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_manual(values = c('grey80',"darkgreen"), labels = c('bottom 70th quantile', 'top 30th quantile')) + 
    scale_alpha_discrete(guide = "none", range = c(0.6, 1), labels = c('bottom 70th quantile', 'top 30th quantile')) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.15, 0.95),
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

var.quant <- ggplot() +
  geom_raster(subset(VAR, !is.na(var.quant)), mapping = aes(x, y, fill = var.quant, alpha = var.quant)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_manual(values = c('grey80', 'dodgerblue3'), labels = c('top 70th quantile', 'bottom 30th quantile')) + 
  scale_alpha_discrete(guide = "none", range = c(0.6, 1), labels = c('bottom 70th quantile', 'top 30th quantile')) +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.15, 0.95),
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

cv.quant <- ggplot() +
  geom_raster(subset(VAR, !is.na(cv.quant)), mapping = aes(x, y, fill = cv.quant, alpha = cv.quant)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_manual(values = c("grey80", "#5b3f64"), labels = c('top 70th quantile', 'bottom 30th quantile')) + 
  scale_alpha_discrete(guide = "none", range = c(0.6, 1), labels = c('bottom 70th quantile', 'top 30th quantile')) +
  annotation_scale(location = "bl",
                   width_hint = 0.2,
                   pad_x = unit(0.7, "in"),
                   pad_y = unit(0.35, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         height = unit(0.9, "cm"),
                         width = unit(0.75, "cm"),
                         pad_y = unit(0.35, "in"),
                         pad_x = unit(0.1, "in"),
                         style = north_arrow_orienteering) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.15, 0.95),
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=7, family = "sans", face = "bold"),
        plot.margin = unit(c(0.2,-0.5,0.1,0.2), "cm")) #top, right, bottom, left 

#violin plots for within vs outside parks

#some data wrangling

VAR2 <- VAR
VAR2$park <- "2"

VAR.q <- rbind(VAR, VAR2)
VAR.q$park <- as.factor(VAR.q$park)

vio.mean <- ggplot(VAR.q, aes(park, mean, fill = park)) +
  geom_violin(adjust = 2) + 
  scale_fill_manual(values = c("#A7C957", "darkgreen", 'gray70')) + 
  ylab("Mean NDVI") +
    scale_x_discrete(labels=c("0" = "Outside Parks", "1" = "Within Parks", "2" = "Canada")) + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

vio.var <- ggplot(VAR.q, aes(park, var, fill = park)) +
  geom_violin(adjust = 2) + 
  scale_fill_manual(values = c('lightskyblue2', 'dodgerblue3', 'gray70')) +
  ylab("Variance in NDVI") +
  scale_x_discrete(labels=c("0" = "Outside Parks", "1" = "Within Parks", "2" = "Canada")) + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

vio.cv <- ggplot(VAR.q, aes(park, cv/100, fill = park)) +
  geom_violin(adjust = 2) + 
  scale_fill_manual(values = c("#ccb3d1", "#5b3f64", 'gray70')) +
  ylab("Coefficient of Variation") +
  scale_x_discrete(labels=c("0" = "Outside Parks", "1" = "Within Parks", "2" = "Canada")) + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

ggarrange(mean.quant, vio.mean, var.quant, vio.var, cv.quant, vio.cv,
          nrow = 3, ncol = 2, labels = "auto")

ggsave('figure4.png',
       units = "in",
       bg = "transparent",
       width = 6.86,
       height = 8.5,
       dpi = 600)

#plot model residuals (appendix) ----------------------------------------------

#spatial log mean of residuals
logmeanres <- ggplot() +
  geom_raster(subset(VAR, !is.na(mean_res)), mapping = aes(x, y, fill = log(mean_res + 2*abs(min(mean_res))))) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn("log Mean Residuals", colours = colors) + #limits = c(0, 0.15)) +
  annotation_scale(location = "bl",
                   width_hint = 0.2,
                   pad_x = unit(1.2, "in"),
                   pad_y = unit(0.35, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         height = unit(0.9, "cm"),
                         width = unit(0.75, "cm"),
                         pad_y = unit(0.35, "in"),
                         pad_x = unit(0.4, "in"),
                         style = north_arrow_orienteering) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 18,
                               barheight = 0.5, title = "log Mean Residuals",  direction = "horizontal"))

#mean of residuals for each area - plots faster than all residuals
mrdoy <- ggplot(RES) + 
  geom_point(aes(doy, mean, alpha = 0.1)) +
  xlab("Day of Year") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

mryear <- ggplot(RES) + 
  geom_point(aes(year, mean, alpha = 0.1)) +
  xlab("Year") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

mrlat <- ggplot(VAR) + 
  geom_point(aes(y, mean, alpha = 0.1)) +
  xlab("Latitude") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

mrlong <- ggplot(VAR) + 
  geom_point(aes(x, mean, alpha = 0.1)) +
  xlab("Longitude") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

mrelev <- ggplot(ELEV) + 
  geom_point(aes(elevation, mean, alpha = 0.1)) +
  xlab("Elevation (m)") +
  ylab("Mean of Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

bcf <- ggarrange(mrlat, mrlong, mrelev, nrow = 3, labels = c("b", "c", "f"))
de <- ggarrange(mrdoy, mryear, ncol = 2, labels = c("d", "e"))
ade <- ggarrange(logmeanres, de, nrow = 2, labels = c("a", ""), heights = c(1, 0.5))
ggarrange(ade, bcf, ncol = 2, widths = c(1, 0.5))

ggsave('residuals.png',
       units = "in",
       bg = "transparent",
       width = 6.86, 
       height = 6.86,
       dpi = 600)

#scatterplot of mean vs var
ggplot() +
  geom_point

#histogram of model residuals
hist(r$res, main = "", xlab = "Modelled Residuals")



#coefficient of variation (appendix) ---------------------------------------------------------

#spatial variance trends
cv <- 
  ggplot() +
  geom_raster(VAR, mapping = aes(x, y, fill = cv/100)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn('Coefficient of Variation', colours = colors.cv) + #limits = c(-0.54, 3.69)) +
  annotation_scale(location = "bl",
                   width_hint = 0.2,
                   pad_x = unit(1.2, "in"),
                   pad_y = unit(0.35, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         height = unit(0.9, "cm"),
                         width = unit(0.75, "cm"),
                         pad_y = unit(0.35, "in"),
                         pad_x = unit(0.4, "in"),
                         style = north_arrow_orienteering) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.4, 0.2), "cm"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) + 
  guides(fill = guide_colorbar(title.position = "top", ticks.colour = NA, barwidth = 18,
                               barheight = 0.5, title = "Coefficient of Variation",  direction = "horizontal"))

#plot smooths

#separate data by within/outside parks
doyparkyes <- CV.DOY[CV.DOY$park == 1,]
doyparkno <- CV.DOY[CV.DOY$park == 0,]

yearparkyes <- CV.YEAR[CV.YEAR$park == 1,]
yearparkno <- CV.YEAR[CV.YEAR$park == 0,]

parkycv <- 
  ggplot() +
  geom_point(yearparkyes, mapping = aes(year, cv/100), size = 0.3, color = "#5b3f64") +
  geom_point(yearparkno, mapping = aes(year, cv/100), size = 0.3, color = "#ccb3d1") +
  geom_smooth(yearparkyes, mapping = aes(year, cv/100, colour = "Within Parks"), span = 0.25, se = FALSE) +
  geom_smooth(yearparkno, mapping = aes(year, cv/100, colour = "Outside Parks"),  span = 0.25, se = FALSE) +
  scale_colour_manual(name = "", values = c("#ccb3d1", "#5b3f64")) +
  xlab("Year") +
  ylab("Coefficient of Variation") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size = 6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.25, 0.2))

#residual trends by day of year
parkdoycv <- 
  ggplot() +
  geom_point(doyparkyes, mapping = aes(doy, cv/100), size = 0.3, color = "#5b3f64") +
  geom_point(doyparkno, mapping = aes(doy, cv/100), size = 0.3, color = "#ccb3d1") +
  geom_smooth(doyparkyes, mapping = aes(doy, cv/100, colour = "Within Parks"), span = 0.15, se = FALSE) +
  geom_smooth(doyparkno, mapping = aes(doy, cv/100, colour = "Outside Parks"),  span = 0.15, se = FALSE) +
  scale_colour_manual(name = "", values = c("#ccb3d1", "#5b3f64")) +
  xlab("Day of Year") +
  ylab("Coefficient of Variation") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm"), #top, left, bottom, right
        legend.position = "none")

#plot mean trends by ecozone

#boxplot to plot mean in different ecozones
boxcv <- 
  ggplot(edata, aes(x = layer, y = cv/100, fill = park)) +
  geom_boxplot(outlier.size = 0.3, lwd = 0.2) +
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values = c("#ccb3d1", "#5b3f64")) +
  scale_x_discrete(labels=c("2" = "Northern Arctic", "1" = "Arctic Cordillera", "3" = "Southern Arctic", "9" = "Boreal Plain",
                            "6" = "Boreal Shield", "15" = "Hudson Plain", "14" = "Montane Cordillera", "8" = "Mixedwood Plain",
                            "13" = "Pacific Maritime", "10" = "Prairies", "11" = "Taiga Cordillera", "4" = "Taiga Plain",
                            "5" = "Taiga Shield", "7" = "Atlantic Maritime", "12" = "Boreal Cordillera")) +
  xlab("Ecozone") +
  ylab("Coefficient of Variation") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 9, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 9, family = "sans", face = "bold"),
        axis.text.y = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1, size = 8, family = "sans"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.7, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.1,0.2), "cm"), #top, left, bottom, right
        legend.position = c(0.1, 0.9))

bc <- ggarrange(parkycv, parkdoycv, nrow = 2, labels = c("b", "c"))
abc <- ggarrange(cv, bc, ncol = 2, labels = c("a", ""), widths = c(0.75, 0.5))
ggarrange(abc, boxcv, nrow = 2, labels = c("", "d"), heights = c(0.75, 0.5))

ggsave('cv.png',
       units = "in",
       width = 6.86, 
       height = 8.52,
       bg = "transparent",
       dpi = 600)


