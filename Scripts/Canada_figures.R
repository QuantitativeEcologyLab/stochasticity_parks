library('dplyr')
library('ggplot2')
library('ggpubr')
library('lubridate') #for decimal_date function
library('mgcv') #for predict.gam
library('sf') #for importing ecozone data
library('terra') #for basemaps
library('rgeoboundaries') #for canada shapefile
library('basemaps') #for figure basemaps
library('tidyterra')
library('ggspatial')

#colors <- c('lightcyan', 'lightskyblue2', 'steelblue1', 'dodgerblue3', 'royalblue4', 'midnightblue')
#colors <- rev(c("#22052d","#3e2248","#5b3f64","#775c7f","#93799a","#b096b6","#ccb3d1"))
colors <- c("#f7b267","#f79d65","#f4845f","#f27059","#f25c54", "#d7263d")

NDVI_cols <- rev(c("#0f2902", "#1d3900","#193401","#274009","#2e4511",
                   "#3d4f21", "#485921","#536321","#69761f","#868924",
                   "#8d8e37","#aaa263","#b5a975","#c2b58c","#c7b995",
                   "#cdbf9f","#e3d6c6","#e7dbce"))

#import data -----------------------------------------------

#complete dataset with predicted values for each observation
p <- readRDS("Canada/predictions_march6.rds")

#canada shapefile
canada <- geoboundaries("Canada")
canada <- st_transform(canada, "EPSG:4326")

#final dataset with predicted values
p <- readRDS('Canada/Data_annotated/50km/predictions_april10.rds')

#variance of each individual coordinate
results <- readRDS("Canada/Data_annotated/50km/results_march20.rds")

#canadian ecozones, merge arctic ecozones into one
ecozones <- sf::st_read("Canada/ecozones/ecozone_shp/Ecozones/ecozones.shp")
ecozones <- st_transform(ecozones, crs = "EPSG:4526")
arctics <- c("Arctic Cordillera", "Northern Arctic", "Southern Arctic")
ecozones$ZONE_NAME <- replace(ecozones$ZONE_NAME, ecozones$ZONE_NAME %in% arctics, "Arctic")

#wrangle protected area data ----------------------------------
parks <- sf::st_read("Canada/PAs/ProtectedConservedArea.gdb")
parks <- st_transform(parks, crs = "EPSG:4326")

parks <- parks[st_is(parks, c("POLYGON", "MULTIPOLYGON")),] #remove parks with complex geometries (cannot be converted to polygon)

parks <- st_make_valid(parks)
parks <- st_intersection(parks, canada) #crop parks to only ones in terrestrial canada

#remove additional marine and proposed protected areas
marine <- parks$TYPE_E[grepl("Marine", parks$TYPE_E)]
proposed <- parks$TYPE_E[grepl("Proposed", parks$TYPE_E)]
#st_geometry(parks) <- "GEOMETRY"
parks <- parks[(parks$TYPE_E %in% marine) == F,]
parks <- parks[(parks$TYPE_E %in% proposed) == F,]

parks <- parks[(parks$O_AREA_HA >= 5000),] #limit to only parks larger than 50km

#plot model residuals----------------------------------------------------------

#residuals vs day of year
rdoy <- ggplot(p) + 
  geom_point(aes(doy, resid)) +
  xlab("Day of Year") +
  ylab("Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

#residuals vs year
ryear <- ggplot(p) + 
  geom_point(aes(year, resid)) +
  xlab("Year") +
  ylab("Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

#residuals vs latitude
rlat <- ggplot(p) +
  geom_point(aes(lat, resid)) +
  xlab("Latitude") +
  ylab("Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

#residuals vs longitude
rlong <- ggplot(p) + 
  geom_point(aes(long, resid)) +
  xlab("Longitude") +
  ylab("Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

#residuals vs quality assurance
rqa <- ggplot(p) + 
  geom_point(aes(QA, resid)) +
  xlab("Quality Assurance") +
  ylab("Residuals") + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

#residuals vs elevation
relev <- ggplot(p) +
  geom_point(aes(elevation, resid)) +
  xlab("Elevation") +
  ylab("Residuals") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

ggarrange(rdoy, ryear, rlat, rlong, rqa, relev,
          nrow = 2, ncol = 3) 

ggsave("residuals.png",
       width = 10,
       height = 7,
       dpi = 300)

#check residuals
appraise(m) &
  theme_classic() &
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

ggsave("model.appraise.png",
       width = 10,
       height = 7,
       dpi = 300)   

#plot general trends---------------------------------------------------------------

#variance across space
ggplot() +
  geom_raster(subset(results, !is.na(var)), mapping = aes(long, lat, fill = var)) +
  geom_sf(data = parks, fill = NA, color = "black") + 
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn('Variance', colours = colors, limits = c(0, 0.15)) +
  annotation_scale(location = "bl", 
                   width_hint = 0.5,
                   pad_x = unit(0.18, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(0.9, 0.75))

ggsave("variance.png",
       width = 9.02,
       height = 7.38,
       dpi = 300)

#mean trends across space

ggplot() +
  geom_raster(subset(results, !is.na(var)), mapping = aes(long, lat, fill = mean)) +
  geom_sf(data = parks, fill = NA, color = "black") + 
  geom_sf(data = canada, fill = NA, color = "black") + 
  scale_fill_gradientn('Mean', colours = NDVI_cols) + #limits = c(0, 0.15)) +
  annotation_scale(location = "bl", 
                   width_hint = 0.5,
                   pad_x = unit(0.18, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(datum = "ESRI:102001") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(0.9, 0.75))

ggsave("mean.png",
       width = 9.02,
       height = 7.38,
       dpi = 300)


#plot all ecozones

ggplot() +
  geom_raster(data = p, aes(x = long, y = lat, fill = ecozone)) +
  scale_fill_discrete(name = "Ecozones") + #limits = c(0, 1)) +
  geom_sf(data = canada, fill = NA, color = "black") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_text("Ecozone"),
        legend.key.size = unit(0.5, "cm"))

ggsave("ecozones.png",
       width = 7.70,
       height = 5.76,
       dpi = 300,
       bg = "transparent")


#general trends within/outside parks---------------------------------------------------

#generate predicted data for mean trends across time
parkyes <- p[p$park == "TRUE",]
parkno <- p[p$park == "FALSE",]

parkymean <- ggplot() +
  geom_smooth(parkyes, mapping = aes(year, mu, colour = "Within Parks"), fill = "#A7C957") +
  geom_smooth(parkno, mapping = aes(year, mu, colour = "Outside Parks"), fill = "#CCD5AE") +
  scale_colour_manual(name = "", values=c("#A7C957", "darkgreen")) +
  xlab("Year") +
  ylab("Mean NDVI") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text("Year"),
        axis.title.y = element_text("Mean"),
        legend.position = c(0.25, 0.8),
        legend.background = element_blank())

parkdoymean <- ggplot() +
  geom_smooth(parkyes, mapping = aes(doy, mu, colour = "Within Parks"), fill = "#A7C957") +
  geom_smooth(parkno, mapping = aes(doy, mu, colour = "Outside Parks"), fill = "#CCD5AE") +
  scale_colour_manual(name = "", values=c("#A7C957", "darkgreen")) +
  xlab("Day of Year") +
  ylab("Mean NDVI") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text("Year"),
        axis.title.y = element_text("Mean"),
        legend.position = c(0.25, 0.8),
        legend.background = element_blank())

#boxplot comparing variance within/outside parks

boxvar <- ggplot(results, aes(x = park, y = var, fill = park)) +
  geom_boxplot() +
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values=c("#f7b267", "#d7263d")) +
  xlab("Ecozone") +
  ylab("Variability") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.5, 0.8))


plot1 <- ggarrange(parkymean, parkdoymean,
                   nrow = 2,
                   labels = "auto")

ggarrange(plot1, boxvar,
          # widths = c(1, 2),
          ncol = 2,
          labels = c("", "c"))

ggsave("parktrends.png",
       width = 8,
       height = 5.5,
       dpi = 450,
       bg = "transparent")

#boxplot comparing variance within/outside parks by ecozones

ggplot(results, aes(x = ecozone, y = var, fill = park)) +
  geom_boxplot() +
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values=c("#f7b267", "#d7263d")) +
  scale_x_discrete(labels=c("ARC" = "Arctic", "ATM" = "Atlantic Maritime", "BRC" = "Boreal Cordillera", "BRP" = "Boreal Plain",
                            "BRS" = "Boreal Shield", "HUD" = "Hudson Plain", "MTC" = "Montane Cordillera", "MXP" = "Mixedwood Plain",
                            "PCM" = "Pacific Maritime", "PRA" = "Prairies", "TGC" = "Taiga Cordillera", "TGP" = "Taiga Plain",
                            "TGS" = "Taiga Shield")) +
  xlab("Ecozone") +
  ylab("Variability") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=45, vjust = 0.95, hjust = 1),
        legend.position = c(0.25, 0.8))

ggsave("ecoparktrends.png",
       width = 7.5,
       height = 4.5,
       dpi = 450,
       bg = "transparent")

#plot trends by ecozone --------------------------------------------------------

#change names in ecozone data to match working dataset
arctics <- c("Arctic Cordillera", "Northern Arctic", "Southern Arctic")
ecozones$ZONE_NAME <- replace(ecozones$ZONE_NAME, ecozones$ZONE_NAME %in% arctics, "Arctic")
ecozones$zone[ecozones$ZONE_NAME == "Arctic"] <- "ARC"
ecozones$zone[ecozones$ZONE_NAME == "Taiga Cordillera"] <- "TGC"
ecozones$zone[ecozones$ZONE_NAME == "Taiga Plain"] <- "TGP"
ecozones$zone[ecozones$ZONE_NAME == "Taiga Shield"] <- "TGS"
ecozones$zone[ecozones$ZONE_NAME == "Hudson Plain"] <- "HUD"
ecozones$zone[ecozones$ZONE_NAME == "Boreal Cordillera"] <- "BRC"
ecozones$zone[ecozones$ZONE_NAME == "Boreal PLain"] <- "BRP"
ecozones$zone[ecozones$ZONE_NAME == "Boreal Shield"] <- "BRS"
ecozones$zone[ecozones$ZONE_NAME == "Prairie"] <- "PRA"
ecozones$zone[ecozones$ZONE_NAME == "Montane Cordillera"] <- "MTC"
ecozones$zone[ecozones$ZONE_NAME == "Pacific Maritime"] <- "PCM"
ecozones$zone[ecozones$ZONE_NAME == "Atlantic Maritime"] <- "ATM"
ecozones$zone[ecozones$ZONE_NAME == "MixedWood Plain"] <- "MXP"

zone.plot <- function(zone = 'ARC', name = "Arctic", height = c(2, 1)){
  
  #subset data to single ecozone
  zn <- results[results$ecozone == zone,]
  
  extent <- st_bbox(c(xmin = min(zn$long), xmax = max(zn$long), ymax = max(zn$lat), ymin = min(zn$lat)), crs = st_crs(4326))
  
  bm <- basemap_terra(ext = extent, map_service = "esri", map_type = "world_physical_map") %>%
    as("SpatRaster")
  
  bm <- project(bm, crs("EPSG:4326"))
  
  #plot variance across ecozone
  varplot <- ggplot() +
    geom_spatraster_rgb(data = bm) + 
    geom_raster(data = subset(zn, !is.na(var)), aes(long, lat, fill = var), interpolate = TRUE) +
    geom_sf(data = parks, fill = NA, color = "black") + 
    geom_sf(data = canada, fill = NA, color = "black") + 
    coord_sf(xlim = c(min(zn$long), max(zn$long)), ylim = c(max(zn$lat), min(zn$lat)), expand = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn('Variance', colours = colors, limits = c(0, 0.2)) +
    annotation_scale(location = "tl", 
                     width_hint = 0.5,
                     pad_x = unit(0.18, "in")) +
    annotation_north_arrow(location = "tl",
                           which_north = "true", 
                           pad_y = unit(0.5, "in"), 
                           style = north_arrow_fancy_orienteering) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.key.size = unit(1.5, "cm"))
  
  mean <- filter(p, ecozone == zone)
  
  #subset data for within/outside parks
  znparkyes <- mean[mean$park == "TRUE",]
  znparkno <- mean[mean$park == "FALSE",]
  
  #plot yearly trends
  parkymean <- ggplot() +
    geom_smooth(znparkyes, mapping = aes(year, mu, colour = "Within Parks"), fill = "#A7C957") +
    geom_smooth(znparkno, mapping = aes(year, mu, colour = "Outside Parks"), fill = "#CCD5AE") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Year") +
    ylab("Mean NDVI") +
    scale_colour_manual(name = "Legend", values=c("#A7C957", "darkgreen")) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text("Year"),
          axis.title.y = element_text("Mean"),
          legend.position = c(0.25, 0.9),
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(size = 8))
  
  #plot daily trends over a year
  parkdoymean <- ggplot() +
    geom_smooth(znparkyes, mapping = aes(doy, mu, colour = "Within Parks"), fill = "#A7C957") +
    geom_smooth(znparkno, mapping = aes(doy, mu, colour = "Outside Parks"), fill = "#CCD5AE") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Day of Year") +
    ylab("Mean NDVI") +
    scale_colour_manual(name = "Legend", values=c("#A7C957", "darkgreen")) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text("Year"),
          axis.title.y = element_text("Mean"),
          legend.position = c(0.25, 0.9),
          legend.key = element_blank(),
          legend.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))
  
  #plot average variance for ecozone
  ecoboxvar <- ggplot(zn, aes(x = ecozone, y = var, fill = park)) +
    geom_boxplot() +
    scale_fill_manual(labels = c("Outside parks", "Within parks"), values=c("#f7b267", "#d7263d")) +
    xlab("Ecozone") +
    ylab("Variance") +
    scale_x_discrete(labels = name) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.25, 0.9),
          legend.text = element_text(size = 8))
  
  #join figures together
  plot <- ggarrange(parkymean, parkdoymean, ecoboxvar, ncol = 3,
                    labels = c("b", "c", "d"))
  
  final <- ggarrange(varplot, plot,
                     nrow = 2,
                     labels = c("a", ""),
                     heights = height)
  
  annotate_figure(final, top = text_grob(name,
                                         color = "black", face = "bold", size = 20))
  
}

zone.plot(zone = "ATM", name = "Atlantic Maritime", height = c(2, 1))

ggsave("atlanticmaritime.png",
       width = 8,
       height = 7.50,
       dpi = 300)

