library('dplyr')
library('ggplot2') #for figures
library('ggpubr') #for arranging plots
library('sf') #for importing ecozone data
library('terra') #for basemaps
library('rgeoboundaries') #for canada shapefile

#import data - models + final dataset-----------------------------------------------

#complete dataset with predicted values for each observation
p <- readRDS("Canada/predictions_march6.rds")

#variance of each individual coordinate
results <- readRDS("Canada/results_march20.rds")

#canadian ecozones, merge arctic ecozones into one
ecozones <- sf::st_read("Canada/ecozones/ecozone_shp/Ecozones/ecozones.shp")
ecozones <- st_transform(ecozones, crs = "EPSG:4526")
arctics <- c("Arctic Cordillera", "Northern Arctic", "Southern Arctic")
ecozones$ZONE_NAME <- replace(ecozones$ZONE_NAME, ecozones$ZONE_NAME %in% arctics, "Arctic")

#canadian protected areas, limit to only terrestrial larger than 50 sq km
parks <- sf::st_read("Canada/PAs/ProtectedConservedArea.gdb")
parks <- st_transform(parks, crs = "EPSG:4526")
marine <- parks$TYPE_E[grepl("Marine", parks$TYPE_E)] #extract marine protected areas
proposed <- parks$TYPE_E[grepl("Proposed", parks$TYPE_E)] #extract proposed protected areas
st_geometry(parks) <- "GEOMETRY"
terr.parks <- parks[(parks$TYPE_E %in% marine) == F,]
terr.parks <- terr.parks[(terr.parks$TYPE_E %in% proposed) == F,]
lg.parks <- terr.parks[(terr.parks$O_AREA_HA >= 5000),] #limit to only parks larger than 50km

#canada shapefile
canada <- geoboundaries("Canada")
canada <- st_transform(canada, "EPSG:4326")

colors <- c('lightcyan', 'lightskyblue2', 'steelblue1', 'dodgerblue3', 'royalblue4', 'midnightblue')
NDVI_cols <- colorRampPalette(rev(c("#0f2902", "#1d3900","#193401","#274009","#2e4511",
                                    "#3d4f21", "#485921","#536321","#69761f","#868924",
                                    "#8d8e37","#aaa263","#b5a975","#c2b58c","#c7b995",
                                    "#cdbf9f","#e3d6c6","#e7dbce")))
eco_cols <- c("#d3d3d3","#005f73","#0a9396","#94d2bd", "#def7ec","#e9d8a6","#b58463", "#f9c74f","#d4a373",
              "#ca6702","#bb3e03","#df8080","#9b2226")

#check model residuals----------------------------------------------------------

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

#plot general trends + data---------------------------------------------------------------

#mean productivity across space
ggplot() +
  geom_raster(subset(results, !is.na(mu)), mapping = aes(long, lat, fill = mu), interpolate = TRUE) +
  geom_sf(data = lg.parks, fill = NA, color = "black") + #add parks
  geom_sf(data = canada, fill = NA, color = "black") + #add canada borders
  scale_fill_gradientn('Variance', colours = colors) +
  annotation_scale(location = "bl", 
                   width_hint = 0.5,
                   pad_x = unit(0.18, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         pad_x = unit(0.2, "in"), 
                         pad_y = unit(4, "in"),
                         style = north_arrow_fancy_orienteering) +
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


#variance across space
ggplot() +
  geom_raster(subset(results, !is.na(var)), mapping = aes(long, lat, fill = var), interpolate = TRUE) +
  geom_sf(data = lg.parks, fill = NA, color = "black") + #add parks
  geom_sf(data = canada, fill = NA, color = "black") + #add canada borders
  scale_fill_gradientn('Variance', colours = colors) +
  annotation_scale(location = "bl", 
                   width_hint = 0.5,
                   pad_x = unit(0.18, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         pad_x = unit(0.2, "in"), 
                         pad_y = unit(4, "in"),
                         style = north_arrow_fancy_orienteering) +
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

#plot all ecozones
ggplot() +
  geom_sf(data = ecozones, aes(fill = ZONE_NAME), color = "white") +
  scale_fill_discrete(name = "Ecozones", type = eco_cols) + 
  geom_sf(data = canada, fill = NA, color = "black") + #add canada outline
  geom_sf(data = lg.parks, fill = NA, color = "black") + #add protected areas outlines
  coord_sf(crs = "EPSG:4326") +
  annotation_scale(location = "bl", 
                   width_hint = 0.2,
                   pad_x = unit(0.18, "in")) +
  annotation_north_arrow(location = "bl",
                         which_north = "true", 
                         pad_x = unit(0.5, "in"),
                         pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_text("Ecozone"),
        legend.key.size = unit(0.5, "cm"))

ggsave("ecoparks.png",
       width = 7.70,
       height = 5.76,
       dpi = 300,
       bg = "transparent")


#general trends within/outside parks---------------------------------------------------

#subset predicted data within vs outside parks
parkyes <- p[p$park == "TRUE",]
parkno <- p[p$park == "FALSE",]

#mean trends over time within/outside parks
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

#mean daily trends over a year within/outside parks
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
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values=c("lightskyblue1", "dodgerblue3")) +
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

plot2 <- ggarrange(plot1, boxvar,
          ncol = 2,
          labels = c("", "c"))

ggarrange(plot2, ecoboxvar, nrow = 2, labels = c("", "d"))

ggsave("parktrends.png",
       width = 6.8,
       height = 8.75,
       dpi = 450,
       bg = "transparent")

#boxplot comparing variance within/outside parks by ecozone
ecoboxvar <- ggplot(results, aes(x = ecozone, y = var, fill = park)) +
  geom_boxplot() +
  scale_fill_manual(name = "", labels = c("Outside parks", "Within parks"), values=c("lightskyblue1", "dodgerblue3")) +
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

#function to create figure of specified ecozone
zone.plot <- function(zone = 'ARC', name = "Arctic", height = c(2, 1)){
  
  #subset data to single ecozone
  zn <- results[results$ecozone == zone,]
  
  #plot variance across ecozone
  varplot <- ggplot() +
    geom_raster(data = subset(zn, !is.na(var)), aes(long, lat, fill = var), interpolate = TRUE) +
    geom_sf(data = parks, fill = NA, color = "black") + 
    geom_sf(data = canada) + 
    coord_sf(xlim = c(min(zn$long), max(zn$long)), ylim = c(max(zn$lat), min(zn$lat)), expand = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn('Variance', colours = colors, limits = c(0, 0.2)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.key.size = unit(1.5, "cm"))
  
  mean <- filter(p, ecozone == zone)
  
  #subset predicted mean data for within/outside parks
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
    scale_fill_manual(labels = c("Outside parks", "Within parks"), values=c("lightskyblue1", "dodgerblue3")) +
    xlab("Ecozone") +
    ylab("Variance") +
    scale_x_discrete(labels = name) +
    #scale_y_discrete(expand = c(0, 0)) +
    #scale_fill_brewer(palette = "Dark2") +
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

zone.plot(zone = "TGC", name = "Taiga Cordillera", height = c(2, 1))

ggsave("taigacordillera.png",
       width = 8,
       height = 7.50,
       dpi = 300)

