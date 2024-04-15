library('dplyr')
library('ggplot2')
library('ggpubr')
library('lubridate') #for decimal_date function
library('mgcv') #for predict.gam
library('sf') #for importing ecozone data - NOT WORKING
library('rgeoboundaries') #for canada shapefile

#import data - models + final dataset-----------------------------------------------

#complete dataset with predicted values for each observation
p <- readRDS("Canada/predictions_march6.rds")

#variance of each individual coordinate
results <- readRDS("Canada/results_march20.rds")

ecozones <- sf::st_read("Canada/ecozones/ecozone_shp/Ecozones/ecozones.shp")
ecozones <- st_transform(ecozones, crs = "EPSG:4326")
arctics <- c("Arctic Cordillera", "Northern Arctic", "Southern Arctic")
ecozones$ZONE_NAME <- replace(ecozones$ZONE_NAME, ecozones$ZONE_NAME %in% arctics, "Arctic")

parks <- sf::st_read("Canada/PAs/ProtectedConservedArea.gdb")
parks <- st_transform(parks, crs = "EPSG:4326")
marine <- parks$TYPE_E[grepl("Marine", parks$TYPE_E)] #extract marine protected areas
proposed <- parks$TYPE_E[grepl("Proposed", parks$TYPE_E)] #extract proposed protected areas
st_geometry(parks) <- "GEOMETRY"
terr.parks <- parks[(parks$TYPE_E %in% marine) == F,]
terr.parks <- terr.parks[(terr.parks$TYPE_E %in% proposed) == F,]
lg.parks <- terr.parks[(terr.parks$O_AREA_HA >= 5000),] #limit to only parks larger than 50km

canada <- geoboundaries("Canada")

colors <- c('lightcyan', 'lightskyblue2', 'steelblue1', 'dodgerblue3', 'royalblue4', 'midnightblue')
NDVI_cols <- colorRampPalette(rev(c("#0f2902", "#1d3900","#193401","#274009","#2e4511",
                                    "#3d4f21", "#485921","#536321","#69761f","#868924",
                                    "#8d8e37","#aaa263","#b5a975","#c2b58c","#c7b995",
                                    "#cdbf9f","#e3d6c6","#e7dbce")))

#check model residuals----------------------------------------------------------

rdoy <- ggplot(p) + #overpredicting in summer - more knots?
  geom_point(aes(doy, resid)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

ryear <- ggplot(p) + 
  geom_point(aes(year, resid)) 

rlat <- ggplot(p) +
  geom_point(aes(lat, resid)) 

rlong <- ggplot(p) + 
  geom_point(aes(long, resid))

rqa <- ggplot(p) + 
  geom_point(aes(QA, resid))

relev <- ggplot(p) +
  geom_point(aes(elevation, resid))

#doesnt work for some reason?
ggarrange(rdoy, ryear, rlat, rlong, rqa,
          nrow = 2, ncol = 3) #+ 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

#plot general trends---------------------------------------------------------------

#variance across space
ggplot() +
  geom_raster(data = subset(results, !is.na(var)), aes(long, lat, fill = var)) +
  scale_fill_gradientn('Variance', colours = colors) + #limits = c(0, 1)) +
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


#trends in mean

ymean <- ggplot(p, aes(year, mu)) +
  geom_smooth(color = 'darkgreen', fill = "#A7C957") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        axis.text.x = element_text("Year"),
        axis.text.y = element_text("Mean"))

doymean <- ggplot(p, aes(doy, mu)) +
  geom_smooth(color = 'darkgreen', fill = "#A7C957") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

canadamean <- ggplot() +
  geom_raster(data = results, aes(long, lat, fill = mean)) +
  scale_fill_gradientn('Mean', colours = ndvi_pal) + #limits = c(0, 1)) +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(0.9, 0.75))

trends <- ggarrange(ymean, doymean,
                    ncol = 2,
                    labels = c("b", "c"))

ggarrange(canadamean, trends,
          nrow = 2,
          labels = c("a", ""))

ggsave("trends.png",
       width = 6.47,
       height = 9.16,
       dpi = 300,
       bg = "transparent")

#plot all ecozones

eco <- ecozones[st_intersection(ecozones, canada) %>% lengths > 0,]

eco <- ecozones[st_difference(st_union(ecozones),st_union(canada)) %>% lengths > 0,]

#better but still not great, colors are ugly
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

#plot all parks

#change projection, lay on top of mean or var plot?
png(filename = "parks.png",
    width = 2000,
    height = 1000,
    res = 300)

plot(lg.parks$GEOMETRY)

dev.off()

png(filename = "candaparks.png",
    width = 2000,
    height = 1000,
    res = 300)
par(bg = NA)

plot(canada$geometry)

dev.copy(png, "canadaparks.png")
dev.off()

#plot CV across canada

latlong <- unique(p[,c(3, 8:10)])
res2 <- list()

for(i in 1:nrow(latlong)){
  
  point <- latlong[i,] #extract single point
  
  matches <- p[which(p$lat == point$lat & p$long == point$long),] #find all data related to single point in original dataset
  
  res <- data.frame(cv2 = var(matches$mu)/mean(matches$mu))
  
  res2[[i]] <- res
  
}

res2 <- do.call(rbind, res2)
res2 <- cbind(p, res2)

ggplot() +
  geom_raster(data = results, aes(long, lat, fill = cv)) +
  scale_fill_gradientn('Coefficient of Variation', colours = colors) +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(0.9, 0.75))

#general trends within/outside parks---------------------------------------------------

#generate predicted data for mean trends across time
parkyes <- p[p$park == "TRUE",]
parkno <- p[p$park == "FALSE",]

parkymean <- ggplot() +
  geom_smooth(parkyes, mapping = aes(year, mu, colour = "Within Parks"), fill = "#A7C957") +
  geom_smooth(parkno, mapping = aes(year, mu, colour = "Outside Parks"), fill = "#CCD5AE") +
  scale_colour_manual(name = "Legend", values=c("#A7C957", "darkgreen")) +
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
  scale_colour_manual(name = "Legend", values=c("#A7C957", "darkgreen")) +
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

boxvar <- ggplot(results, aes(x = ecozone, y = var, fill = park)) +
  geom_boxplot() +
  scale_fill_manual(name = "Legend", labels = c("Outside parks", "Within parks"), values=c("lightskyblue1", "dodgerblue3")) +
  xlab("Ecozone") +
  ylab("Variability") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.25, 0.8))

plot1 <- ggarrange(parkymean, parkdoymean,
                   nrow = 2,
                   labels = "auto")

ggarrange(plot1, boxvar,
          widths = c(1, 2),
          ncol = 2,
          labels = c("", "c"))


ggsave("parktrends.png",
       width = 10.63,
       height = 5.5,
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

#map needs some cleaning up but works! add basemap + park polygons
zone.plot <- function(zone = 'ARC', name = "Arctic", height = c(2, 1)){
  
  zn <- results[results$ecozone == zone,]
  
  bbox <- st_bbox(c(xmin = min(zn$long), xmax = max(zn$long), ymax = max(zn$lat), ymin = min(zn$lat)), crs = st_crs(4326))
  
  varplot <- ggplot() +
    geom_raster(data = subset(zn, !is.na(var)), aes(long, lat, fill = var)) +
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
  
  znparkyes <- mean[mean$park == "TRUE",]
  znparkno <- mean[mean$park == "FALSE",]
  
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
  
  plot <- ggarrange(parkymean, parkdoymean, ecoboxvar, ncol = 3,
                    labels = c("b", "c", "d"))
  
  final <- ggarrange(varplot, plot,
            nrow = 2,
            labels = c("a", ""),
            heights = height)
  
  annotate_figure(final, top = text_grob(name,
                                        color = "black", face = "bold", size = 20))
  
}

zone.plot(zone = "ARC", name = "Arctic", height = c(2, 1))

ggsave("pacificmaritime.png",
       width = 8,
       height = 7.50,
       dpi = 300)

#plots for presentation-----------------------------------------------------------
library('terra')
library('gridExtra')

single.rast <- rast("Canada/NDVI/NOAA_Files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810703_c20170609173647.nc")

png(filename = "NDVI_raster_single2.png",
    width = 2000,
    height = 1000,
    res = 300)

plot(single.rast$NDVI)

dev.off()

df <- p[c(1:5),c(1:10)]

png("dataframe.png", height= 800, width= 2700, res = 300)
g <-tableGrob(df)
grid.arrange(g)
dev.off()


#outdated --------------------------------------------------------------------

mm <- readRDS("Canada/models/bam_mean_1GB_mar3.rds")

mv <- readRDS("Canada/models/bam_var_1GB_mar4.rds")

parks <- sf::st_read("Canada/PAs/ProtectedConservedArea.gdb")

data <- readRDS("Canada/data_eco.rds")
daysthin <- seq(from = 1, to = 366, by = 5)
data.thin <- data[data$doy %in% daysthin,] #thin data to resolution of every 5 days
d <- mutate(data.thin, NDVI_scaled = (NDVI + 1) / 2) #scale NDVI to fit beta distribution (0 to 1)

#plot general trends---------------------------------------------------------------

DAYS <- unique(d$date)
RES <- list()

#create a loop to predict the new mean data from the model
for(i in 1:length(DAYS)){
  newd <- filter(d, date == DAYS[i])
  
  #add the predicted values back to the data frame with all the data
  preds <-
    bind_cols(newd,
              as.data.frame(predict.gam(mean_ndvi, newdata = newd, type = 'response')) %>%
                rename(mu = 'predict.gam(mean_ndvi, newdata = newd, type = "response")'))
  
  if(!is.na(preds$NDVI[1])){
    RES[[i]] <- data.frame(day = unique(preds$date),
                           mu = mean(preds$mu))
  }
  
}

RESULTSM <- do.call(rbind, RES)

DAYS <- unique(p$date)
RES <- list()

#create a loop to predict the new variance data from the model
for(i in 1:length(DAYS)){
  newd <- filter(p, date == DAYS[i])
  
  #add the predicted values back to the data frame with all the data
  preds <-
    bind_cols(newd,
              as.data.frame(predict.gam(var_ndvi, newdata = newd, type = 'response')) %>%
                rename(phi = 'predict.gam(var_ndvi, newdata = newd, type = "response")'))
  
  if(!is.na(preds$phi[1])){
    RES[[i]] <- data.frame(day = unique(preds$date),
                           phi = mean(preds$phi))
  }
  
}

RESULTSV <- do.call(rbind, RES)

res <- bind_cols(RESULTSM, RESULTSV) %>%
  mutate(sigma2 = phi * (1 - mu) * mu, #variance
         mu = mu * 2 - 1, #rescale ndvi
         var = sigma2 * 4) #scale variance appropriately

res <- res[,-3]
names(res)[1] <- "day"

res <- mutate(res,
              dec_date = decimal_date(day),
              year = year(day),
              doy = yday(day))

ymean <- ggplot(res, aes(year, mu)) +
  geom_smooth(color = 'darkgreen') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        axis.text.x = element_text("Year"),
        axis.text.y = element_text("Mean"))

yvar <- ggplot(res, aes(year, var)) +
  geom_smooth(color = 'dodgerblue3') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

doymean <- ggplot(res, aes(doy, mu)) +
  geom_smooth(color = 'darkgreen') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

doyvar <- ggplot(res, aes(doy, var)) +
  geom_smooth(color = 'dodgerblue3') +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

ggarrange(ymean, yvar, doymean, doyvar, 
          nrow = 2,
          ncol = 2)

#plot all ecozones

eco <- ecozones[st_intersection(ecozones, canada) %>% lengths > 0,]

eco <- ecozones[st_difference(st_union(ecozones),st_union(canada)) %>% lengths > 0,]

ggplot() +
  geom_sf(data = eco, aes(fill = ZONE_NAME)) +
  scale_fill_discrete(name = "Ecozones") + #limits = c(0, 1)) +
  geom_sf(data = canada, fill = NA, color = "black") + 
  geom_
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

#plot trends within/outside parks---------------------------------------------------

#generating predictions takes a while ://
#working - clean up formatting and add legends

#within parks
parkyes <- d[d$park == "TRUE",]

DAYS <- unique(d$date)
RES.parkyes <- list()

#generate predictions within parks
for(i in 1:length(DAYS)){
  newd <- filter(parkyes, date == DAYS[i])
  
  #add the predicted values back to the data frame with all the data
  preds <-
    bind_cols(newd,
              as.data.frame(predict.gam(m, newdata = newd, type = 'response')) %>%
                rename(mu = V1, #mean
                       phi = V2) %>% #scale
                mutate(sigma2 = phi * (1 - mu) * mu, #variance
                       mu = mu * 2 - 1, #rescale ndvi
                       sigma2 = sigma2 * 4)) #scale variance appropriately
  
  if(!is.na(preds$NDVI[1])){
    RES.parkyes[[i]] <- data.frame(day = unique(preds$date),
                                   mu = mean(preds$mu),
                                   var = mean(preds$sigma2))
  }
  
  
}

RESULTS.parkyes <- do.call(rbind, RES.parkyes)
RESULTS.parkyes <- mutate(RESULTS.parkyes,
                          dec_date = decimal_date(day),
                          year = year(day),
                          doy = yday(day))

#outside parks
parkno <- d[d$park == "FALSE",]
RES.parkno <- list()

#generate predictions outside parks
for(i in 1:length(DAYS)){
  newd <- filter(parkno, date == DAYS[i])
  
  #add the predicted values back to the data frame with all the data
  preds <-
    bind_cols(newd,
              as.data.frame(predict.gam(m, newdata = newd, type = 'response')) %>%
                rename(mu = V1, #mean
                       phi = V2) %>% #scale
                mutate(sigma2 = phi * (1 - mu) * mu, #variance
                       mu = mu * 2 - 1, #rescale ndvi
                       sigma2 = sigma2 * 4)) #scale variance appropriately
  
  if(!is.na(preds$NDVI[1])){
    RES.parkno[[i]] <- data.frame(day = unique(preds$date),
                                  mu = mean(preds$mu),
                                  var = mean(preds$sigma2))
  }
  
  
}

RESULTS.parkno <- do.call(rbind, RES.parkno)

RESULTS.parkyes <- mutate(RESULTS.parkyes,
                          dec_date = decimal_date(day),
                          year = year(day),
                          doy = yday(day))

RESULTS.parkno <- mutate(RESULTS.parkno,
                         dec_date = decimal_date(day),
                         year = year(day),
                         doy = yday(day))

parkymean <- ggplot() +
  geom_smooth(RESULTS.parkyes, mapping = aes(year, mu), color = 'darkgreen', fill = "#a7c957") +
  geom_smooth(RESULTS.parkno, mapping = aes(year, mu), color = '#a7c957', fill = "#ccd5ae") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        axis.text.x = element_text("Year"),
        axis.text.y = element_text("Mean"))

parkyvar <- ggplot() +
  geom_smooth(RESULTS.parkyes, mapping = aes(year, var), color = 'dodgerblue3', fill = "lightskyblue1") +
  geom_smooth(RESULTS.parkno, mapping = aes(year, var), color = 'lightskyblue1', fill = "#e0fbfc") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

parkdoymean <- ggplot() +
  geom_smooth(RESULTS.parkyes, mapping = aes(doy, mu), color = 'darkgreen', fill = "#a7c957") +
  geom_smooth(RESULTS.parkno, mapping = aes(doy, mu), color = '#a7c957', fill = "#ccd5ae") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

parkdoyvar <- ggplot() +
  geom_smooth(RESULTS.parkyes, mapping = aes(doy, var), color = 'dodgerblue3', fill = "lightskyblue1") +
  geom_smooth(RESULTS.parkno, mapping = aes(doy, var), color = 'lightskyblue1', fill = "#e0fbfc") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        legend.key.size = unit(1.5, "cm"))

ggarrange(parkymean, parkyvar, parkdoymean, parkdoyvar, 
          nrow = 2,
          ncol = 2, 
          labels = "auto")

ggsave("parktrends.png",
       width = ,
       height = ,
       dpi = ,
       bg = transparent)

#plot maps of mean & variance through time - by ecozone --------------------------

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

#map needs some cleaning up but works! add canada reference?
zone.mean.plot <- function(zone = 'ARC', year = 2010, doy = 100, res = 0.5){
  
  #create new coordinates to generate continuous predictions
  shp <- expand.grid(long = seq(ext(ecozones[ecozones$zone == zone,])[1], ext(ecozones[ecozones$zone == zone,])[2], by = res),
                     lat = seq(ext(ecozones[ecozones$zone == zone,])[3], ext(ecozones[ecozones$zone == zone,])[4], by = res))
  
  #convert to sf
  shp <- st_as_sf(shp, coords = c("long", "lat"), crs = "EPSG:4326")
  
  #crop coordinates to only ones within the given ecozone
  shp <- shp[ecozones[ecozones$zone == zone,], op = st_within]
  
  #convert back to dataframe
  shp.df <- shp %>%
    st_transform(4326) %>%
    st_coordinates() %>%
    as.data.frame()
  
  names(shp.df)[1] <- "long"
  names(shp.df)[2] <- "lat"
  
  shp.df$year <- year
  shp.df$doy <- doy
  shp.df$park <- as.factor(FALSE)
  shp.df$ecozone <- zone
  
  preds <- bind_cols(shp.df,
                     as.data.frame(predict(mm, newdata = shp.df, type = 'response')) %>%
                       rename(mu = 'predict(mm, newdata = shp.df, type = "response")') %>% # mean parameters
                       mutate(mu = mu * 2 - 1)) # rescale to [-1, 1])
  
  #plot the estimated mean and variance from the models onto a map 
  
  ggplot() +
    geom_raster(data = preds, aes(long, lat, fill = mu)) +
    scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.key.size = unit(1.5, "cm"))
  
}

meanplot <- zone.mean.plot(zone = 'MXP', year = 2010, doy = 100, res = 0.01) #SLAYED

#maps for variance
zone.var.plot <- function(zone = 'ARC', year = 2010, doy = 100, res = 0.5){
  
  #create new coordinates to generate continuous predictions
  shp <- expand.grid(long = seq(ext(ecozones[ecozones$zone == 'ARC',])[1], ext(ecozones[ecozones$zone == 'ARC',])[2], by = 0.01),
                     lat = seq(ext(ecozones[ecozones$zone == 'ARC',])[3], ext(ecozones[ecozones$zone == 'ARC',])[4], by = 0.01))
  
  #convert to sf
  shp <- st_as_sf(shp, coords = c("long", "lat"), crs = "EPSG:4326")
  
  #crop coordinates to only ones within the given ecozone
  shp <- shp[ecozones[ecozones$zone == 'ARC',], op = st_within]
  
  #convert back to dataframe
  shp.df <- shp %>%
    st_transform(4326) %>%
    st_coordinates() %>%
    as.data.frame()
  
  names(shp.df)[1] <- "long"
  names(shp.df)[2] <- "lat"
  
  shp.df$year <- 2010
  shp.df$doy <- 100
  shp.df$park <- as.factor(FALSE)
  shp.df$ecozone <- 'ARC'
  
  preds.m <- bind_cols(shp.df,
                       as.data.frame(predict(mm, newdata = shp.df, type = 'response')) %>%
                         rename(mu = 'predict(mm, newdata = shp.df, type = "response")')) 
  
  preds.mv <- bind_cols(preds.m,
                        as.data.frame(predict(mv, newdata = shp.df, type = 'response')) %>%
                          rename(phi = 'predict(mv, newdata = shp.df, type = "response")')) 
  
  preds <- mutate(preds.mv, 
                  sigma2 = phi * (1 - mu) * mu, #variance
                  mu = mu * 2 - 1, #rescale ndvi
                  var = sigma2 * 4) #scale variance appropriately
  
  #plot the estimated variance from the models onto a map 
  
  ggplot() +
    geom_raster(data = preds, aes(long, lat, fill = var)) +
    scale_fill_gradientn('Variance', colours = colors) + #limits = c(0, 1)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.key.size = unit(1.5, "cm"))
  
}

varplot <- zone.var.plot(zone = 'MXP', year = 2010, doy = 100, res = 0.01) #ATE

plot <- ggarrange(meanplot, varplot, ncol = 2,
                  labels = "auto")

annotate_figure(plot, top = text_grob("Mixedwood Plain", 
                                      color = "black", face = "bold", size = 25))

ggsave("mixedwoodplain.png",
       width = 13.86,
       height = 6.22,
       dpi = 300,
       bg = "transparent")

#plot maps of mean & variance through time - Canada -------------------------------

#function to easily plot different days and years
canada.plot <- function(year = 2010, doy = 100) {
  
  latlong <- d[,c(3, 8:10)]
  latlong <- distinct(latlong)
  
  latlong$year <- year
  latlong$doy <- doy
  latlong$park <- as.factor(FALSE)
  
  preds.m <- bind_cols(latlong,
                       as.data.frame(predict(mm, newdata = latlong, type = 'response')) %>%
                         rename(mu = 'predict(mm, newdata = latlong, type = "response")')) 
  
  preds.mv <- bind_cols(preds.m,
                        as.data.frame(predict(mv, newdata = latlong, type = 'response')) %>%
                          rename(phi = 'predict(mv, newdata = latlong, type = "response")')) 
  
  preds <- mutate(preds.mv, 
                  sigma2 = phi * (1 - mu) * mu, #variance
                  mu = mu * 2 - 1, #rescale ndvi
                  var = sigma2 * 4) #scale variance appropriately
  
  #plot the estimated mean and variance from the models onto a map 
  
  mean <- ggplot() +
    geom_raster(data = preds, aes(long, lat, fill = mu)) +
    geom_sf(data = canada, fill = NA, color = "black") +
    scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.key.size = unit(1.5, "cm"))
  
  var <- ggplot() +
    geom_raster(data = preds, aes(long, lat, fill = var)) +
    geom_sf(data = canada, fill = NA, color = "black") + 
    scale_fill_gradientn('Variance', colours = colors) + #limits = c(-1, 1)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.key.size = unit(1.5, "cm"))
  
  ggarrange(mean, var, ncol = 2)
  
}

can <- canada.plot(year = 2070, doy = 100)

annotate_figure(can, top = text_grob("April 2070", 
                                     color = "black", face = "bold", size = 25))

ggsave("canada2070.png",
       width = 13.86,
       height = 6.22,
       dpi = 300,
       bg = "transparent")
