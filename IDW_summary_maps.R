# Summary Maps Cancer Crabs

# Everett Rzeszowski
# April 22, 2025
# updated May 8, 2025

rm(list = ls())
DaisyDuck <- c("cowplot", "gstat", "ggpubr", "gratia", "FNN",
               "insurancerating",
               "marmap", "mgcv", "MuMIn", "patchwork",
               "regclass", "sf", "sp", "strucchange", 
               "terra", "tidyterra", "tidyverse", "wesanderson")
invisible(lapply(DaisyDuck, library, character.only = T))

# Shapefiles / spatial layers
NOAA_stat_areas <- st_read("~/Documents/R Projects/Shapefiles/NEFSC_GIS/Statistical_Areas_2010.shp") |>
  dplyr::filter(Id > 460 & Id < 641) 

load("~/Documents/R Projects/Shapefiles/NMFS_Lobster_StationsNLengths_1968_202404_UnAbridged_BigUnits_Detailed_250211.Rdata")
Configs$NMFS_Strata -> Strata
Strata_2 <- st_as_sf(Strata)

b2 = getNOAA.bathy(lon1 = -76.2, lon2 = -64, lat1 = 35, lat2 = 45, 
                   resolution = 1)
bf2 = tibble(fortify.bathy(b2))

countries <- map_data("world")
states <- map_data("state")

# Define stock areas
Jonah_crab_stocks_rough <- data.frame(Stock = c(rep("IGOM", times = 4), 
                                                rep("OGOM", times = 4), 
                                                rep("ISNE", times = 3), 
                                                rep("OSNE", times = 37), 
                                                rep("CAN", times = 10)),
                                      Id = c(511, 512, 513, 514,
                                             515, 521, 522, 561,
                                             538, 539, 611, 
                                             562, 525, 526, 537, 533, 534, 541, 542, 543, seq(612, 639, by = 1),
                                             463, 462, 468, 461, 467, 466, 465, 464, 552, 551)
)
Jonah_crab_stocks_rough |> left_join(NOAA_stat_areas) -> Jonah_crab_stocks_rough
NOAA_stat_areas |> left_join(Jonah_crab_stocks_rough) |> st_as_sf() -> NOAA_stat_areas
NOAA_stat_areas |> filter(! Id %in% c(623, 633, 468, 543, 542, 541, 534, 533, 624, 629, 639, 638, 634, 637, 628, 633, 640, 469)) -> good_habitat
NOAA_stat_areas |> filter(  Id %in% c(623, 633, 468, 543, 542, 541, 534, 533, 624, 629, 639, 638, 634, 637, 628, 633, 640, 469)) -> poor_habitat

sf_use_s2(FALSE)
NOAA_stat_areas |> group_by(Stock) |> #filter(Id != 631) |>
  summarize(geometry = sf::st_union(geometry)) |>
  ungroup() -> Stocks

Stocks2 <- st_read("~/Documents/R Projects/Cancer_HSI2/NOAA_Edit_Shape/editShape.shp")
  

# springJonah_rast <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringJonah.tif")
# fallJonah_rast <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallJonah.tif")
# springRock_rast <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringRock.tif")
# fallRock_rast <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallRock.tif")

springJonah_rast_crop <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringJonah_crop.tif")
fallJonah_rast_crop <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallJonah_crop.tif")
springRock_rast_crop <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringRock_crop.tif")
fallRock_rast_crop <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallRock_crop.tif")

# remove areas with no tows to best of ability
# springJonah_rast <- terra::mask(springJonah_rast, good_habitat)
# fallJonah_rast <- terra::mask(fallJonah_rast, good_habitat)
# springRock_rast <- terra::mask(springRock_rast, good_habitat)
# fallRock_rast <- terra::mask(fallRock_rast, good_habitat)

# springJonah_rast_crop <- terra::mask(springJonah_rast_crop, good_habitat)
# fallJonah_rast_crop <- terra::mask(fallJonah_rast_crop, good_habitat)
# springRock_rast_crop <- terra::mask(springRock_rast_crop, good_habitat)
# fallRock_rast_crop <- terra::mask(fallRock_rast_crop, good_habitat)

springJonah_rast_crop <- terra::mask(springJonah_rast_crop, Strata_2)
fallJonah_rast_crop <- terra::mask(fallJonah_rast_crop, Strata_2)
springRock_rast_crop <- terra::mask(springRock_rast_crop, Strata_2)
fallRock_rast_crop <- terra::mask(fallRock_rast_crop, Strata_2)

# u_sj_rast <- terra::app(springJonah_rast, fun = mean, na.rm = F)
# u_fj_rast <- terra::app(fallJonah_rast, fun = mean, na.rm = F)
# u_sr_rast <- terra::app(springRock_rast, fun = mean, na.rm = F)
# u_fr_rast <- terra::app(fallRock_rast, fun = mean, na.rm = F)

u_sj_rast_crop <- terra::app(springJonah_rast_crop, fun = mean, na.rm = F)
u_fj_rast_crop <- terra::app(fallJonah_rast_crop, fun = mean, na.rm = F)
u_sr_rast_crop <- terra::app(springRock_rast_crop, fun = mean, na.rm = F)
u_fr_rast_crop <- terra::app(fallRock_rast_crop, fun = mean, na.rm = F)

Stock_Label <- tibble(Stock = c("CAN", "IGOM", "ISNE", "OGOM", "OSNE"), 
                      lon = c(-66.5, -70, -72, -69, -72),
                      lat = c(43, 43.2, 41.2, 42.2, 39.7))

# plotting stocks there are issues in the shape file
  # issues are somewhere in this 40000 row shitshow

############## FIGURE 2 UNCROPPED MEAN SEASONAL SURVEY DISTRIBUTION ############
# jonah spring
ggplot() +
  geom_raster(data = u_sj_rast_crop, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 0.1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_point(data = springJonah_trawl, aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT)) +
  geom_sf(data = Stonks |> filter(is.na(Stock)), fill = "grey", alpha = 0.85) + 
  geom_sf(data = Stonks, aes(group = Stock, color = Stock), alpha = 0.01, linewidth = 0.5) + 
  geom_sf(data = poor_habitat, fill = "grey", alpha = 0.85) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Jonah crab, mean Spring abundance",
       y = "Latitude",
       x = "",
       fill = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 18) -> sj1

ggsave(filename = "fig2testpoints.png",
       plot = sj1,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 250,
       height = 7,
       width = 10)

# jonah fall
ggplot() +
  geom_raster(data = u_fj_rast, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 0.1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stonks |> filter(is.na(Stock)), fill = "grey", alpha = 0.75) + 
  geom_sf(data = Stonks, aes(group = Stock, color = Stock), alpha = 0.01) + 
  geom_sf(data = poor_habitat, fill = "grey", alpha = 0.85) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Jonah crab, mean Fall abundance",
       y = "Latitude",
       x = "",
       fill = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 18) -> fj1

# rock spring
ggplot() +
  geom_raster(data = u_sr_rast, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 0.1), na.value = "#9e0142") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stonks |> filter(is.na(Stock)), fill = "grey", alpha = 0.75) + 
  geom_sf(data = Stonks, aes(group = Stock, color = Stock), alpha = 0.01) + 
  geom_sf(data = poor_habitat, fill = "grey", alpha = 0.85) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Rock crab, mean Spring abundance",
       y = "Latitude",
       x = "Longitude",
       fill = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 18) -> sr1

ggsave(filename = "fig2testSpringRock.png",
       plot = sr1,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 250,
       height = 7,
       width = 10)

# rock fall
ggplot() +
  geom_raster(data = u_fr_rast, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 0.1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stonks |> filter(is.na(Stock)), fill = "grey", alpha = 0.75) + 
  geom_sf(data = Stonks, aes(group = Stock, color = Stock), alpha = 0.01) + 
  geom_sf(data = poor_habitat, fill = "grey", alpha = 0.85) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Rock crab, mean Fall abundance",
       y = "",
       x = "Longitude",
       fill = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 18) -> fr1

(sj1 + fj1) / (sr1 + fr1) + plot_annotation(tag_levels = "A") -> f2

ggsave(filename = "figure2.png",
       plot = f2,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 14,
       width = 20)

######################## CROPPED DISTRIBUTIONAL MAPS #########################
# jonah spring
ggplot() +
  geom_raster(data = u_sj_rast_crop, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 0.1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = NOAA_stat_areas, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
 # geom_sf(data = Stocks |> filter(is.na(Stock)), fill = "grey", alpha = 0.85) + 
 # geom_sf(data = Stonks, alpha = 0.1, linewidth = 0.5) +
 # geom_sf(data = poor_habitat, fill = "grey", alpha = 0.85) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Jonah crab, mean Spring abundance",
       y = "Latitude",
       x = "",
       fill = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 18) -> sj1_crop

ggsave(filename = "fig2test_crop.png",
       plot = sj1_crop,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 250,
       height = 7,
       width = 10)

# jonah fall
ggplot() +
  geom_raster(data = u_fj_rast_crop, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 0.1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = NOAA_stat_areas, aes(group = Stock, color = Stock), alpha = 0.01) + 
 # geom_sf(data = Stonks |> filter(is.na(Stock)), fill = "grey", alpha = 0.85) + 
 # geom_sf(data = poor_habitat, fill = "grey", alpha = 0.85) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Jonah crab, mean Fall abundance",
       y = "Latitude",
       x = "",
       fill = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 18) -> fj1_crop

# rock spring
ggplot() +
  geom_raster(data = u_sr_rast_crop, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 0.15), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = NOAA_stat_areas, aes(group = Stock, color = Stock), alpha = 0.01) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Rock crab, mean Spring abundance",
       y = "Latitude",
       x = "Longitude",
       fill = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 18) -> sr1_crop

ggsave(filename = "fig2testSpringRockCrop.png",
       plot = sr1_crop,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 250,
       height = 7,
       width = 10)

# rock fall
ggplot() +
  geom_raster(data = u_fr_rast_crop, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 0.1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = NOAA_stat_areas, aes(group = Stock, color = Stock), alpha = 0.01) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Rock crab, mean Fall abundance",
       y = "",
       x = "Longitude",
       fill = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 18) -> fr1_crop

(sj1_crop + fj1_crop) / (sr1_crop + fr1_crop) + plot_annotation(tag_levels = "A") -> f2_crop

ggsave(filename = "figure2_crop.png",
       plot = f2_crop,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 14,
       width = 20)

##############
# Calculate centroid annual population centroid location by weighting interpolated point
Can  <- st_as_sf(Stocks |> filter(Stock == "CAN"))
Isne <- st_as_sf(Stocks |> filter(Stock == "ISNE"))
Osne <- st_as_sf(Stocks |> filter(Stock == "OSNE"))
Igom <- st_as_sf(Stocks |> filter(Stock == "IGOM"))
Ogom <- st_as_sf(Stocks |> filter(Stock == "OGOM"))

# Calculate annual mean Catch per tow by population
cut.sjCan_crop        <- mask(springJonah_rast_crop, Can)
cut.sjIsne_crop       <- mask(springJonah_rast_crop, Isne)
cut.sjOsne_crop       <- mask(springJonah_rast_crop, Osne)
cut.sjIgom_crop       <- mask(springJonah_rast_crop, Igom)
cut.sjOgom_crop       <- mask(springJonah_rast_crop, Ogom)

cut.fjCan_crop        <- mask(fallJonah_rast_crop, Can)
cut.fjIsne_crop       <- mask(fallJonah_rast_crop, Isne)
cut.fjOsne_crop       <- mask(fallJonah_rast_crop, Osne)
cut.fjIgom_crop       <- mask(fallJonah_rast_crop, Igom)
cut.fjOgom_crop       <- mask(fallJonah_rast_crop, Ogom)

cut.srCan_crop        <- mask(springRock_rast_crop, Can)
cut.srIsne_crop       <- mask(springRock_rast_crop, Isne)
cut.srOsne_crop       <- mask(springRock_rast_crop, Osne)
cut.srIgom_crop       <- mask(springRock_rast_crop, Igom)
cut.srOgom_crop       <- mask(springRock_rast_crop, Ogom)

cut.frCan_crop        <- mask(fallRock_rast_crop, Can)
cut.frIsne_crop       <- mask(fallRock_rast_crop, Isne)
cut.frOsne_crop       <- mask(fallRock_rast_crop, Osne)
cut.frIgom_crop       <- mask(fallRock_rast_crop, Igom)
cut.frOgom_crop       <- mask(fallRock_rast_crop, Ogom)
########### Calucalte weighted spatial averages (weighted centroids) for populations
COBiomass <- function(r){
  biomass_center <- vect(data.frame(x = numeric(0), y = numeric(0)),
                         geom = c("x", "y"),
                         crs = crs(r))
  
  for(i in 1:terra::nlyr(r)){
    print(i)
    # Get the coordinates and values
    values <- values(r)[,i]  # abundance (e.g., var1.pred)
    coords <- crds(r)        # cell center coordinates
    
    # Remove NA values; Coords only returns non-NA
    valid <- !is.na(values)
    vals <- values[valid]
    
    # Weighted mean (center of biomass)
    x_center <- sum(coords[,1] * vals) / sum(vals)
    y_center <- sum(coords[,2] * vals) / sum(vals)
    temp <- vect(data.frame(x = x_center, y = y_center), 
                            geom = c("x", "y"),
                            crs = crs(r))
    # Create a point SpatVector for the biomass center
    biomass_center <- rbind(biomass_center, temp)
  }
  
  return(biomass_center)
  }
  
sjCAN_centroids  <- COBiomass(cut.sjCan_crop)
sjIGOM_centroids <- COBiomass(cut.sjIgom_crop)
sjOGOM_centroids <- COBiomass(cut.sjOgom_crop)
sjISNE_centroids <- COBiomass(cut.sjIsne_crop)
sjOSNE_centroids <- COBiomass(cut.sjOsne_crop)
sj_centroids <- COBiomass(springJonah_rast_crop)

fjCAN_centroids  <- COBiomass(cut.fjCan_crop)
fjIGOM_centroids <- COBiomass(cut.fjIgom_crop)
fjOGOM_centroids <- COBiomass(cut.fjOgom_crop)
fjISNE_centroids <- COBiomass(cut.fjIsne_crop)
fjOSNE_centroids <- COBiomass(cut.fjOsne_crop)
fj_centroids <- COBiomass(fallJonah_rast_crop)

srCAN_centroids  <- COBiomass(cut.srCan_crop)
srIGOM_centroids <- COBiomass(cut.srIgom_crop)
srOGOM_centroids <- COBiomass(cut.srOgom_crop)
srISNE_centroids <- COBiomass(cut.srIsne_crop)
srOSNE_centroids <- COBiomass(cut.srOsne_crop)
sr_centroids <- COBiomass(springRock_rast_crop)

frCAN_centroids  <- COBiomass(cut.frCan_crop)
frIGOM_centroids <- COBiomass(cut.frIgom_crop)
frOGOM_centroids <- COBiomass(cut.frOgom_crop)
frISNE_centroids <- COBiomass(cut.frIsne_crop)
frOSNE_centroids <- COBiomass(cut.frOsne_crop)
fr_centroids <- COBiomass(fallRock_rast_crop)

sj_Allcentroids <- bind_rows(as.data.frame(geom(sjCAN_centroids)), 
                             as.data.frame(geom(sjIGOM_centroids)), 
                             as.data.frame(geom(sjOGOM_centroids)), 
                             as.data.frame(geom(sjISNE_centroids)), 
                             as.data.frame(geom(sjOSNE_centroids)),
                             as.data.frame(geom(sj_centroids)),
                             .id = "Stock") |> 
  mutate(Stock = factor(Stock, labels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE", "Total")),
         Year = rep(seq(1968, 2024, by = 1), times = 6)) |>
  filter(!is.nan(x))

fj_Allcentroids <- bind_rows(as.data.frame(geom(fjCAN_centroids)), 
                             as.data.frame(geom(fjIGOM_centroids)), 
                             as.data.frame(geom(fjOGOM_centroids)), 
                             as.data.frame(geom(fjISNE_centroids)), 
                             as.data.frame(geom(fjOSNE_centroids)),
                             as.data.frame(geom(fj_centroids)),
                             .id = "Stock") |> 
  mutate(Stock = factor(Stock, labels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE", "Total")),
         Year = rep(c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)), times = 6)) |>
  filter(!is.nan(x))

sr_Allcentroids <- bind_rows(as.data.frame(geom(srCAN_centroids)), 
                             as.data.frame(geom(srIGOM_centroids)), 
                             as.data.frame(geom(srOGOM_centroids)), 
                             as.data.frame(geom(srISNE_centroids)), 
                             as.data.frame(geom(srOSNE_centroids)),
                             as.data.frame(geom(sr_centroids)),
                             .id = "Stock") |> 
  mutate(Stock = factor(Stock, labels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE", "Total")),
         Year = rep(seq(1968, 2024, by = 1), times = 6)) |>
  filter(!is.nan(x))

fr_Allcentroids <- bind_rows(as.data.frame(geom(frCAN_centroids)), 
                             as.data.frame(geom(frIGOM_centroids)), 
                             as.data.frame(geom(frOGOM_centroids)), 
                             as.data.frame(geom(frISNE_centroids)), 
                             as.data.frame(geom(frOSNE_centroids)),
                             as.data.frame(geom(fr_centroids)),
                             .id = "Stock") |> 
  mutate(Stock = factor(Stock, labels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE", "Total")),
         Year = rep(c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)), times = 6)) |>
  filter(!is.nan(x))

############### Associate Marmap depths with centroid points ###################
sjnn <- get.knnx(bf2 |> select(x, y), 
                 sj_Allcentroids |> filter(!is.nan(x)) |> select(x, y), k = 1)

sj_Allcentroids <- sj_Allcentroids %>%
  mutate(depth = bf2$z[sjnn$nn.index[,1]])

fjnn <- get.knnx(bf2 |> select(x, y), 
                 fj_Allcentroids |> filter(!is.nan(x)) |> select(x, y), k = 1)

fj_Allcentroids <- fj_Allcentroids %>%
  mutate(depth = bf2$z[fjnn$nn.index[,1]])

srnn <- get.knnx(bf2 |> select(x, y), 
                 sr_Allcentroids |> filter(!is.nan(x)) |> select(x, y), k = 1)

sr_Allcentroids <- sr_Allcentroids %>%
  mutate(depth = bf2$z[srnn$nn.index[,1]])

frnn <- get.knnx(bf2 |> select(x, y), 
                 fr_Allcentroids |> filter(!is.nan(x)) |> select(x, y), k = 1)

fr_Allcentroids <- fr_Allcentroids %>%
  mutate(depth = bf2$z[frnn$nn.index[,1]])

##################### Visualize center-of-biomass on maps ######################
ggplot() +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -100, -2000), color = "gray79") +
  geom_sf(data = NOAA_stat_areas, aes(group = Stock, fill = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_point(data = sj_Allcentroids |> filter(Stock == "Total" & Year != 2020 & Year != 2023),
             aes(x = x, y = y, color = as.numeric(Year)), size = 2.5) +
  #scale_color_distiller(palette = "Spectral",  limits = c(1968, 2024), na.value = "white") +
  scale_color_gradientn(colors = wes_palette(name = "Zissou1", type = "continuous"),  limits = c(1968, 2024), na.value = "white") +
  labs(
       subtitle = "Spring Jonah Crab",
       x = "Longitude",
       y = "Latitude",
       color = "Year") +
  theme_classic(base_size = 18) -> sj_centroidMap

  ggplot() +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -100, -2000), color = "gray79") +
    geom_sf(data = NOAA_stat_areas, aes(group = Stock, fill = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_point(data = fj_Allcentroids |> filter(Stock == "Total" & Year != 2020 & Year != 2023),
               aes(x = x, y = y, color = as.numeric(Year)), size = 2.5) +
    #scale_color_distiller(palette = "Spectral",  limits = c(1968, 2024), na.value = "white") +
    scale_color_gradientn(colors = wes_palette(name = "Zissou1", type = "continuous"),  limits = c(1968, 2024), na.value = "white") +
    labs(subtitle = "Fall Jonah Crab",
         x = "Longitude", 
         y = "Latitude",
         color = "Year") +
  theme_classic(base_size = 18) -> fj_centroidMap
  
  ggplot() +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -100, -2000), color = "gray79") +
    geom_sf(data = NOAA_stat_areas, aes(group = Stock, fill = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_point(data = sr_Allcentroids |> filter(Stock == "Total" & Year != 2020 & Year != 2023),
               aes(x = x, y = y, color = as.numeric(Year)), size = 2.5) +
    #scale_color_distiller(palette = "Spectral",  limits = c(1968, 2024), na.value = "white") +
    scale_color_gradientn(colors = wes_palette(name = "Zissou1", type = "continuous"),  limits = c(1968, 2024), na.value = "white") +
    labs(subtitle = "Spring Atl. Rock Crab",
         x = "Longitude", 
         y = "Latitude",
         color = "Year") +
  theme_classic(base_size = 18) -> sr_centroidMap
  
  ggplot() +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -100, -2000), color = "gray79") +
    geom_sf(data = NOAA_stat_areas, aes(group = Stock, fill = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_point(data = fr_Allcentroids |> filter(Stock == "Total" & Year != 2020 & Year != 2023),
               aes(x = x, y = y, color = as.numeric(Year)), size = 2.5) +
    #scale_color_distiller(palette = "Spectral",  limits = c(1968, 2024), na.value = "white") +
    scale_color_gradientn(colors = wes_palette(name = "Zissou1", type = "continuous"),  limits = c(1968, 2024), na.value = "white") +
    labs(subtitle = "Fall Atl. Rock Crab",
         x = "Longitude", 
         y = "Latitude",
         color = "Year") +
  theme_classic(base_size = 18) -> fr_centroidMap
  
  sj_Allcentroids |> filter(Year != 2020 & Year != 2023 & Stock == "Total" & depth > -500) |> 
    ggplot(aes(x = Year, y = y)) +
    geom_line(linewidth = 2) + 
    geom_smooth(method = "gam") +
    coord_cartesian(ylim = c(37, 43)) +
    labs(title = "Center-of-Biomass",
         x = "Year",
         y = "Latitude (°N)") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none") -> sj_COB_ts
  
  fj_Allcentroids |> filter(Year != 2020 & Year != 2023 & Stock == "Total") |> 
    ggplot(aes(x = Year, y = y)) +
    geom_line(linewidth = 2) + 
    geom_smooth(method = "gam") +
    coord_cartesian(ylim = c(37, 43)) +
    labs(title = "Center-of-Biomass",
         x = "Year",
         y = "Latitude (°N)") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none") -> fj_COB_ts
  
  sr_Allcentroids |> filter(Year != 2020 & Year != 2023 & Stock == "Total" ) |> 
    ggplot(aes(x = Year, y = y)) +
    geom_line(linewidth = 2) + 
    geom_smooth(method = "gam") +
    coord_cartesian(ylim = c(37, 43)) +
    labs(title = "Center-of-Biomass",
         x = "Year",
         y = "Latitude (°N)") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none") -> sr_COB_ts
  
  fr_Allcentroids |> filter(Year != 2020 & Year != 2023 & Stock == "Total" ) |> 
    ggplot(aes(x = Year, y = y)) +
    geom_line(linewidth = 2) + 
    geom_smooth(method = "gam") +
    coord_cartesian(ylim = c(37, 43)) +
    labs(title = "Center-of-Biomass",
         x = "Year",
         y = "Latitude (°N)") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none") -> fr_COB_ts

(sj_centroidMap + fj_centroidMap) / (sr_centroidMap + fr_centroidMap) + 
    plot_layout(guides = "collect") +
    plot_annotation(title = "Annual Species Center-of-Biomass",
                    tag_levels = "A", 
                    theme = theme(title = element_text(size = 17))) -> COB_maps
  
ggdraw() +
  draw_plot(COB_maps) +
  draw_plot(sj_COB_ts, x = 0.13, y = 0.76, width = 0.12, height  = 0.14) +
  draw_plot(sr_COB_ts, x = 0.13, y = 0.28, width = 0.12, height  = 0.14) +
  draw_plot(fj_COB_ts, x = 0.55, y = 0.76, width = 0.12, height  = 0.14) +
  draw_plot(fr_COB_ts, x = 0.55, y = 0.28, width = 0.12, height  = 0.14) -> COB_map2
    
ggsave(filename = "COB_map.png",
       plot = COB_map2,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 150,
       height = 14,
       width = 16)

############ Center of biomass time series ############
# (sj_COB_ts + fj_COB_ts) / (sr_COB_ts + fr_COB_ts) +
#   plot_annotation(title = "Seasonal Centers-of-Biomass", 
#                   tag_levels = "A", 
#                   theme = theme(title = element_text(size = 17))) -> COB_ts
# 
# ggsave(filename = "COB_lat_ts.png",
#        plot = COB_ts,
#        path = "~/Documents/R Projects/Cancer_HSI2/Figures",
#        dpi = 750,
#        height = 7,
#        width = 10)

################# Center-of-Biomass (depth) ###########
sj_Allcentroids |> filter(Year != 2020 & Year != 2023 & Stock == "Total" & depth > -500) |> 
  ggplot(aes(x = Year, y = depth)) +
  #geom_rect(data = sj_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  #scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2) + 
  geom_smooth(method = "gam") +
  #coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Spring Jonah Crab",
       x = "Year",
       y = "Center-of-Biomass (m)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

fj_Allcentroids |> filter(Year != 2020 & Year != 2023 & Stock == "Total") |> 
  ggplot(aes(x = Year, y = depth)) +
  #geom_rect(data = sj_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  #scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2) + 
  geom_smooth(method = "gam") +
  #coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Fall Jonah Crab",
       x = "Year",
       y = "Center-of-Biomass (m)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

sr_Allcentroids |> filter(Year != 2020 & Year != 2023 & Stock == "Total" ) |> 
  ggplot(aes(x = Year, y = depth)) +
  #geom_rect(data = sj_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  #scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2) + 
  geom_smooth(method = "gam") +
  #coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Spring Rock Crab",
       x = "Year",
       y = "Center-of-Biomass (m)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

fr_Allcentroids |> filter(Year != 2020 & Year != 2023 & Stock == "Total" ) |> 
  ggplot(aes(x = Year, y = depth)) +
  #geom_rect(data = sj_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  #scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2) + 
  geom_smooth(method = "gam") +
  #coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Fall Rock Crab",
       x = "Year",
       y = "Center-of-Biomass (m)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

########### TIME SERIES uncropped ################
########### Spring jonah ##############
# Calculate annual mean Catch per tow by population
cut.sjCan       <- mask(springJonah_rast_crop, Can)
cut.sjIsne      <- mask(springJonah_rast_crop, Isne)
cut.sjOsne      <- mask(springJonah_rast_crop, Osne)
cut.sjIgom       <- mask(springJonah_rast_crop, Igom)
cut.sjOgom       <- mask(springJonah_rast_crop, Ogom)

mean.cut.sjCan       <- terra::global(cut.sjCan, fun = mean, na.rm = T)
sjCan2             <- tibble(mean = mean.cut.sjCan$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("CAN"))

mean.cut.sjIsne       <- terra::global(cut.sjIsne, fun = mean, na.rm = T)
sjIsne2            <- tibble(mean = mean.cut.sjIsne$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("ISNE"))

mean.cut.sjOsne       <- terra::global(cut.sjOsne, fun = mean, na.rm = T)
sjOsne2            <- tibble(mean = mean.cut.sjOsne$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("OSNE"))

mean.cut.sjIgom       <- terra::global(cut.sjIgom, fun = mean, na.rm = T)
sjIgom2            <- tibble(mean = mean.cut.sjIgom$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("IGOM"))

mean.cut.sjOgom       <- terra::global(cut.sjOgom, fun = mean, na.rm = T)
sjOgom2            <- tibble(mean = mean.cut.sjOgom$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("OGOM"))

# Changepoint analysis on each timeseries
CANsjts <- sjCan2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjCANbp <- breakpoints(ts(CANsjts)~1, breaks = 6)

ISNEsjts <- sjIsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjISNEbp <- breakpoints(ts(ISNEsjts)~1, breaks = 6)

OSNEsjts <- sjOsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjOSNEbp <- breakpoints(ts(OSNEsjts)~1, breaks = 6)

IGOMsjts <- sjIgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjIGOMbp <- breakpoints(ts(IGOMsjts)~1, breaks = 6)

OGOMsjts <- sjOgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjOGOMbp <- breakpoints(ts(OGOMsjts)~1, breaks = 6)
summary(sjOGOMbp) 

# combine results 
bind_rows(sjCan2, sjIsne2, sjOsne2, sjIgom2, sjOgom2) |>
  mutate(Stock = factor(Stock, levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE"))) -> sjTSagg

summary(sjCANbp)#$RSS[2,]
summary(sjISNEbp) 
summary(sjOSNEbp) 
summary(sjIGOMbp) 
summary(sjOGOMbp) 

sj_rect <- tibble(xm = c(1968, sjCan2$year[31], sjCan2$year[45], 1968, sjIsne2$year[13], 1968, sjOsne2$year[14],
                         1968, sjIgom2$year[29], sjIgom2$year[45], 1968, sjOgom2$year[45]),
                  xM = c(sjCan2$year[31], sjCan2$year[45], 2024, sjIsne2$year[13], 2024, sjOsne2$year[14], 2024,
                         sjIgom2$year[29], sjIgom2$year[45], 2024, sjOgom2$year[45], 2024),
                  ym = rep(0),
                  yM = rep(0.25),
                  regime = c(1, 2, 3, 1, 2, 1, 2, 
                             1, 2, 3, 1, 2),
                  Stock = factor(c("CAN", "CAN", "CAN", "ISNE", "ISNE", "OSNE", "OSNE","IGOM", "IGOM", "IGOM", "OGOM", "OGOM"), 
                            levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE")))

sjTSagg |> filter(year != 2020 & year != 2023) |> 
  ggplot() +
  facet_grid(.~Stock) +
  geom_rect(data = sj_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2, aes(x = year, y = mean)) + 
  geom_smooth(method = "gam", aes(x = year, y = mean)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Spring Jonah Crab",
       x = "Year",
       y = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") -> sjPlot


########## Fall Jonah #############
cut.fjCan       <- mask(fallJonah_rast_crop, Can)
cut.fjIsne      <- mask(fallJonah_rast_crop, Isne)
cut.fjOsne      <- mask(fallJonah_rast_crop, Osne)
cut.fjIgom       <- mask(fallJonah_rast_crop, Igom)
cut.fjOgom       <- mask(fallJonah_rast_crop, Ogom)

mean.cut.fjCan       <- terra::global(cut.fjCan, fun = mean, na.rm = T)
fjCan2             <- tibble(mean = mean.cut.fjCan$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("CAN"))

mean.cut.fjIsne       <- terra::global(cut.fjIsne, fun = mean, na.rm = T)
fjIsne2            <- tibble(mean = mean.cut.fjIsne$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("ISNE"))

mean.cut.fjOsne       <- terra::global(cut.fjOsne, fun = mean, na.rm = T)
fjOsne2            <- tibble(mean = mean.cut.fjOsne$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("OSNE"))

mean.cut.fjIgom       <- terra::global(cut.fjIgom, fun = mean, na.rm = T)
fjIgom2            <- tibble(mean = mean.cut.fjIgom$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("IGOM"))

mean.cut.fjOgom       <- terra::global(cut.fjOgom, fun = mean, na.rm = T)
fjOgom2            <- tibble(mean = mean.cut.fjOgom$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("OGOM"))

# Changepoint analysis on each timeseries
CANfjts <- fjCan2 |> filter(year > 1968) |> select(mean)
fjCANbp <- breakpoints(ts(CANfjts)~1, breaks = 6)
summary(fjCANbp) 

ISNEfjts <- fjIsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
fjISNEbp <- breakpoints(ts(ISNEfjts)~1, breaks = 6)
summary(fjISNEbp) 

OSNEfjts <- fjOsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
fjOSNEbp <- breakpoints(ts(OSNEfjts)~1, breaks = 6)
summary(fjOSNEbp) 

IGOMfjts <- fjIgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
fjIGOMbp <- breakpoints(ts(IGOMfjts)~1, breaks = 6)
summary(fjIGOMbp) 

OGOMfjts <- fjOgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
fjOGOMbp <- breakpoints(ts(OGOMfjts)~1, breaks = 6)
summary(fjOGOMbp) 

# combine results 
bind_rows(fjCan2, fjIsne2, fjOsne2, fjIgom2, fjOgom2) |>
  mutate(Stock = factor(Stock, levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE"))) -> fjTSagg

summary(fjCANbp)#$RSS[2,]
summary(fjISNEbp) 
summary(fjOSNEbp) 
summary(fjIGOMbp) 
summary(fjOGOMbp) 

fj_rect <- tibble(xm = c(1968, fjCan2$year[30], 1968, fjIsne2$year[42], 1968, fjOsne2$year[40],
                         1968, fjIgom2$year[30], fjIgom2$year[43], 1968, fjOgom2$year[30]),
                  xM = c(fjCan2$year[30], 2024, fjIsne2$year[42], 2024, fjOsne2$year[40], 2024,
                         fjIgom2$year[30], fjIgom2$year[43], 2024, fjOgom2$year[30], 2024),
                  ym = rep(0),
                  yM = rep(0.25),
                  regime = c(1, 2, 1, 2, 1, 2, 
                             1, 2, 3, 1, 2),
                  Stock = factor(c("CAN", "CAN", "ISNE", "ISNE", "OSNE", "OSNE","IGOM", "IGOM", "IGOM", "OGOM", "OGOM"), 
                                 levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE")))
fjTSagg |> #filter(year != 2020 & year != 2023) |> 
  ggplot() +
  facet_grid(.~Stock) +
  geom_rect(data = fj_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2, aes(x = year, y = mean)) + 
  geom_smooth(method = "gam", aes(x = year, y = mean)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Fall Jonah Crab",
       x = "Year",
       y = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") -> fjPlot

########### Spring Rock ##############
# Calculate annual mean Catch per tow by population
cut.srCan       <- mask(springRock_rast_crop, Can)
cut.srIsne      <- mask(springRock_rast_crop, Isne)
cut.srOsne      <- mask(springRock_rast_crop, Osne)
cut.srIgom      <- mask(springRock_rast_crop, Igom)
cut.srOgom      <- mask(springRock_rast_crop, Ogom)

mean.cut.srCan       <- terra::global(cut.srCan, fun = mean, na.rm = T)
srCan2             <- tibble(mean = mean.cut.srCan$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("CAN"))

mean.cut.srIsne       <- terra::global(cut.srIsne, fun = mean, na.rm = T)
srIsne2            <- tibble(mean = mean.cut.srIsne$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("ISNE"))

mean.cut.srOsne       <- terra::global(cut.srOsne, fun = mean, na.rm = T)
srOsne2            <- tibble(mean = mean.cut.srOsne$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("OSNE"))

mean.cut.srIgom       <- terra::global(cut.srIgom, fun = mean, na.rm = T)
srIgom2            <- tibble(mean = mean.cut.srIgom$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("IGOM"))

mean.cut.srOgom       <- terra::global(cut.srOgom, fun = mean, na.rm = T)
srOgom2            <- tibble(mean = mean.cut.srOgom$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("OGOM"))

# Changepoint analysis on each timeseries
CANsrts <- srCan2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srCANbp <- breakpoints(ts(CANsrts)~1, breaks = 6)
summary(srCANbp) 

ISNEsrts <- srIsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srISNEbp <- breakpoints(ts(ISNEsrts)~1, breaks = 6)
summary(srISNEbp) 

OSNEsrts <- srOsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srOSNEbp <- breakpoints(ts(OSNEsrts)~1, breaks = 6)
summary(srOSNEbp) 

IGOMsrts <- srIgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srIGOMbp <- breakpoints(ts(IGOMsrts)~1, breaks = 6)
summary(srIGOMbp) 

OGOMsrts <- srOgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srOGOMbp <- breakpoints(ts(OGOMsrts)~1, breaks = 6)
summary(srOGOMbp) 

# combine results 
bind_rows(srCan2, srIsne2, srOsne2, srIgom2, srOgom2) |>
  mutate(Stock = factor(Stock, levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE"))) -> srTSagg


summary(srCANbp)$RSS[2,]
summary(srISNEbp) 
summary(srOSNEbp) 
summary(srIGOMbp) $RSS[2,]
summary(srOGOMbp) 

sr_rect <- tibble(xm = c(1968, 1968, srIsne2$year[8], srIsne2$year[16], 1968, srOsne2$year[8], srOsne2$year[16],
                         1968, 1968, sjOgom2$year[39]),
                  xM = c(2024, sjIsne2$year[8], sjIsne2$year[16], 2024, srOsne2$year[8], srOsne2$year[16], 2024,
                          2024, sjOgom2$year[39], 2024),
                  ym = rep(0),
                  yM = rep(0.25),
                  regime = c(1, 1, 2, 3, 1, 2, 3,  
                             1, 1, 2),
                  Stock = factor(c("CAN", "ISNE", "ISNE", "ISNE", "OSNE", "OSNE", "OSNE", "IGOM", "OGOM", "OGOM"), 
                                 levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE")))


srTSagg |> filter(year != 2020 & year != 2023) |> 
  ggplot() +
  facet_grid(.~Stock) +
  geom_rect(data = sr_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2, aes(x = year, y = mean)) + 
  geom_smooth(method = "gam", aes(x = year, y = mean)) +
  
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Spring Atlantic Rock Crab",
       x = "Year",
       y = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") ->srPlot


########## Fall Rock #############
cut.frCan       <- mask(fallRock_rast_crop, Can)
cut.frIsne      <- mask(fallRock_rast_crop, Isne)
cut.frOsne      <- mask(fallRock_rast_crop, Osne)
cut.frIgom      <- mask(fallRock_rast_crop, Igom)
cut.frOgom      <- mask(fallRock_rast_crop, Ogom)

mean.cut.frCan       <- terra::global(cut.frCan, fun = mean, na.rm = T)
frCan2             <- tibble(mean = mean.cut.frCan$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("CAN"))

mean.cut.frIsne       <- terra::global(cut.frIsne, fun = mean, na.rm = T)
frIsne2            <- tibble(mean = mean.cut.frIsne$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("ISNE"))

mean.cut.frOsne       <- terra::global(cut.frOsne, fun = mean, na.rm = T)
frOsne2            <- tibble(mean = mean.cut.frOsne$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("OSNE"))

mean.cut.frIgom       <- terra::global(cut.frIgom, fun = mean, na.rm = T)
frIgom2            <- tibble(mean = mean.cut.frIgom$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("IGOM"))

mean.cut.frOgom       <- terra::global(cut.frOgom, fun = mean, na.rm = T)
frOgom2            <- tibble(mean = mean.cut.frOgom$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("OGOM"))

# Changepoint analysis on each timeseries
CANfrts <- frCan2 |> filter(year > 1968) |> select(mean)
frCANbp <- breakpoints(ts(CANfrts)~1, breaks = 6)
summary(frCANbp) 

ISNEfrts <- frIsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
frISNEbp <- breakpoints(ts(ISNEfrts)~1, breaks = 6)
summary(frISNEbp) 

OSNEfrts <- frOsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
frOSNEbp <- breakpoints(ts(OSNEfrts)~1, breaks = 6)
summary(frOSNEbp) 

IGOMfrts <- frIgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
frIGOMbp <- breakpoints(ts(IGOMfrts)~1, breaks = 6)
summary(frIGOMbp) 

OGOMfrts <- frOgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
frOGOMbp <- breakpoints(ts(OGOMfrts)~1, breaks = 6)
summary(frOGOMbp) 

# combine results 
bind_rows(frCan2, frIsne2, frOsne2, frIgom2, frOgom2) |>
  mutate(Stock = factor(Stock, levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE"))) -> frTSagg


summary(frCANbp)
summary(frISNEbp) 
summary(frOSNEbp) 
summary(frIGOMbp) $RSS[2,]
summary(frOGOMbp) 

fr_rect <- tibble(xm = c(1968, 1968, 1968, 1968,
                         1968, sjOgom2$year[29], 2024),
                  xM = c(2024, 2024, 2024, 2024,
                         sjOgom2$year[29], 2024, 2024),
                  ym = rep(0),
                  yM = rep(0.25),
                  regime = c(1, 1, 1,
                             1, 1, 2, 3),
                  Stock = factor(c("CAN", "ISNE", "OSNE", "IGOM", "OGOM", "OGOM", "OGOM"), 
                                 levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE")))


frTSagg |> #filter(year != 2020 & year != 2023) |> 
  ggplot() +
  facet_grid(.~Stock) +
  geom_rect(data = fr_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2, aes(x = year, y = mean)) + 
  geom_smooth(method = "gam", aes(x = year, y = mean)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Fall Atlantic Rock Crab",
       x = "Year",
       y = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") -> frPlot

sjPlot / fjPlot / srPlot / frPlot + 
  plot_annotation(title = "Annual Mean Interpolated Catch by Stock Area", 
                  tag_levels = "A", 
                  theme = theme(title = element_text(size = 17)))




########### TIME SERIES Cropped ################
########### Spring jonah ##############
# Calculate annual mean Catch per tow by population
cut.sjCan       <- mask(springJonah_rast_crop, Can)
cut.sjIsne      <- mask(springJonah_rast_crop, Isne)
cut.sjOsne      <- mask(springJonah_rast_crop, Osne)
cut.sjIgom       <- mask(springJonah_rast_crop, Igom)
cut.sjOgom       <- mask(springJonah_rast_crop, Ogom)

mean.cut.sjCan       <- terra::global(cut.sjCan, fun = mean, na.rm = T)
sjCan2             <- tibble(mean = mean.cut.sjCan$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("CAN"))

mean.cut.sjIsne       <- terra::global(cut.sjIsne, fun = mean, na.rm = T)
sjIsne2            <- tibble(mean = mean.cut.sjIsne$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("ISNE"))

mean.cut.sjOsne       <- terra::global(cut.sjOsne, fun = mean, na.rm = T)
sjOsne2            <- tibble(mean = mean.cut.sjOsne$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("OSNE"))

mean.cut.sjIgom       <- terra::global(cut.sjIgom, fun = mean, na.rm = T)
sjIgom2            <- tibble(mean = mean.cut.sjIgom$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("IGOM"))

mean.cut.sjOgom       <- terra::global(cut.sjOgom, fun = mean, na.rm = T)
sjOgom2            <- tibble(mean = mean.cut.sjOgom$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("OGOM"))

# Changepoint analysis on each timeseries
CANsjts <- sjCan2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjCANbp <- breakpoints(ts(CANsjts)~1, breaks = 6)

ISNEsjts <- sjIsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjISNEbp <- breakpoints(ts(ISNEsjts)~1, breaks = 6)

OSNEsjts <- sjOsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjOSNEbp <- breakpoints(ts(OSNEsjts)~1, breaks = 6)

IGOMsjts <- sjIgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjIGOMbp <- breakpoints(ts(IGOMsjts)~1, breaks = 6)

OGOMsjts <- sjOgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
sjOGOMbp <- breakpoints(ts(OGOMsjts)~1, breaks = 6)
summary(sjOGOMbp) 

# combine results 
bind_rows(sjCan2, sjIsne2, sjOsne2, sjIgom2, sjOgom2) |>
  mutate(Stock = factor(Stock, levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE"))) -> sjTSagg

summary(sjCANbp)#$RSS[2,]
summary(sjISNEbp) 
summary(sjOSNEbp) 
summary(sjIGOMbp) 
summary(sjOGOMbp) 

sj_rect <- tibble(xm = c(1968, sjCan2$year[31], sjCan2$year[45], 1968, sjIsne2$year[13], 1968, sjOsne2$year[14],
                         1968, sjIgom2$year[29], sjIgom2$year[45], 1968, sjOgom2$year[45]),
                  xM = c(sjCan2$year[31], sjCan2$year[45], 2024, sjIsne2$year[13], 2024, sjOsne2$year[14], 2024,
                         sjIgom2$year[29], sjIgom2$year[45], 2024, sjOgom2$year[45], 2024),
                  ym = rep(0),
                  yM = rep(0.25),
                  regime = c(1, 2, 3, 1, 2, 1, 2, 
                             1, 2, 3, 1, 2),
                  Stock = factor(c("CAN", "CAN", "CAN", "ISNE", "ISNE", "OSNE", "OSNE","IGOM", "IGOM", "IGOM", "OGOM", "OGOM"), 
                                 levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE")))

sjTSagg |> filter(year != 2020 & year != 2023) |> 
  ggplot() +
  facet_grid(.~Stock) +
  geom_rect(data = sj_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2, aes(x = year, y = mean)) + 
  geom_smooth(method = "gam", aes(x = year, y = mean)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Spring Jonah Crab",
       x = "Year",
       y = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") -> sjPlot


########## Fall Jonah #############
cut.fjCan       <- mask(fallJonah_rast_crop, Can)
cut.fjIsne      <- mask(fallJonah_rast_crop, Isne)
cut.fjOsne      <- mask(fallJonah_rast_crop, Osne)
cut.fjIgom       <- mask(fallJonah_rast_crop, Igom)
cut.fjOgom       <- mask(fallJonah_rast_crop, Ogom)

mean.cut.fjCan       <- terra::global(cut.fjCan, fun = mean, na.rm = T)
fjCan2             <- tibble(mean = mean.cut.fjCan$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("CAN"))

mean.cut.fjIsne       <- terra::global(cut.fjIsne, fun = mean, na.rm = T)
fjIsne2            <- tibble(mean = mean.cut.fjIsne$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("ISNE"))

mean.cut.fjOsne       <- terra::global(cut.fjOsne, fun = mean, na.rm = T)
fjOsne2            <- tibble(mean = mean.cut.fjOsne$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("OSNE"))

mean.cut.fjIgom       <- terra::global(cut.fjIgom, fun = mean, na.rm = T)
fjIgom2            <- tibble(mean = mean.cut.fjIgom$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("IGOM"))

mean.cut.fjOgom       <- terra::global(cut.fjOgom, fun = mean, na.rm = T)
fjOgom2            <- tibble(mean = mean.cut.fjOgom$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("OGOM"))

# Changepoint analysis on each timeseries
CANfjts <- fjCan2 |> filter(year > 1968) |> select(mean)
fjCANbp <- breakpoints(ts(CANfjts)~1, breaks = 6)
summary(fjCANbp) 

ISNEfjts <- fjIsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
fjISNEbp <- breakpoints(ts(ISNEfjts)~1, breaks = 6)
summary(fjISNEbp) 

OSNEfjts <- fjOsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
fjOSNEbp <- breakpoints(ts(OSNEfjts)~1, breaks = 6)
summary(fjOSNEbp) 

IGOMfjts <- fjIgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
fjIGOMbp <- breakpoints(ts(IGOMfjts)~1, breaks = 6)
summary(fjIGOMbp) 

OGOMfjts <- fjOgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
fjOGOMbp <- breakpoints(ts(OGOMfjts)~1, breaks = 6)
summary(fjOGOMbp) 

# combine results 
bind_rows(fjCan2, fjIsne2, fjOsne2, fjIgom2, fjOgom2) |>
  mutate(Stock = factor(Stock, levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE"))) -> fjTSagg

summary(fjCANbp)#$RSS[2,]
summary(fjISNEbp) 
summary(fjOSNEbp) 
summary(fjIGOMbp) 
summary(fjOGOMbp) 

fj_rect <- tibble(xm = c(1968, fjCan2$year[30], 1968, fjIsne2$year[42], 1968, fjOsne2$year[40],
                         1968, fjIgom2$year[30], fjIgom2$year[43], 1968, fjOgom2$year[30]),
                  xM = c(fjCan2$year[30], 2024, fjIsne2$year[42], 2024, fjOsne2$year[40], 2024,
                         fjIgom2$year[30], fjIgom2$year[43], 2024, fjOgom2$year[30], 2024),
                  ym = rep(0),
                  yM = rep(0.25),
                  regime = c(1, 2, 1, 2, 1, 2, 
                             1, 2, 3, 1, 2),
                  Stock = factor(c("CAN", "CAN", "ISNE", "ISNE", "OSNE", "OSNE","IGOM", "IGOM", "IGOM", "OGOM", "OGOM"), 
                                 levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE")))
fjTSagg |> #filter(year != 2020 & year != 2023) |> 
  ggplot() +
  facet_grid(.~Stock) +
  geom_rect(data = fj_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2, aes(x = year, y = mean)) + 
  geom_smooth(method = "gam", aes(x = year, y = mean)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Fall Jonah Crab",
       x = "Year",
       y = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") -> fjPlot

########### Spring Rock ##############
# Calculate annual mean Catch per tow by population
cut.srCan       <- mask(springRock_rast_crop, Can)
cut.srIsne      <- mask(springRock_rast_crop, Isne)
cut.srOsne      <- mask(springRock_rast_crop, Osne)
cut.srIgom      <- mask(springRock_rast_crop, Igom)
cut.srOgom      <- mask(springRock_rast_crop, Ogom)

mean.cut.srCan       <- terra::global(cut.srCan, fun = mean, na.rm = T)
srCan2             <- tibble(mean = mean.cut.srCan$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("CAN"))

mean.cut.srIsne       <- terra::global(cut.srIsne, fun = mean, na.rm = T)
srIsne2            <- tibble(mean = mean.cut.srIsne$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("ISNE"))

mean.cut.srOsne       <- terra::global(cut.srOsne, fun = mean, na.rm = T)
srOsne2            <- tibble(mean = mean.cut.srOsne$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("OSNE"))

mean.cut.srIgom       <- terra::global(cut.srIgom, fun = mean, na.rm = T)
srIgom2            <- tibble(mean = mean.cut.srIgom$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("IGOM"))

mean.cut.srOgom       <- terra::global(cut.srOgom, fun = mean, na.rm = T)
srOgom2            <- tibble(mean = mean.cut.srOgom$mean,
                             year = seq(1968, 2024, by = 1),
                             Stock = rep("OGOM"))

# Changepoint analysis on each timeseries
CANsrts <- srCan2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srCANbp <- breakpoints(ts(CANsrts)~1, breaks = 6)
summary(srCANbp) 

ISNEsrts <- srIsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srISNEbp <- breakpoints(ts(ISNEsrts)~1, breaks = 6)
summary(srISNEbp) 

OSNEsrts <- srOsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srOSNEbp <- breakpoints(ts(OSNEsrts)~1, breaks = 6)
summary(srOSNEbp) 

IGOMsrts <- srIgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srIGOMbp <- breakpoints(ts(IGOMsrts)~1, breaks = 6)
summary(srIGOMbp) 

OGOMsrts <- srOgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
srOGOMbp <- breakpoints(ts(OGOMsrts)~1, breaks = 6)
summary(srOGOMbp) 

# combine results 
bind_rows(srCan2, srIsne2, srOsne2, srIgom2, srOgom2) |>
  mutate(Stock = factor(Stock, levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE"))) -> srTSagg


summary(srCANbp)$RSS[2,]
summary(srISNEbp) 
summary(srOSNEbp) 
summary(srIGOMbp) $RSS[2,]
summary(srOGOMbp) 

sr_rect <- tibble(xm = c(1968, 1968, srIsne2$year[8], srIsne2$year[16], 1968, srOsne2$year[8], srOsne2$year[16],
                         1968, 1968, sjOgom2$year[39]),
                  xM = c(2024, sjIsne2$year[8], sjIsne2$year[16], 2024, srOsne2$year[8], srOsne2$year[16], 2024,
                         2024, sjOgom2$year[39], 2024),
                  ym = rep(0),
                  yM = rep(0.25),
                  regime = c(1, 1, 2, 3, 1, 2, 3,  
                             1, 1, 2),
                  Stock = factor(c("CAN", "ISNE", "ISNE", "ISNE", "OSNE", "OSNE", "OSNE", "IGOM", "OGOM", "OGOM"), 
                                 levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE")))


srTSagg |> filter(year != 2020 & year != 2023) |> 
  ggplot() +
  facet_grid(.~Stock) +
  geom_rect(data = sr_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2, aes(x = year, y = mean)) + 
  geom_smooth(method = "gam", aes(x = year, y = mean)) +
  
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Spring Atlantic Rock Crab",
       x = "Year",
       y = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") ->srPlot


########## Fall Rock #############
cut.frCan       <- mask(fallRock_rast_crop, Can)
cut.frIsne      <- mask(fallRock_rast_crop, Isne)
cut.frOsne      <- mask(fallRock_rast_crop, Osne)
cut.frIgom      <- mask(fallRock_rast_crop, Igom)
cut.frOgom      <- mask(fallRock_rast_crop, Ogom)

mean.cut.frCan       <- terra::global(cut.frCan, fun = mean, na.rm = T)
frCan2             <- tibble(mean = mean.cut.frCan$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("CAN"))

mean.cut.frIsne       <- terra::global(cut.frIsne, fun = mean, na.rm = T)
frIsne2            <- tibble(mean = mean.cut.frIsne$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("ISNE"))

mean.cut.frOsne       <- terra::global(cut.frOsne, fun = mean, na.rm = T)
frOsne2            <- tibble(mean = mean.cut.frOsne$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("OSNE"))

mean.cut.frIgom       <- terra::global(cut.frIgom, fun = mean, na.rm = T)
frIgom2            <- tibble(mean = mean.cut.frIgom$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("IGOM"))

mean.cut.frOgom       <- terra::global(cut.frOgom, fun = mean, na.rm = T)
frOgom2            <- tibble(mean = mean.cut.frOgom$mean,
                             year = c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1)),
                             Stock = rep("OGOM"))

# Changepoint analysis on each timeseries
CANfrts <- frCan2 |> filter(year > 1968) |> select(mean)
frCANbp <- breakpoints(ts(CANfrts)~1, breaks = 6)
summary(frCANbp) 

ISNEfrts <- frIsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
frISNEbp <- breakpoints(ts(ISNEfrts)~1, breaks = 6)
summary(frISNEbp) 

OSNEfrts <- frOsne2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
frOSNEbp <- breakpoints(ts(OSNEfrts)~1, breaks = 6)
summary(frOSNEbp) 

IGOMfrts <- frIgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
frIGOMbp <- breakpoints(ts(IGOMfrts)~1, breaks = 6)
summary(frIGOMbp) 

OGOMfrts <- frOgom2 |> filter(year > 1968 & year != 2020 & year != 2023) |> select(mean)
frOGOMbp <- breakpoints(ts(OGOMfrts)~1, breaks = 6)
summary(frOGOMbp) 

# combine results 
bind_rows(frCan2, frIsne2, frOsne2, frIgom2, frOgom2) |>
  mutate(Stock = factor(Stock, levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE"))) -> frTSagg

summary(frCANbp)
summary(frISNEbp) 
summary(frOSNEbp) 
summary(frIGOMbp) $RSS[2,]
summary(frOGOMbp) 

fr_rect <- tibble(xm = c(1968, 1968, 1968, 1968,
                         1968, sjOgom2$year[29], 2024),
                  xM = c(2024, 2024, 2024, 2024,
                         sjOgom2$year[29], 2024, 2024),
                  ym = rep(0),
                  yM = rep(0.25),
                  regime = c(1, 1, 1,
                             1, 1, 2, 3),
                  Stock = factor(c("CAN", "ISNE", "OSNE", "IGOM", "OGOM", "OGOM", "OGOM"), 
                                 levels = c("CAN", "IGOM", "OGOM", "ISNE", "OSNE")))

frTSagg |> #filter(year != 2020 & year != 2023) |> 
  ggplot() +
  facet_grid(.~Stock) +
  geom_rect(data = fr_rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
  scale_fill_distiller(palette = "Greys", direction = -1) + 
  geom_line(linewidth = 2, aes(x = year, y = mean)) + 
  geom_smooth(method = "gam", aes(x = year, y = mean)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  labs(title = "Fall Atlantic Rock Crab",
       x = "Year",
       y = expression(log[10]*"(Catch per tow + 1)")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") -> frPlot

sjPlot / fjPlot / srPlot / frPlot + 
  plot_annotation(title = "Annual Mean Interpolated Catch by Stock Area", 
                  tag_levels = "A", 
                  theme = theme(title = element_text(size = 17))) -> Catch_ts

ggsave(filename = "Catch_ts.png",
       plot = Catch_ts, 
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 16,
       width = 16)
