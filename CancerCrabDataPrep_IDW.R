# Jonah crab and rock crab HSI project

# Everett Rzeszowski
# Friday April 4, 2024

rm(list = ls())

DaisyDuck <- c("gstat", "ggpubr", "gratia", "insurancerating",
               "marmap", "mgcv", "MuMIn", 
               "regclass", "sf", "sp", "strucchange",
               "terra", "tidyterra", "tidyverse", "wesanderson")
invisible(lapply(DaisyDuck, library, character.only = T))

#### DATA IMPORTS #####
# import shape files of interest
NOAA_stat_areas <- st_read("~/Documents/R Projects/Shapefiles/NEFSC_GIS/Statistical_Areas_2010.shp") |>
  dplyr::filter(Id > 450 & Id < 615) 
GOM_stock <- st_read("~/Documents/R Projects/Shapefiles/Lobster_stock_areas copy/stockarealobgom.shp")
GBK_stock <- st_read("~/Documents/R Projects/Shapefiles/Lobster_stock_areas copy/stockarealobnotcanada.shp")
SNE_stock <- st_read("~/Documents/R Projects/Shapefiles/Lobster_stock_areas copy/stockarealobsne.shp")
NE_stocks <- bind_rows(GOM_stock, GBK_stock, SNE_stock, 
                       .id = "source")
NE_stocks$source <- factor(NE_stocks$source, labels = c('GOM','GBK','SNE'))
NE_stocks <- NE_stocks[1:7,]

# exploratory map
b2 = getNOAA.bathy(lon1 = -74.2, lon2 = -65, lat1 = 39, lat2 = 45, 
                   resolution = 1)
bf2 = tibble(fortify.bathy(b2))

countries <- map_data("world")

# import trawl survey data
spring_jonah <- read_csv("~/Documents/R Projects/NEFSC_TRAWL_SURVEY/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVCAT.csv") |>
  filter(SVSPP == 312) 
spring_rock <- read_csv("~/Documents/R Projects/NEFSC_TRAWL_SURVEY/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVCAT.csv") |>
  filter(SVSPP == 313)
spring_meta  <- read_csv("~/Documents/R Projects/NEFSC_TRAWL_SURVEY/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/22561_UNION_FSCS_SVSTA.csv")

fall_jonah <- read_csv("~/Documents/R Projects/NEFSC_TRAWL_SURVEY/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVCAT_RAW.csv") |>
  filter(SVSPP == 312) 
fall_rock <- read_csv("~/Documents/R Projects/NEFSC_TRAWL_SURVEY/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVCAT_RAW.csv") |>
  filter(SVSPP == 313)
fall_meta  <- read_csv("~/Documents/R Projects/NEFSC_TRAWL_SURVEY/22560_NEFSCFallFisheriesIndependentBottomTrawlData/22560_UNION_FSCS_SVSTA_RAW.csv")

ves_conv_Jonah <- read_csv("~/Documents/R Projects/Connectivity_Manuscript_Figures/Data/NEFSC_conversion_factors.csv") |>
  filter(SVSPP == 312)
ves_conv_Rock <- read_csv("~/Documents/R Projects/Connectivity_Manuscript_Figures/Data/NEFSC_conversion_factors.csv") |>
  filter(SVSPP == 313)

################ DATA WRANGLING ################
# data combinations and edits
spring_jonah |>
  full_join(y = spring_meta) |>
  mutate(EXPCATCHNUM = if_else(is.na(EXPCATCHNUM), 0, EXPCATCHNUM),
         EXPCATCHWT = if_else(is.na(EXPCATCHWT), 0, EXPCATCHWT),
         AVGDEPTH = AVGDEPTH * -1) |> # adjust for tows that had no lobster at all
  filter(!is.na(BEGLAT),
         !is.na(EST_YEAR),
         AVGDEPTH > -650,
         AVGDEPTH < -1) |>
  mutate(EXPCATCHWT = if_else(EST_YEAR < 1985,  # Apply conversion factors from Pinksy / OceanAdapt/compile.R
                              EXPCATCHWT * ves_conv_Jonah$DCF_WT,  # * 1 , conversion factor for SRVESSLE == DE is also * 1
                              EXPCATCHWT),
         EXPCATCHWT = if_else(EST_YEAR > 1972 & EST_YEAR < 1982, 
                              EXPCATCHWT / ves_conv_Jonah$GCF_WT,
                              EXPCATCHWT),
         EXPCATCHWT = if_else(SVVESSEL %in% c("HB", "PC"), 
                              EXPCATCHWT / 1, # taken from "big_spring" data.frame from github: Pinskylab/OceanAdapt/compile.R lime 1308
                              EXPCATCHWT)) |>
  group_by(CRUISE6, CRUISE, STRATUM, TOW, STATION, STATUS_CODE, ID) |>
  summarize(TOTAL_CATCH_WT = sum(EXPCATCHWT)) |>
  left_join(y = spring_meta, multiple = "any") -> spring_jonah_agg

spring_rock |>
  full_join(y = spring_meta) |>
  mutate(EXPCATCHNUM = if_else(is.na(EXPCATCHNUM), 0, EXPCATCHNUM),
         EXPCATCHWT = if_else(is.na(EXPCATCHWT), 0, EXPCATCHWT),
         AVGDEPTH = AVGDEPTH * -1) |> # adjust for tows that had no lobster at all
  filter(!is.na(BEGLAT),
         !is.na(EST_YEAR),
         AVGDEPTH > -650,
         AVGDEPTH < -1) |>
  mutate(EXPCATCHWT = if_else(EST_YEAR < 1985,  # Apply conversion factors from Pinksy / OceanAdapt/compile.R
                              EXPCATCHWT * ves_conv_Rock$DCF_WT,  # * 1 , conversion factor for SRVESSLE == DE is also * 1
                              EXPCATCHWT),
         EXPCATCHWT = if_else(EST_YEAR > 1972 & EST_YEAR < 1982, 
                              EXPCATCHWT / ves_conv_Rock$GCF_WT,
                              EXPCATCHWT),
         EXPCATCHWT = if_else(SVVESSEL %in% c("HB", "PC"), 
                              EXPCATCHWT / 3.864, # taken from "big_spring" data.frame from github: Pinskylab/OceanAdapt/compile.R lime 1308
                              EXPCATCHWT)) |>
  group_by(CRUISE6, CRUISE, STRATUM, TOW, STATION, STATUS_CODE, ID) |> 
  summarize(TOTAL_CATCH_WT = sum(EXPCATCHWT)) |>
  left_join(y = spring_meta, multiple = "any")-> spring_rock_agg

fall_jonah |>
  full_join(y = fall_meta) |>
  mutate(EXPCATCHNUM = if_else(is.na(EXPCATCHNUM), 0, EXPCATCHNUM),
         EXPCATCHWT = if_else(is.na(EXPCATCHWT), 0, EXPCATCHWT),
         AVGDEPTH = AVGDEPTH * -1,
         EST_YEAR = as.numeric(EST_YEAR)) |> # adjust for tows that had no lobster at all
  filter(!is.na(BEGLAT),
         !is.na(EST_YEAR),
         AVGDEPTH > -650,
         AVGDEPTH < -1) |>
  mutate(EXPCATCHWT = if_else(EST_YEAR < 1985,  # Apply conversion factors from Pinksy / OceanAdapt/compile.R
                              EXPCATCHWT * ves_conv_Jonah$DCF_WT,  # * 1 , conversion factor for SRVESSLE == DE is also * 1
                              EXPCATCHWT),
         EXPCATCHWT = if_else(EST_YEAR > 1972 & EST_YEAR < 1982, 
                              EXPCATCHWT / ves_conv_Jonah$GCF_WT,
                              EXPCATCHWT),
         EXPCATCHWT = if_else(SVVESSEL %in% c("HB", "PC"), 
                              EXPCATCHWT / 1, # taken from "big_spring" data.frame from github: Pinskylab/OceanAdapt/compile.R lime 1308
                              EXPCATCHWT)) |>
  group_by(CRUISE6, CRUISE, STRATUM, TOW, STATION, STATUS_CODE, ID) |> 
  summarize(TOTAL_CATCH_WT = sum(EXPCATCHWT)) |>
  left_join(y = fall_meta, multiple = "any") |>
  mutate(AVGDEPTH = AVGDEPTH * -1,
         EST_YEAR = as.numeric(EST_YEAR)) |> # adjust for tows that had no lobster at all
  filter(!is.na(BEGLAT),
         !is.na(EST_YEAR),
         AVGDEPTH > -650,
         AVGDEPTH < -1) -> fall_jonah_agg

fall_rock |>
  full_join(y = fall_meta) |>
  mutate(EXPCATCHNUM = if_else(is.na(EXPCATCHNUM), 0, EXPCATCHNUM),
         EXPCATCHWT = if_else(is.na(EXPCATCHWT), 0, EXPCATCHWT),
         AVGDEPTH = AVGDEPTH * -1,
         EST_YEAR = as.numeric(EST_YEAR)) |> # adjust for tows that had no lobster at all
  filter(!is.na(BEGLAT),
         !is.na(EST_YEAR),
         AVGDEPTH > -650,
         AVGDEPTH < -1) |>
  mutate(EXPCATCHWT = if_else(EST_YEAR < 1985,  # Apply conversion factors from Pinksy / OceanAdapt/compile.R
                              EXPCATCHWT * ves_conv_Rock$DCF_WT,  # * 1 , conversion factor for SRVESSLE == DE is also * 1
                              EXPCATCHWT),
         EXPCATCHWT = if_else(EST_YEAR > 1972 & EST_YEAR < 1982, 
                              EXPCATCHWT / ves_conv_Rock$GCF_WT,
                              EXPCATCHWT),
         EXPCATCHWT = if_else(SVVESSEL %in% c("HB", "PC"), 
                              EXPCATCHWT / 2.479, # taken from "big_spring" data.frame from github: Pinskylab/OceanAdapt/compile.R lime 1308
                              EXPCATCHWT)) |>
  group_by(CRUISE6, CRUISE, STRATUM, TOW, STATION, STATUS_CODE, ID) |> 
  summarize(TOTAL_CATCH_WT = sum(EXPCATCHWT)) |>
  left_join(y = fall_meta, multiple = "any")  |>
  mutate(AVGDEPTH = AVGDEPTH * -1,
         EST_YEAR = as.numeric(EST_YEAR)) |> # adjust for tows that had no lobster at all
  filter(!is.na(BEGLAT),
         !is.na(EST_YEAR),
         AVGDEPTH > -650,
         AVGDEPTH < -1) -> fall_rock_agg

############# EXPLORATORY PLOTS ############
spring_jonah_agg |> 
  ggplot(aes(x = log10(TOTAL_CATCH_WT + 1))) +
  geom_histogram() + theme_classic()

fall_jonah_agg |>
  ggplot(aes(x = log10(TOTAL_CATCH_WT + 1))) +
  geom_histogram() + theme_classic()

spring_rock_agg |>
  ggplot(aes(x = log10(TOTAL_CATCH_WT + 1))) +
  geom_histogram() + theme_classic()

fall_rock_agg |>
  ggplot(aes(x = log10(TOTAL_CATCH_WT + 1))) +
  geom_histogram() + theme_classic()

#### without zeros
spring_jonah_agg |> filter(TOTAL_CATCH_WT > 0) |>
  ggplot(aes(x = log10(TOTAL_CATCH_WT + 1))) +
  geom_histogram() + theme_classic()

fall_jonah_agg |> filter(TOTAL_CATCH_WT > 0) |>
  ggplot(aes(x = log10(TOTAL_CATCH_WT + 1))) +
  geom_histogram() + theme_classic()

spring_rock_agg |> filter(TOTAL_CATCH_WT > 0) |>
  ggplot(aes(x = log10(TOTAL_CATCH_WT + 1))) +
  geom_histogram() + theme_classic()

fall_rock_agg |> filter(TOTAL_CATCH_WT > 0) |>
  ggplot(aes(x = log10(TOTAL_CATCH_WT + 1))) +
  geom_histogram() + theme_classic()

#### associating with potential HSI components
###### temp
spring_jonah_agg |> 
  ggplot(aes(x = BOTTEMP, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() + 
  coord_cartesian(xlim = c(0, 30)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Jonah crab: Spring",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Bottom Temp. (°C)") +
  theme_classic() -> jsT

fall_jonah_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = BOTTEMP, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 30)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Jonah crab: Fall",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Bottom Temp. (°C)") +
  theme_classic() -> jfT

spring_rock_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = BOTTEMP, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 30)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Rock crab: Spring",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Bottom Temp. (°C)") +
  theme_classic() -> rsT

fall_rock_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = BOTTEMP, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 30)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Rock crab: Fall",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Bottom Temp. (°C)") +
  theme_classic() -> rfT

(jsT + jfT) / (rsT + rfT) + plot_layout(guides = "collect")

######## salinity
spring_jonah_agg |> 
  ggplot(aes(x = BOTSALIN, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(22, 37)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Jonah crab: Spring",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Bottom Salinity (ppt)") +
  theme_classic() -> jsS

fall_jonah_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = BOTSALIN, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(22, 37)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Jonah crab: Fall",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Bottom Salinity (ppt)") +
  theme_classic() -> jfS

spring_rock_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = BOTSALIN, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(22, 37)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Rock crab: Spring",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Bottom Salinity (ppt)") +
  theme_classic() -> rsS

fall_rock_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = BOTSALIN, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(22, 37)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Rock crab: Fall",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Bottom Salinity (ppt)") +
  theme_classic() -> rfS

(jsS + jfS) / (rsS + rfS) + plot_layout(guides = "collect")

######### depth
spring_jonah_agg |> 
  ggplot(aes(x = AVGDEPTH * -1, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(-620, 0)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Jonah crab: Spring",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Depth (m)") +
  theme_classic() -> jsD

fall_jonah_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = AVGDEPTH, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(-620, 0)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Jonah crab: Fall",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Depth (m)") +
  theme_classic() -> jfD

spring_rock_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = AVGDEPTH * -1 , y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(-620, 0)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Rock crab: Spring",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Depth (m)") +
  theme_classic() -> rsD

fall_rock_agg |> filter(EST_YEAR >= 1968) |>
  ggplot(aes(x = AVGDEPTH, y = log10(TOTAL_CATCH_WT + 1), color = EST_YEAR)) +
  geom_point() +
  coord_cartesian(xlim = c(-620, 0)) +
  scale_color_distiller(palette = "Spectral", na.value = "grey") +
  labs(title = "Rock crab: Fall",
       y = expression("log"[10]~"(Total Catch Weight [kg])"),
       x = "Depth (m)") +
  theme_classic() -> rfD

(jsD + jfD) / (rsD + rfD) + plot_layout(guides = "collect")

############# IDW SECTION ##################
#### JONAH CRAB ####
#### SPRING ########
# spring_jonah_agg |> 
#   ungroup() |>
#   filter(EST_YEAR == 1977) |>
#   select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, total_catch_wt = TOTAL_CATCH_WT) -> points
# 
# # make spatial points data frame
# coords <- points[ , c("lon","lat")]
# data <- points[ , 3]
# crs    <- CRS(crs(NOAA_stat_areas))
# 
# spdf <- SpatialPointsDataFrame(coords = coords,
#                                data = data,
#                                proj4string = crs)
# grd <- as.data.frame(spsample(spdf, "regular", n = 100000))
# names(grd)       <- c("lon", "lat")
# coordinates(grd) <- c("lon", "lat")
# gridded(grd)     <- TRUE  # Create SpatialPixel object
# fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# crs(grd) <- crs(spdf)

bbox_vals <- c(xmin = -77.2, xmax = -65, ymin = 34.9, ymax = 44.8)
res <- 0.1

x_range <- seq(bbox_vals["xmin"], bbox_vals["xmax"], by = res)
y_range <- seq(bbox_vals["ymin"], bbox_vals["ymax"], by = res)

n_cols <- length(x_range)
n_rows <- length(y_range)

gt <- GridTopology(cellcentre.offset = c(bbox_vals["xmin"] + res/2, bbox_vals["ymin"] + res/2),
                   cellsize = c(res, res),
                   cells.dim = c(n_cols - 1, n_rows - 1))  # subtract 1 to account for centers

grd <- SpatialGrid(gt, proj4string = CRS(crs(NOAA_stat_areas)))
# grd_df <- SpatialPixelsDataFrame(grd, data = data.frame(dummy = rep(NA, length(grd))))

SPRINGJONAH_IDW <-list() #log base 10
for(i in 1968:2024){
  print(i)
  spring_jonah_agg |> 
    ungroup() |>
    filter(EST_YEAR == i) |>
    select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, total_catch_wt = TOTAL_CATCH_WT) -> points
  
  # make spatial points data frame
  coords <- points[ , c("lon","lat")]
  data <- points[ , 3]
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data,
                                 proj4string = crs)
  
  trawl.idw1 <- idw(log10(total_catch_wt + 1) ~ 1, spdf, 
                    newdata = grd,
                    idp = 2)
  
  
  trawl.idw2 <- SpatialGridDataFrame(grid = grd,
                                     data = data.frame(pred = trawl.idw1$var1.pred),
                                     proj4string = crs)
  
  trawl.rast       <- rast(trawl.idw2)
  
  ggplot() +
    geom_raster(data = trawl.rast, aes(x = x, y = y, fill = pred)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000)) +
    geom_sf(data = Stocks, linewidth = 0.5, color = 'gray5', fill = "gray8", alpha = 0.01) +
    geom_point(alpha = 0.5, data = spring_jonah_agg |> filter(EST_YEAR == i), 
               aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = log10(TOTAL_CATCH_WT + 1))) +
    scale_color_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "firebrick3") +
    coord_sf(xlim = c(-73.75, -65.5), ylim= c(39.25, 45)) +
    theme_classic() -> p

  ggsave(filename = paste0("SpringJonah_log10_",i,".png", sep = ""),
         plot = p,
         dpi = 750,
         height = 10,
         width = 10,
         path = "~/Documents/R Projects/Cancer_HSI2/IDW Plots/SpringJonah_IDWv2")

  SPRINGJONAH_IDW <- c(SPRINGJONAH_IDW, trawl.rast)
}

fullSpringJonah <- terra::rast(SPRINGJONAH_IDW)
writeRaster(fullSpringJonah, "~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringJonah_Fried.tif",
            overwrite = T)

#### JONAH CRAB ####
#### FALL ########
# fall_jonah_agg |> 
#   ungroup() |>
#   filter(EST_YEAR == 1977) |>
#   select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, total_catch_wt = TOTAL_CATCH_WT) -> points
# 
# # make spatial points data frame
# coords <- points[ , c("lon","lat")]
# data <- points[ , 3]
# crs    <- CRS(crs(NOAA_stat_areas))
# 
# spdf <- SpatialPointsDataFrame(coords = coords,
#                                data = data,
#                                proj4string = crs)
# grd <- as.data.frame(spsample(spdf, "regular", n = 100000))
# names(grd)       <- c("lon", "lat")
# coordinates(grd) <- c("lon", "lat")
# gridded(grd)     <- TRUE  # Create SpatialPixel object
# fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# crs(grd) <- crs(spdf)

FALLJONAH_IDW <- list() #log base 10
years <- c(seq(1968, 2019, by = 1), seq(2021, 2024, by = 1))
for(i in years){
  print(i)
  fall_jonah_agg |> 
    ungroup() |>
    filter(EST_YEAR == i) |>
    select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, total_catch_wt = TOTAL_CATCH_WT) -> points
  
  # make spatial points data frame
  coords <- points[ , c("lon","lat")]
  data <- points[ , 3]
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data,
                                 proj4string = crs)
  
  trawl.idw1 <- idw(log10(total_catch_wt + 1) ~ 1, spdf, 
                    newdata = grd,
                    idp = 2)
  
  
  trawl.idw2 <- SpatialGridDataFrame(grid = grd,
                                     data = data.frame(pred = trawl.idw1$var1.pred),
                                     proj4string = crs)
  
  trawl.rast       <- rast(trawl.idw2)

  ggplot() +
    geom_raster(data = trawl.rast, aes(x = x, y = y, fill = pred)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000)) +
    geom_sf(data = Stocks, linewidth = 0.5, color = 'gray5', fill = "gray8", alpha = 0.01) +
    geom_point(alpha = 0.5, data = fall_jonah_agg |> filter(EST_YEAR == i),
               aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = log10(TOTAL_CATCH_WT + 1))) +
    scale_color_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "firebrick3") +
    coord_sf(xlim = c(-73.75, -65.5), ylim= c(39.25, 45)) +
    theme_classic() -> p

  ggsave(filename = paste0("FallJonah_log10_",i,".png", sep = ""),
         plot = p,
         dpi = 750,
         height = 10,
         width = 10,
         path = "~/Documents/R Projects/Cancer_HSI2/IDW Plots/FallJonah_IDWv2")
  
  FALLJONAH_IDW <- c(FALLJONAH_IDW, trawl.rast)
}

fullFallJonah <- terra::rast(FALLJONAH_IDW)
writeRaster(fullFallJonah, "~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallJonah_Fried.tif",
            overwrite = T)

#### ROCK CRAB ####
#### SPRING ########
# spring_rock_agg |> 
#   ungroup() |>
#   filter(EST_YEAR == 1977) |>
#   select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, total_catch_wt = TOTAL_CATCH_WT) -> points
# 
# # make spatial points data frame
# coords <- points[ , c("lon","lat")]
# data <- points[ , 3]
# crs    <- CRS(crs(NOAA_stat_areas))
# 
# spdf <- SpatialPointsDataFrame(coords = coords,
#                                data = data,
#                                proj4string = crs)
# grd <- as.data.frame(spsample(spdf, "regular", n = 100000))
# names(grd)       <- c("lon", "lat")
# coordinates(grd) <- c("lon", "lat")
# gridded(grd)     <- TRUE  # Create SpatialPixel object
# fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# crs(grd) <- crs(spdf)

SPRINGROCK_IDW <- list() #log base 10
for(i in 1968:2024){
  print(i)
  spring_rock_agg |> 
    ungroup() |>
    filter(EST_YEAR == i) |>
    select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, total_catch_wt = TOTAL_CATCH_WT) -> points
  
  # make spatial points data frame
  coords <- points[ , c("lon","lat")]
  data <- points[ , 3]
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data,
                                 proj4string = crs)
  
  trawl.idw1 <- idw(log10(total_catch_wt + 1) ~ 1, spdf, 
                    newdata = grd,
                    idp = 2)
  
  
  trawl.idw2 <- SpatialGridDataFrame(grid = grd,
                                     data = data.frame(pred = trawl.idw1$var1.pred),
                                     proj4string = crs)
  
  trawl.rast       <- rast(trawl.idw2)

  ggplot() +
    geom_raster(data = trawl.rast, aes(x = x, y = y, fill = pred)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000)) +
    geom_sf(data = Stocks, linewidth = 0.5, color = 'gray5', fill = "gray8", alpha = 0.01) +
    geom_point(alpha = 0.5, data = spring_rock_agg |> filter(EST_YEAR == i), 
               aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = log10(TOTAL_CATCH_WT + 1))) +
    scale_color_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "firebrick3") +
    coord_sf(xlim = c(-73.75, -65.5), ylim= c(39.25, 45)) +
    theme_classic() -> p
  
  ggsave(filename = paste0("SpringRock_log10_",i,".png", sep = ""),
         plot = p,
         dpi = 750,
         height = 10,
         width = 10,
         path = "~/Documents/R Projects/Cancer_HSI2/IDW Plots/SpringRock_IDWv2/")
  
  SPRINGROCK_IDW <- c(SPRINGROCK_IDW, trawl.rast)
}

fullSpringRock <- terra::rast(SPRINGROCK_IDW)
writeRaster(fullSpringRock, "~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringRock_Fried.tif",
            overwrite = T)

#### ROCK CRAB ####
#### FALL ########
# fall_rock_agg |> 
#   ungroup() |>
#   filter(EST_YEAR == 1977) |>
#   select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, total_catch_wt = TOTAL_CATCH_WT) -> points
# 
# # make spatial points data frame
# coords <- points[ , c("lon","lat")]
# data <- points[ , 3]
# crs    <- CRS(crs(NOAA_stat_areas))
# 
# spdf <- SpatialPointsDataFrame(coords = coords,
#                                data = data,
#                                proj4string = crs)
# grd <- as.data.frame(spsample(spdf, "regular", n = 100000))
# names(grd)       <- c("lon", "lat")
# coordinates(grd) <- c("lon", "lat")
# gridded(grd)     <- TRUE  # Create SpatialPixel object
# fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# crs(grd) <- crs(spdf)

FALLROCK_IDW <- list() #log base 10
for(i in years){
  print(i)
  fall_rock_agg |> 
    ungroup() |>
    filter(EST_YEAR == i) |>
    select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, total_catch_wt = TOTAL_CATCH_WT) -> points
  
  # make spatial points data frame
  coords <- points[ , c("lon","lat")]
  data <- points[ , 3]
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data,
                                 proj4string = crs)
  
  trawl.idw1 <- idw(log10(total_catch_wt + 1) ~ 1, spdf, 
                    newdata = grd,
                    idp = 2)
  
  
  trawl.idw2 <- SpatialGridDataFrame(grid = grd,
                                     data = data.frame(pred = trawl.idw1$var1.pred),
                                     proj4string = crs)
  
  trawl.rast       <- rast(trawl.idw2)

  ggplot() +
    geom_raster(data = trawl.rast, aes(x = x, y = y, fill = pred)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000)) +
    geom_sf(data = Stocks, linewidth = 0.5, color = 'gray5', fill = "gray8", alpha = 0.01) +
    geom_point(alpha = 0.5, data = fall_rock_agg |> filter(EST_YEAR == i), 
               aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = log10(TOTAL_CATCH_WT + 1))) +
    scale_color_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "firebrick3") +
    coord_sf(xlim = c(-73.75, -65.5), ylim= c(39.25, 45)) +
    theme_classic() -> p
  
  ggsave(filename = paste0("FallRock_log10_",i,".png", sep = ""),
         plot = p,
         dpi = 750,
         height = 10,
         width = 10,
         path = "~/Documents/R Projects/Cancer_HSI2/IDW Plots/FallRock_IDWv2/")
  
  FALLROCK_IDW <- c(FALLROCK_IDW, trawl.rast)
}

fullFallRock<- terra::rast(FALLROCK_IDW)
writeRaster(fullFallRock, "~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallRock_Fried.tif",
            overwrite = T)

##### ENVIRONMENT #####
# FALL 
# fall_rock_agg |>
#   ungroup() |>
#   filter(EST_YEAR == 1977) |>
#   select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT) -> points
# 
# # make spatial points data frame
# coords <- points[ , c("lon","lat")]
# data <- points[ , 2]
# crs    <- CRS(crs(NOAA_stat_areas))
# 
# spdf <- SpatialPointsDataFrame(coords = coords,
#                                data = data,
#                                proj4string = crs)
# grd <- as.data.frame(spsample(spdf, "regular", n = 100000))
# names(grd)       <- c("lon", "lat")
# coordinates(grd) <- c("lon", "lat")
# gridded(grd)     <- TRUE  # Create SpatialPixel object
# fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# crs(grd) <- crs(spdf)

bbox_vals <- c(xmin = -77.2, xmax = -65, ymin = 34.9, ymax = 44.8)
res <- 0.1

x_range <- seq(bbox_vals["xmin"], bbox_vals["xmax"], by = res)
y_range <- seq(bbox_vals["ymin"], bbox_vals["ymax"], by = res)

n_cols <- length(x_range)
n_rows <- length(y_range)

gt <- GridTopology(cellcentre.offset = c(bbox_vals["xmin"] + res/2, bbox_vals["ymin"] + res/2),
                   cellsize = c(res, res),
                   cells.dim = c(n_cols - 1, n_rows - 1))  # subtract 1 to account for centers

grd <- SpatialGrid(gt, proj4string = CRS(crs(NOAA_stat_areas)))
# grd_df <- SpatialPixelsDataFrame(grd, data = data.frame(dummy = rep(NA, length(grd))))


FALL_TEMP_IDW <- list()
years_env <- c(seq(1992, 2019, by = 1), seq(2021, 2024, by = 1))
for(i in years_env){
  print(i)
  fall_rock_agg |> 
    ungroup() |>
    filter(EST_YEAR == i & !is.na(BOTTEMP)) |>
    select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, BOTTEMP) -> points
  
  # make spatial points data frame
  coords <- points[ , c("lon","lat")]
  data <- points[ , 3]
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data,
                                 proj4string = crs)
  
  temp.idw1 <- idw(BOTTEMP ~ 1, spdf, 
                    newdata = grd,
                    idp = 2)
  
  
  temp.idw2 <- SpatialGridDataFrame(grid = grd,
                                     data = data.frame(pred = temp.idw1$var1.pred),
                                     proj4string = crs)
  
  temp.rast       <- rast(temp.idw2)

  ggplot() +
    geom_raster(data = temp.rast, aes(x = x, y = y, fill = pred)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(0, 20), na.value = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000)) +
    geom_sf(data = Stocks, linewidth = 0.5, color = 'gray5', fill = "gray8", alpha = 0.01) +
    geom_map(data = countries, map = countries, aes(long, lat, map_id = region), fill = "grey") +
    geom_point(alpha = 0.5, data = fall_rock_agg |> filter(EST_YEAR == i & !is.na(BOTTEMP)), 
               aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = BOTTEMP)) +
    scale_color_distiller(palette = "Spectral",  limits = c(0, 20), na.value = "firebrick3") +
    coord_sf(xlim = c(-73.75, -65.5), ylim= c(39.25, 45)) +
    theme_classic() -> p
  
  ggsave(filename = paste0("FallTemp_",i,".png", sep = ""),
         plot = p,
         dpi = 750,
         height = 10,
         width = 10,
         path = "~/Documents/R Projects/Cancer_HSI2/IDW Plots/FallTemp_IDW/")
  
  FALL_TEMP_IDW <- c(FALL_TEMP_IDW, temp.rast)

}

fullFallTemp <- terra::rast(FALL_TEMP_IDW)

writeRaster(fullFallTemp, "~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallTemp.tif",
            overwrite = T)


FALL_SALI_IDW <- list()
for(i in years_env){
  
  fall_rock_agg |> 
    ungroup() |>
    filter(EST_YEAR == i & !is.na(BOTSALIN)) |>
    select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, BOTSALIN) -> points
  
  # make spatial points data frame
  coords <- points[ , c("lon","lat")]
  data <- points[ , 3]
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data,
                                 proj4string = crs)
  
  sali.idw1 <- idw(BOTSALIN~ 1, spdf, 
                   newdata = grd,
                   idp = 2)
  
  
  sali.idw2 <- SpatialGridDataFrame(grid = grd,
                                    data = data.frame(pred = sali.idw1$var1.pred),
                                    proj4string = crs)
  
  sali.rast       <- rast(sali.idw2)
  
  ggplot() +
    geom_raster(data = sali.rast, aes(x = x, y = y, fill = pred)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(28, 36), na.value = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000)) +
    geom_sf(data = Stocks, linewidth = 0.5, color = 'gray5', fill = "gray8", alpha = 0.01) +
    geom_point(alpha = 0.5, data = fall_rock_agg |> filter(EST_YEAR == i), 
               aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = BOTSALIN)) +
    scale_color_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "firebrick3") +
    coord_sf(xlim = c(-73.75, -65.5), ylim= c(39.25, 45)) +
    theme_classic() -> p
  
  ggsave(filename = paste0("FallSali_",i,".png", sep = ""),
         plot = p,
         dpi = 750,
         height = 10,
         width = 10,
         path = "~/Documents/R Projects/Cancer_HSI2/IDW Plots/FallSali_IDW/")
  
  FALL_SALI_IDW <- c(FALL_SALI_IDW, sali.rast)
  
}

fullFallSali <- terra::rast(FALL_SALI_IDW)

writeRaster(fullFallSali, "~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallSali.tif",
            overwrite = T)


# SPRING 
# spring_rock_agg |> 
#   ungroup() |>
#   filter(EST_YEAR == 1977) |>
#   select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT) -> points
# 
# # make spatial points data frame
# coords <- points[ , c("lon","lat")]
# data <- points[ , 3]
# crs    <- CRS(crs(NOAA_stat_areas))
# 
# spdf <- SpatialPointsDataFrame(coords = coords,
#                                data = data,
#                                proj4string = crs)
# grd <- as.data.frame(spsample(spdf, "regular", n = 100000))
# names(grd)       <- c("lon", "lat")
# coordinates(grd) <- c("lon", "lat")
# gridded(grd)     <- TRUE  # Create SpatialPixel object
# fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# crs(grd) <- crs(spdf)

SPRING_TEMP_IDW <- list()
for(i in 1992:2024){
  
  print(i)
  spring_rock_agg |> 
    ungroup() |>
    filter(EST_YEAR == i & !is.na(BOTTEMP)) |>
    select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, BOTTEMP) -> points
  
  # make spatial points data frame
  coords <- points[ , c("lon","lat")]
  data <- points[ , 3]
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data,
                                 proj4string = crs)
  
  temp.idw1 <- idw(BOTTEMP ~ 1, spdf, 
                   newdata = grd,
                   idp = 2)
  
  temp.idw2 <- SpatialGridDataFrame(grid = grd,
                                    data = data.frame(pred = temp.idw1$var1.pred),
                                    proj4string = crs)
  
  temp.rast       <- rast(temp.idw2)
  
  ggplot() +
    geom_raster(data = temp.rast, aes(x = x, y = y, fill = pred)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(0, 20), na.value = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000)) +
    geom_sf(data = Stocks, linewidth = 0.5, color = 'gray5', fill = "gray8", alpha = 0.01) +
    geom_point(alpha = 0.5, data = fall_rock_agg |> filter(EST_YEAR == i), 
               aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = BOTTEMP)) +
    scale_color_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "firebrick3") +
    coord_sf(xlim = c(-73.75, -65.5), ylim= c(39.25, 45)) +
    theme_classic() -> p
  
  ggsave(filename = paste0("SpringTemp_",i,".png", sep = ""),
         plot = p,
         dpi = 750,
         height = 10,
         width = 10,
         path = "~/Documents/R Projects/Cancer_HSI2/IDW Plots/SpringTemp_IDW/")
  
  SPRING_TEMP_IDW <- c(SPRING_TEMP_IDW, temp.rast)

}

fullSpringTemp <- terra::rast(SPRING_TEMP_IDW)
writeRaster(fullSpringTemp, "~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringTemp.tif",
            overwrite = T)

SPRING_SALI_IDW <- list()
for(i in 1992:2024){
  
  print(i)
  spring_rock_agg |> 
    ungroup() |>
    filter(EST_YEAR == i & !is.na(BOTSALIN)) |>
    select(lon = DECDEG_BEGLON, lat = DECDEG_BEGLAT, BOTSALIN) -> points
  
  # make spatial points data frame
  coords <- points[ , c("lon","lat")]
  data <- points[ , 3]
  spdf <- SpatialPointsDataFrame(coords = coords,
                                 data = data,
                                 proj4string = crs)
  
  sali.idw1 <- idw(BOTSALIN~ 1, spdf, 
                   newdata = grd,
                   idp = 2)
  
  sali.idw2 <- SpatialGridDataFrame(grid = grd,
                                    data = data.frame(pred = sali.idw1$var1.pred),
                                    proj4string = crs)
  
  sali.rast       <- rast(sali.idw2)
  
  ggplot() +
    geom_raster(data = sali.rast, aes(x = x, y = y, fill = pred)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(28, 36), na.value = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000)) +
    geom_sf(data = Stocks, linewidth = 0.5, color = 'gray5', fill = "gray8", alpha = 0.01) +
    geom_point(alpha = 0.5, data = fall_rock_agg |> filter(EST_YEAR == i), 
               aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = BOTSALIN)) +
    scale_color_distiller(palette = "Spectral",  limits = c(0, 2), na.value = "firebrick3") +
    coord_sf(xlim = c(-73.75, -65.5), ylim= c(39.25, 45)) +
    theme_classic() -> p
  
  ggsave(filename = paste0("SpringSali_",i,".png", sep = ""),
         plot = p,
         dpi = 750,
         height = 10,
         width = 10,
         path = "~/Documents/R Projects/Cancer_HSI2/IDW Plots/SpringSali_IDW/")
  
  SPRING_SALI_IDW <- c(SPRING_SALI_IDW, sali.rast)
  
}

fullSpringSali <- terra::rast(SPRING_SALI_IDW)
writeRaster(fullSpringSali, "~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringSali.tif")
