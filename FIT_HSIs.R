# Fit HSIs by Stock Area, Season, Sex and with / without Time Period

# Everett Rzeszowski 
# April 21, 2025
################## #################### ############################# ##########
# 1) Sexes combined
# For each stock unit identify generate seasonal mean catch per tow time series (interpolated vs. raw)

# a) Fit HSIs on time series
# b) Identify changepoints in the time series -> fit subsequent HSIs

# 2) Seperate sexes and do the same as above

# Approach: 
# 1) Prepare aggregate trawl survey data and shapefile overlays in /CancerCrabDataPrep_IDW.R
# Save all files to prevent need to constantly re-run
# 2) Develop function to develop time series
# 3) Develop function to fit HSIs
# 4) Export all figures 

############ import data ###########
rm(list = ls())
DaisyDuck <- c("cowplot", "gstat", "ggpubr", "gratia", "insurancerating",
               "marmap", "mgcv", "MuMIn", "patchwork",
               "regclass", "sf", "sp", "strucchange", 
               "terra", "tidyterra", "tidyverse", "wesanderson")
invisible(lapply(DaisyDuck, library, character.only = T))

# IDW CPUE Surfaces
springJonah_rast <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringJonah_Fried.tif")
fallJonah_rast <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallJonah_Fried.tif")
springRock_rast <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringRock_Fried.tif")
fallRock_rast <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallRock_Fried.tif")

# IDW Environment


# Aggregated trawl data 
springJonah_trawl <- read_csv("~/Documents/R Projects/Cancer_HSI2/Data_Aggregate_HSI/SpringJonah.csv")
fallJonah_trawl <- read_csv("~/Documents/R Projects/Cancer_HSI2/Data_Aggregate_HSI/FallJonah.csv")
springRock_trawl <- read_csv("~/Documents/R Projects/Cancer_HSI2/Data_Aggregate_HSI/SpringRock.csv")
fallRock_trawl <- read_csv("~/Documents/R Projects/Cancer_HSI2/Data_Aggregate_HSI/FallRock.csv")

# Shapefiles / spatial layers
NOAA_stat_areas <- st_read("~/Documents/R Projects/Shapefiles/NEFSC_GIS/Statistical_Areas_2010.shp") |>
  dplyr::filter(Id > 460 & Id < 641) 

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

# Add stock units to trawl data
springJonah_trawl |> left_join(Jonah_crab_stocks_rough, by = c("AREA" = "Id")) -> springJonah_trawl
fallJonah_trawl |> left_join(Jonah_crab_stocks_rough, by = c("AREA" = "Id")) -> fallJonah_trawl
springRock_trawl |> left_join(Jonah_crab_stocks_rough, by = c("AREA" = "Id")) -> springRock_trawl
fallRock_trawl |> left_join(Jonah_crab_stocks_rough, by = c("AREA" = "Id")) -> fallRock_trawl

############################ MAP STOCKS ########################################
ggplot() +
  geom_contour(aes(x = x, y = y, z = z), 
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000, -2000), ) + 
  geom_sf(data = NOAA_stat_areas, aes(group = Stock, color = Stock, fill = Stock), alpha = 0.5) + 
  geom_map(data = countries, map = countries, aes(map_id = region), fill = "white", color = "gray") +
  geom_map(data = states, map = states, aes(map_id = region), fill = "white", color = "gray") +
  coord_sf(xlim = c(-76.2, -65), ylim = c(35.5, 45)) +
  labs(title = "Approximate Jonah Crab Stocks", x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 16)

####################################### Write Functions ########################
# initial function doesn't have options for sex, stock unit, season etc. 
# you have to use filtered data

# Two approaches to this as well, we are doing the second for now
# 1) do all VIF tests, varaible class division, and SI calculation in function for varying variables
# 2) do variable class division and SI calc for same variables

fit_tsBreaks <- function(data, years, breaks = 5){
  
  ifelse(class(data) == "SpatRaster",
         
         terra::global(data, fun = mean, na.rm = T) -> data,
         
         data |> filter(EST_YEAR %in% years) |>
           mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
                  logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
           group_by(EST_YEAR) |>
           summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) -> data)
  
  data |> select(mean) -> ts_data
  
  bp_ts <- breakpoints(ts(ts_data) ~ 1, breaks = breaks)
  print(summary(bp_ts))
  
  m1 <- lm(ts(ts_data)~breakfactor(bp_ts, breaks = 1))
  m2 <- lm(ts(ts_data)~breakfactor(bp_ts, breaks = 2))
  m3 <- lm(ts(ts_data)~breakfactor(bp_ts, breaks = 3))
  m4 <- lm(ts(ts_data)~breakfactor(bp_ts, breaks = 4))
  m5 <- lm(ts(ts_data)~breakfactor(bp_ts, breaks = 5))
  
  print(paste0("AIC:", AIC(m1, m2, m3, m4, m5)$AIC))
  print(paste0("BIC:", BIC(m1, m2, m3, m4, m5)$BIC))
  print(paste0("AICc:", MuMIn::AICc(m1, m2, m3, m4, m5)$AICc))
}

fit_HSI <- function(data, Season = "Fall"){
  # make class divisions in data
  data$fisher_TEMP <- fisher(data$BOTTEMP, n = 20)
  data$fisher_SALI <- fisher(data$BOTSALIN, n = 20)
  data$fisher_DEPT <- fisher(data$AVGDEPTH, n = 21, diglab = 10)
  data$fisher_LONG <- fisher(data$DECDEG_BEGLON, n = 21)
  
  # latitude class division varies seasonally 
  ifelse(Season == "Fall",
         data$fisher_LATI <- fisher(data$DECDEG_BEGLAT, n = 20),
         data$fisher_LATI <- fisher(data$DECDEG_BEGLAT, n = 21)
  )
  
  # calculate Suitability Indices
  data |>
    group_by(fisher_TEMP) |>
    summarize(TEMP_CPUE = mean(log10(TOTAL_CATCH_WT + 1))) |>
    mutate(TEMP_SI = ((TEMP_CPUE - min(TEMP_CPUE)) / (max(TEMP_CPUE) - min(TEMP_CPUE)))) -> TEMP_SI
  
  data |>
    group_by(fisher_SALI) |>
    summarize(SALI_CPUE = mean(log10(TOTAL_CATCH_WT + 1))) |>
    filter(!is.na(fisher_SALI)) |>
    mutate(SALI_SI = ((SALI_CPUE - min(SALI_CPUE)) / (max(SALI_CPUE) - min(SALI_CPUE)))) -> SALI_SI
  
  data |>
    group_by(fisher_DEPT) |>
    summarize(DEPT_CPUE = mean(log10(TOTAL_CATCH_WT + 1))) |>
    mutate(DEPT_SI = ((DEPT_CPUE - min(DEPT_CPUE)) / (max(DEPT_CPUE) - min(DEPT_CPUE)))) -> DEPT_SI
  
  data |>
    group_by(fisher_LATI) |>
    summarize(LATI_CPUE = mean(log10(TOTAL_CATCH_WT + 1))) |>
    mutate(LATI_SI = ((LATI_CPUE - min(LATI_CPUE)) / (max(LATI_CPUE) - min(LATI_CPUE)))) -> LATI_SI
  
  data |>
    group_by(fisher_LONG) |>
    summarize(LONG_CPUE = mean(log10(TOTAL_CATCH_WT + 1))) |>
    mutate(LONG_SI = ((LONG_CPUE - min(LONG_CPUE)) / (max(LONG_CPUE) - min(LONG_CPUE)))) -> LONG_SI
  
  # bind and return SIs
  data_SI <- bind_cols(TEMP_SI[1:20,], SALI_SI[1:20,], DEPT_SI[1:20,], LATI_SI[1:20,], LONG_SI[1:20,])
  return(data_SI)
}

get_midpoint <- function(interval) {
  
  nums <- as.numeric(unlist(regmatches(interval, gregexpr("[0-9.]+", interval))))
  mean(nums)
}

get_minpoint <- function(interval) {
  
  nums <- as.numeric(unlist(regmatches(interval, gregexpr("[0-9.]+", interval))))
  min(nums)
}

get_maxpoint <- function(interval) {
  
  nums <- as.numeric(unlist(regmatches(interval, gregexpr("[0-9.]+", interval))))
  max(nums)
}

################### Filter data for HSIs / TSs #################################
# Trawl raw
springJonah_trawl |> filter(Stock == "ISNE") -> springJonah_trawl_ISNE
springJonah_trawl |> filter(Stock == "OSNE") -> springJonah_trawl_OSNE
springJonah_trawl |> filter(Stock == "IGOM") -> springJonah_trawl_IGOM
springJonah_trawl |> filter(Stock == "OGOM") -> springJonah_trawl_OGOM
springJonah_trawl |> filter(Stock == "CAN") -> springJonah_trawl_CAN

fallJonah_trawl |> filter(Stock == "ISNE") -> fallJonah_trawl_ISNE
fallJonah_trawl |> filter(Stock == "OSNE") -> fallJonah_trawl_OSNE
fallJonah_trawl |> filter(Stock == "IGOM") -> fallJonah_trawl_IGOM
fallJonah_trawl |> filter(Stock == "OGOM") -> fallJonah_trawl_OGOM
fallJonah_trawl |> filter(Stock == "CAN") -> fallJonah_trawl_CAN


springRock_trawl |> filter(Stock == "ISNE") -> springRock_trawl_ISNE
springRock_trawl |> filter(Stock == "OSNE") -> springRock_trawl_OSNE
springRock_trawl |> filter(Stock == "IGOM") -> springRock_trawl_IGOM
springRock_trawl |> filter(Stock == "OGOM") -> springRock_trawl_OGOM
springRock_trawl |> filter(Stock == "CAN") -> springRock_trawl_CAN

fallRock_trawl |> filter(Stock == "ISNE") -> fallRock_trawl_ISNE
fallRock_trawl |> filter(Stock == "OSNE") -> fallRock_trawl_OSNE
fallRock_trawl |> filter(Stock == "IGOM") -> fallRock_trawl_IGOM
fallRock_trawl |> filter(Stock == "OGOM") -> fallRock_trawl_OGOM 
fallRock_trawl |> filter(Stock == "CAN") -> fallRock_trawl_CAN

# IDW 
ISNE <- st_as_sf(NOAA_stat_areas |> filter(Stock == "ISNE"))
OSNE <- st_as_sf(NOAA_stat_areas |> filter(Stock == "OSNE"))
IGOM <- st_as_sf(NOAA_stat_areas |> filter(Stock == "IGOM"))
OGOM <- st_as_sf(NOAA_stat_areas |> filter(Stock == "OGOM"))
CAN <- st_as_sf(NOAA_stat_areas |> filter(Stock == "CAN"))

springJonah_rast_ISNE       <- mask(springJonah_rast, ISNE)
springJonah_rast_OSNE       <- mask(springJonah_rast, OSNE)
springJonah_rast_IGOM       <- mask(springJonah_rast, IGOM)
springJonah_rast_OGOM       <- mask(springJonah_rast, OGOM)
springJonah_rast_CAN       <- mask(springJonah_rast, CAN)

fallJonah_rast_ISNE       <- mask(fallJonah_rast, ISNE)
fallJonah_rast_OSNE       <- mask(fallJonah_rast, OSNE)
fallJonah_rast_IGOM       <- mask(fallJonah_rast, IGOM)
fallJonah_rast_OGOM       <- mask(fallJonah_rast, OGOM)
fallJonah_rast_CAN       <- mask(fallJonah_rast, CAN)

springRock_rast_ISNE       <- mask(springRock_rast, ISNE)
springRock_rast_OSNE       <- mask(springRock_rast, OSNE)
springRock_rast_IGOM       <- mask(springRock_rast, IGOM)
springRock_rast_OGOM       <- mask(springRock_rast, OGOM)
springRock_rast_CAN       <- mask(springRock_rast, CAN)

fallRock_rast_ISNE       <- mask(fallRock_rast, ISNE)
fallRock_rast_OSNE       <- mask(fallRock_rast, OSNE)
fallRock_rast_IGOM       <- mask(fallRock_rast, IGOM)
fallRock_rast_OGOM       <- mask(fallRock_rast, OGOM)
fallRock_rast_CAN       <- mask(fallRock_rast, CAN)

############################# Fit Breakpoint Timeseries ########################
fit_tsBreaks(springJonah_trawl_ISNE, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(springJonah_trawl_OSNE, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(springJonah_trawl_IGOM, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(springJonah_trawl_OGOM, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(springJonah_trawl_CAN, years = c(seq(1968, 2019), 2021, 2022, 2024))

################## Visualize timeseries from raw trawl data ####################
################## Spring Jonah Trawl  ################## 
(springJonah_trawl_ISNE |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
   mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
          logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
   group_by(EST_YEAR) |>
   summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
   ggplot() +
   #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
   #scale_fill_distiller(palette = "Greys", direction = -1) + 
   geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
   geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
   coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
   labs(title = "Jonah crab - raw",
        subtitle = "Spring, ISNE",
        x = "Year",
        y = expression(log[10]*"(Catch per tow + 1)")) +
   theme_classic(base_size = 14) +
   theme(legend.position = "none") -> psjISNE_raw)

(springJonah_trawl_OSNE |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - raw",
         subtitle = "Spring, OSNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none") -> psjOSNE_raw)

(springJonah_trawl_IGOM |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - raw",
         subtitle = "Spring, IGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none") -> psjIGOM_raw)

(springJonah_trawl_OGOM |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - raw",
         subtitle = "Spring, OGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 18) +
    theme(legend.position = "none") -> psjOGOM_raw)

################## Fall Jonah Trawl ################## 
fit_tsBreaks(fallJonah_trawl_ISNE, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(fallJonah_trawl_OSNE, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(fallJonah_trawl_IGOM, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(fallJonah_trawl_OGOM, years = c(seq(1968, 2019), 2021, 2022, 2024))

(fallJonah_trawl_ISNE |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - raw",
         subtitle = "Fall, ISNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfjISNE_raw)

(fallJonah_trawl_OSNE |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - raw",
         subtitle = "Fall, OSNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfjOSNE_raw)

(fallJonah_trawl_IGOM |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - raw",
         subtitle = "Fall, IGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfjIGOM_raw)

(fallJonah_trawl_OGOM |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - raw",
         subtitle = "Fall, OGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfjOGOM_raw)

################## Spring Rock Trawl ################## 
fit_tsBreaks(springRock_trawl_ISNE, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(springRock_trawl_OSNE, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(springRock_trawl_IGOM, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(springRock_trawl_OGOM, years = c(seq(1968, 2019), 2021, 2022, 2024))

(springRock_trawl_ISNE |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - raw",
         subtitle = "Spring, ISNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psrISNE_raw)

(springRock_trawl_OSNE |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - raw",
         subtitle = "Spring, OSNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psrOSNE_raw)

(springRock_trawl_IGOM |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - raw",
         subtitle = "Spring, IGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psrIGOM_raw)

(springRock_trawl_OGOM |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - raw",
         subtitle = "Spring, OGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psrOGOM_raw)

################## Fall Rock Trawl ################## 
fit_tsBreaks(fallRock_trawl_ISNE, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(fallRock_trawl_OSNE, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(fallRock_trawl_IGOM, years = c(seq(1968, 2019), 2021, 2022, 2024))
fit_tsBreaks(fallRock_trawl_OGOM, years = c(seq(1968, 2019), 2021, 2022, 2024))

(fallRock_trawl_ISNE |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.15), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - raw",
         subtitle = "Fall, ISNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfrISNE_raw)

(fallRock_trawl_OSNE |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.15), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - raw",
         subtitle = "Fall, OSNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfrOSNE_raw)

(fallRock_trawl_IGOM |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.15), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - raw",
         subtitle = "Fall, IGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfrIGOM_raw)

(fallRock_trawl_OGOM |> filter(EST_YEAR >= 1968 & EST_YEAR != 2020 & EST_YEAR != 2023) |>
    mutate(TOTAL_CATCH_WT = if_else(is.na(TOTAL_CATCH_WT), 0, TOTAL_CATCH_WT),
           logTOTAL_CATCH_WT = log10(TOTAL_CATCH_WT + 1)) |>
    group_by(EST_YEAR) |>
    summarise(mean = mean(logTOTAL_CATCH_WT, na.rm = T)) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = EST_YEAR, y = mean)) + 
    geom_smooth(method = "gam", aes(x = EST_YEAR, y = mean)) +
    coord_cartesian(ylim = c(0, 0.15), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - raw",
         subtitle = "Fall, OGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfrOGOM_raw)

########################## Spring Jonah Raster #################################
mean.sj.ISNE          <- terra::global(springJonah_rast_ISNE, fun = mean, na.rm = T)
sj.ISNE               <- tibble(mean = mean.sj.ISNE$mean,
                                year = seq(1968, 2024, by = 1))
mean.sj.OSNE          <- terra::global(springJonah_rast_OSNE, fun = mean, na.rm = T)
sj.OSNE               <- tibble(mean = mean.sj.OSNE$mean,
                                year = seq(1968, 2024, by = 1))
mean.sj.IGOM          <- terra::global(springJonah_rast_IGOM, fun = mean, na.rm = T)
sj.IGOM               <- tibble(mean = mean.sj.IGOM$mean,
                                year = seq(1968, 2024, by = 1))
mean.sj.OGOM          <- terra::global(springJonah_rast_OGOM, fun = mean, na.rm = T)
sj.OGOM               <- tibble(mean = mean.sj.OGOM$mean,
                                year = seq(1968, 2024, by = 1))

(sj.ISNE |> filter(year >= 1968 & year != 2020 & year != 2023) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - interpolated",
         subtitle = "Spring, ISNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psjISNE_rast)

(sj.OSNE|> filter(year >= 1968 & year != 2020 & year != 2023) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - interpolated",
         subtitle = "Spring, OSNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psjOSNE_rast)

(sj.IGOM |> filter(year >= 1968 & year != 2020 & year != 2023)|>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - interpolated",
         subtitle = "Spring, IGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psjIGOM_rast)

(sj.OGOM |> filter(year >= 1968 & year != 2020 & year != 2023) |> 
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - interpolated",
         subtitle = "Spring, OGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psjOGOM_rast)

#############
mean.fj.ISNE          <- terra::global(fallJonah_rast_ISNE, fun = mean, na.rm = T)
fj.ISNE               <- tibble(mean = mean.fj.ISNE$mean,
                                year = c(seq(1968, 2023)))
mean.fj.OSNE          <- terra::global(fallJonah_rast_OSNE, fun = mean, na.rm = T)
fj.OSNE               <- tibble(mean = mean.fj.OSNE$mean,
                                year = seq(1968, 2023, by = 1))
mean.fj.IGOM          <- terra::global(fallJonah_rast_IGOM, fun = mean, na.rm = T)
fj.IGOM               <- tibble(mean = mean.fj.IGOM$mean,
                                year = seq(1968, 2023, by = 1))
mean.fj.OGOM          <- terra::global(fallJonah_rast_OGOM, fun = mean, na.rm = T)
fj.OGOM               <- tibble(mean = mean.fj.OGOM$mean,
                                year = seq(1968, 2023, by = 1))

(fj.ISNE |> filter(year >= 1968 & year != 2020 & year != 2023) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - interpolated",
         subtitle = "Fall, ISNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfjISNE_rast)

(fj.OSNE|> filter(year >= 1968 & year != 2020 & year != 2023) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - interpolated",
         subtitle = "Fall, OSNE",
         x = "Year",
         y = bquote("log10(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfjOSNE_rast)

(fj.IGOM |> filter(year >= 1968 & year != 2020 & year != 2023)|>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - interpolated",
         subtitle = "Fall, IGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfjIGOM_rast)

(fj.OGOM |> filter(year >= 1968 & year != 2020 & year != 2023) |> 
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.25), xlim = c(1968, 2024)) +
    labs(title = "Jonah crab - interpolated",
         subtitle = "Fall, OGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfjOGOM_rast)

#############
mean.sr.ISNE          <- terra::global(springRock_rast_ISNE, fun = mean, na.rm = T)
sr.ISNE               <- tibble(mean = mean.sr.ISNE$mean,
                                year = seq(1968, 2024, by = 1))
mean.sr.OSNE          <- terra::global(springRock_rast_OSNE, fun = mean, na.rm = T)
sr.OSNE               <- tibble(mean = mean.sr.OSNE$mean,
                                year = seq(1968, 2024, by = 1))
mean.sr.IGOM          <- terra::global(springRock_rast_IGOM, fun = mean, na.rm = T)
sr.IGOM               <- tibble(mean = mean.sr.IGOM$mean,
                                year = seq(1968, 2024, by = 1))
mean.sr.OGOM          <- terra::global(springRock_rast_OGOM, fun = mean, na.rm = T)
sr.OGOM               <- tibble(mean = mean.sr.OGOM$mean,
                                year = seq(1968, 2024, by = 1))

(sr.ISNE |> filter(year >= 1968 & year != 2020 & year != 2023) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.3), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - interpolated",
         subtitle = "Spring, ISNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psrISNE_rast)

(sr.OSNE|> filter(year >= 1968 & year != 2020 & year != 2023) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.3), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - interpolated",
         subtitle = "Spring, OSNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psrOSNE_rast)

(sr.IGOM |> filter(year >= 1968 & year != 2020 & year != 2023)|>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.3), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - interpolated",
         subtitle = "Spring, IGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psrIGOM_rast)

(sr.OGOM |> filter(year >= 1968 & year != 2020 & year != 2023) |> 
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.3), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - interpolated",
         subtitle = "Spring, OGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> psrOGOM_rast)

#############
mean.fr.ISNE          <- terra::global(fallRock_rast_ISNE, fun = mean, na.rm = T)
fr.ISNE               <- tibble(mean = mean.fr.ISNE$mean,
                                year = seq(1968, 2023, by = 1))
mean.fr.OSNE          <- terra::global(fallRock_rast_OSNE, fun = mean, na.rm = T)
fr.OSNE               <- tibble(mean = mean.fr.OSNE$mean,
                                year = seq(1968, 2023, by = 1))
mean.fr.IGOM          <- terra::global(fallRock_rast_IGOM, fun = mean, na.rm = T)
fr.IGOM               <- tibble(mean = mean.fr.IGOM$mean,
                                year = seq(1968, 2023, by = 1))
mean.fr.OGOM          <- terra::global(fallRock_rast_OGOM, fun = mean, na.rm = T)
fr.OGOM               <- tibble(mean = mean.fr.OGOM$mean,
                                year = seq(1968, 2023, by = 1))

(fr.ISNE |> filter(year >= 1968 & year != 2020 & year != 2023) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.2), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - interpolated",
         subtitle = "Fall, ISNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfrISNE_rast)

(fr.OSNE|> filter(year >= 1968 & year != 2020 & year != 2023) |>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.2), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - interpolated",
         subtitle = "Fall, OSNE",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfrOSNE_rast)

(fr.IGOM |> filter(year >= 1968 & year != 2020 & year != 2023)|>
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.2), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - interpolated",
         subtitle = "Fall, IGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfrIGOM_rast)

(fr.OGOM |> filter(year >= 1968 & year != 2020 & year != 2023) |> 
    ggplot() +
    #geom_rect(data = rect, aes(xmin = xm, ymin = ym, xmax = xM, ymax = yM, fill = regime), alpha = 0.2) +
    #scale_fill_distiller(palette = "Greys", direction = -1) + 
    geom_line(linewidth = 2, aes(x = year, y = mean)) + 
    geom_smooth(method = "gam", aes(x = year, y = mean)) +
    coord_cartesian(ylim = c(0, 0.2), xlim = c(1968, 2024)) +
    labs(title = "Rock crab - interpolated",
         subtitle = "Fall, OGOM",
         x = "Year",
         y = expression(log[10]*"(Catch per tow + 1)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none") -> pfrOGOM_rast)

######################## SOME SUMMARY FIGURES ##################################
(psjISNE_raw + psjOSNE_raw) / (psjIGOM_raw + psjOGOM_raw) + plot_annotation(tag_levels = "A")
(pfjISNE_raw + pfjOSNE_raw) / (pfjIGOM_raw + pfjOGOM_raw) + plot_annotation(tag_levels = "A")
(psrISNE_raw + psrOSNE_raw) / (psrIGOM_raw + psrOGOM_raw) + plot_annotation(tag_levels = "A")
(pfrISNE_raw + pfrOSNE_raw) / (pfrIGOM_raw + pfrOGOM_raw) + plot_annotation(tag_levels = "A")

(psjISNE_rast + psjOSNE_rast) / (psjIGOM_rast + psjOGOM_rast) + plot_annotation(tag_levels = "A")
(pfjISNE_rast + pfjOSNE_rast) / (pfjIGOM_rast + pfjOGOM_rast) + plot_annotation(tag_levels = "A")
(psrISNE_rast + psrOSNE_rast) / (psrIGOM_rast + psrOGOM_rast) + plot_annotation(tag_levels = "A")
(pfrISNE_rast + pfrOSNE_rast) / (pfrIGOM_rast + pfrOGOM_rast) + plot_annotation(tag_levels = "A")

(psjISNE_raw + psjISNE_rast) / (psjOSNE_raw + psjOSNE_rast)
(psjIGOM_raw + psjIGOM_rast) / (psjOGOM_raw + psjOGOM_rast)

(pfjISNE_raw + pfjISNE_rast) / (pfjOSNE_raw + pfjOSNE_rast)
(pfjIGOM_raw + pfjIGOM_rast) / (pfjOGOM_raw + pfjOGOM_rast)

(psrISNE_raw + psrISNE_rast) / (psrOSNE_raw + psrOSNE_rast)
(psrIGOM_raw + psrIGOM_rast) / (psrOGOM_raw + psrOGOM_rast)

(pfrISNE_raw + pfrISNE_rast) / (pfrOSNE_raw + pfrOSNE_rast)
(pfrIGOM_raw + pfrIGOM_rast) / (pfrOGOM_raw + pfrOGOM_rast)

######################################### Fitting HSIs #########################
# step 1) determine catch value / CPUE calculation -- vessel conversion factors used, CPUE is overkill without adding other surveys
# step 2) VIF test for variables to include
# step 3 and 4) divide variables and calculate SIs (this is done in fc fit_HSI())
############## ########## SPRING JONAH############## ############## ############## 
# SPRING JONAH ISNE
SJ_ISNE1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springJonah_trawl_ISNE)

summary(SJ_ISNE1)
performance::check_collinearity(SJ_ISNE1)

# SPRING JONAH OSNE
SJ_OSNE1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springJonah_trawl_OSNE)

summary(SJ_OSNE1)
performance::check_collinearity(SJ_OSNE1)

# SPRING JONAH IGOM 
SJ_IGOM1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springJonah_trawl_IGOM)

summary(SJ_IGOM1)
performance::check_collinearity(SJ_IGOM1)

# SPRING JONAH OGOM
SJ_OGOM1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springJonah_trawl_OGOM)

summary(SJ_OGOM1)
performance::check_collinearity(SJ_OGOM1)

# SPRING JONAH CAN
SJ_CAN1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springJonah_trawl_CAN)

summary(SJ_CAN1)
performance::check_collinearity(SJ_CAN1)

SJ_CAN3 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + AVGDEPTH, data = springJonah_trawl_CAN)
performance::check_collinearity(SJ_CAN3)

############################ FIT SPRING JONAH HSIs #############################
SJ_HSI_ISNE <- fit_HSI(springJonah_trawl_ISNE, Season = "Spring")
SJ_HSI_OSNE <- fit_HSI(springJonah_trawl_OSNE, Season = "Spring")
SJ_HSI_IGOM <- fit_HSI(springJonah_trawl_IGOM, Season = "Spring")
SJ_HSI_OGOM <- fit_HSI(springJonah_trawl_OGOM, Season = "Spring")
SJ_HSI_CAN <- fit_HSI(springJonah_trawl_CAN, Season = "Spring")

SpringJonahHSI <- bind_rows(SJ_HSI_ISNE, SJ_HSI_OSNE, SJ_HSI_IGOM, SJ_HSI_OGOM, SJ_HSI_CAN, .id = "Stock") |>
  mutate(Stock = factor(Stock, levels = c(5, 4, 3, 2, 1), labels = c("CAN", "IGOM", "OGOM", "OSNE", "ISNE"))) |>
  group_by(row_number()) |>
  mutate(midTemp = get_midpoint(fisher_TEMP),
         midSali = get_midpoint(fisher_SALI), 
         midDept = get_midpoint(fisher_DEPT),
         midLong = get_midpoint(fisher_LONG),
         midLati = get_midpoint(fisher_LATI)) |> 
  ungroup() 

(SpringJonahHSI |> filter(Stock != "ISNE") |>
  ggplot(aes(x = midTemp, y = TEMP_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(2, 25)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Temperature",
       x = "Temperature (°C)", y = "Suitability Index") +
  theme_classic(base_size = 18) +
    theme(legend.position = "none") -> sj_TSI)

(SpringJonahHSI |> filter(Stock != "ISNE") |>
  ggplot(aes(x = midSali, y = SALI_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(24.5, 36)) +
    geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
    labs(title = "Salinity",
       x = "Salinity (psu)", y = "") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") -> sj_SSI)

SpringJonahHSI |> filter(Stock != "ISNE") |>
  ggplot(aes(x = midDept, y = DEPT_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(0, 500)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Depth",
       x = "Depth (m)", y = "") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") -> sj_DSI

((sj_TSI + sj_SSI + sj_DSI) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
                  title = "Spring Jonah Crab",
                  theme = theme(title = element_text(size = 18))) -> SJ_HSI)
############## ############## FALL JONAH ############## ############## ############## 
# FALL JONAH ISNE
FJ_ISNE1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallJonah_trawl_ISNE)

summary(FJ_ISNE1)
performance::check_collinearity(FJ_ISNE1)

# FALL JONAH OSNE
FJ_OSNE1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallJonah_trawl_OSNE)

summary(FJ_OSNE1)
performance::check_collinearity(FJ_OSNE1)

# FALL JONAH IGOM 
FJ_IGOM1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallJonah_trawl_IGOM)

summary(FJ_IGOM1)
performance::check_collinearity(FJ_IGOM1)
  
FJ_IGOM3 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + AVGDEPTH, data = fallJonah_trawl_IGOM)
performance::check_collinearity(FJ_IGOM3)

# FALL JONAH OGOM
FJ_OGOM1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallJonah_trawl_OGOM)

summary(FJ_OGOM1)
performance::check_collinearity(FJ_OGOM1)

FJ_OGOM3 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + AVGDEPTH, data = fallJonah_trawl_OGOM)
performance::check_collinearity(FJ_OGOM3)

# FALL JONAH CAN
FJ_CAN1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallJonah_trawl_CAN)

summary(FJ_CAN1)
performance::check_collinearity(FJ_CAN1)

FJ_CAN3 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + AVGDEPTH, data = fallJonah_trawl_CAN)
performance::check_collinearity(FJ_CAN3)
############################ FIT FALL JONAH HSIs #############################
FJ_HSI_ISNE <- fit_HSI(fallJonah_trawl_ISNE, Season = "Fall")
FJ_HSI_OSNE <- fit_HSI(fallJonah_trawl_OSNE, Season = "Fall")
FJ_HSI_IGOM <- fit_HSI(fallJonah_trawl_IGOM, Season = "Fall")
FJ_HSI_OGOM <- fit_HSI(fallJonah_trawl_OGOM, Season = "Fall")
FJ_HSI_CAN <- fit_HSI(fallJonah_trawl_CAN, Season = "Fall")


FallJonahHSI <- bind_rows(FJ_HSI_ISNE, FJ_HSI_OSNE, FJ_HSI_IGOM, FJ_HSI_OGOM, FJ_HSI_CAN, .id = "Stock") |>
  mutate(Stock = factor(Stock, levels = c(5, 4, 3, 2, 1), labels = c("CAN", "IGOM", "OGOM", "OSNE", "ISNE"))) |>
  group_by(row_number()) |>
  mutate(midTemp = get_midpoint(fisher_TEMP),
         midSali = get_midpoint(fisher_SALI), 
         midDept = get_midpoint(fisher_DEPT),
         midLong = get_midpoint(fisher_LONG),
         midLati = get_midpoint(fisher_LATI)) |> 
  ungroup() 

(FallJonahHSI |>  filter(Stock != "ISNE") |>
  ggplot(aes(x = midTemp, y = TEMP_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(2, 25)) +
    geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
    labs(title = "Temperature",
       x = "Temperature (°C)", y = "Suitability Index") +
  theme_classic(base_size = 18)  +
  theme(legend.position = "none") -> fj_TSI)

FallJonahHSI |>  filter(Stock != "ISNE") |>
  ggplot(aes(x = midSali, y = SALI_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(24.5, 36)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Salinity",
       x = "Salinity (psu)", y = "") +
  theme_classic(base_size = 18)  +
  theme(legend.position = "none") -> fj_SSI

FallJonahHSI |>  filter(Stock != "ISNE") |>
  ggplot(aes(x = midDept, y = DEPT_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(0, 500)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Depth",
       x = "Depth (m)", y = "") +
  theme_classic(base_size = 18)  +
  theme(legend.position = "none") -> fj_DSI

((fj_TSI + fj_SSI + fj_DSI) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
                  title = "Fall Jonah Crab",
                  theme = theme(title = element_text(size = 20))) -> FJ_HSI)

############## ############## SPRING ROCK############## ############## ############## 
# SPRING ROCK ISNE
SR_ISNE1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springRock_trawl_ISNE)

summary(SR_ISNE1)
performance::check_collinearity(SR_ISNE1)


# SPRING ROCK OSNE
SR_OSNE1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springRock_trawl_OSNE)

summary(SR_OSNE1)
performance::check_collinearity(SR_OSNE1)

# SPRING ROCK IGOM 
SR_IGOM1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springRock_trawl_IGOM)

summary(SR_IGOM1)
performance::check_collinearity(SR_IGOM1)

# SPRING ROCK OGOM
SR_OGOM1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springRock_trawl_OGOM)

summary(SR_OGOM1)
performance::check_collinearity(SR_OGOM1)

# SPRING ROCK CAN
SR_CAN1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = springRock_trawl_CAN)

summary(SR_OGOM1)
performance::check_collinearity(SR_OGOM1)

############################ FIT Spring ROCK HSIs #############################
SR_HSI_ISNE <- fit_HSI(springRock_trawl_ISNE, Season = "Spring")
SR_HSI_OSNE <- fit_HSI(springRock_trawl_OSNE, Season = "Spring")
SR_HSI_IGOM <- fit_HSI(springRock_trawl_IGOM, Season = "Spring")
SR_HSI_OGOM <- fit_HSI(springRock_trawl_OGOM, Season = "Spring")
SR_HSI_CAN <- fit_HSI(springRock_trawl_CAN, Season = "Spring")


SpringRockHSI <- bind_rows(SR_HSI_ISNE, SR_HSI_OSNE, SR_HSI_IGOM, SR_HSI_OGOM, SR_HSI_CAN, .id = "Stock") |>
  mutate(Stock = factor(Stock, levels = c(5, 4, 3, 2, 1), labels = c("CAN", "IGOM", "OGOM", "OSNE", "ISNE"))) |>
  group_by(row_number()) |>
  mutate(midTemp = get_midpoint(fisher_TEMP),
         midSali = get_midpoint(fisher_SALI), 
         midDept = get_midpoint(fisher_DEPT),
         midLong = get_midpoint(fisher_LONG),
         midLati = get_midpoint(fisher_LATI)) |> 
  ungroup() 

SpringRockHSI |>   filter(Stock != "ISNE") |>
  ggplot(aes(x = midTemp, y = TEMP_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(2, 25)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Temperature",
       x = "Temperature (°C)", y = "Suitability Index") +
  theme_classic(base_size = 18)  +
  theme(legend.position = "none") -> sr_TSI

SpringRockHSI |>  filter(Stock != "ISNE") |>
  ggplot(aes(x = midSali, y = SALI_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(24.5, 36)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Salinity",
       x = "Salinity (psu)", y = "") +
  theme_classic(base_size = 18)  +
  theme(legend.position = "none") -> sr_SSI

SpringRockHSI |>  filter(Stock != "ISNE") |>
  ggplot(aes(x = midDept, y = DEPT_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(0, 500)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Depth",
       x = "Depth (m)", y = "") +
  theme_classic(base_size = 18)  +
  theme(legend.position = "none") -> sr_DSI


((sr_TSI + sr_SSI + sr_DSI) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A",
                    title = "Spring Atl. Rock Crab",
                    theme = theme(title = element_text(size = 20))) -> SR_HSI)

############## FALL ROCK ############
# FALL ROCK ISNE
FR_ISNE1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallRock_trawl_ISNE)

summary(FR_ISNE1)
performance::check_collinearity(FR_ISNE1)

# FALL ROCK OSNE
FR_OSNE1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallRock_trawl_OSNE)

summary(FR_OSNE1)
performance::check_collinearity(FR_OSNE1)

# FALL ROCK IGOM 
FR_IGOM1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallRock_trawl_IGOM)

summary(FR_IGOM1)
performance::check_collinearity(FR_IGOM1)

# FALL ROCK OGOM
FR_OGOM1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallRock_trawl_OGOM)


summary(FR_OGOM1)
performance::check_collinearity(FR_OGOM1)

# FALL ROCK CAN
FR_CAN1 <- lm(TOTAL_CATCH_WT ~ BOTTEMP + BOTSALIN + AVGDEPTH, data = fallRock_trawl_CAN)


summary(FR_CAN1)
performance::check_collinearity(FR_CAN1)


############################ FIT FALL ROCK HSIs #############################
FR_HSI_ISNE <- fit_HSI(fallRock_trawl_ISNE, Season = "Fall")
FR_HSI_OSNE <- fit_HSI(fallRock_trawl_OSNE, Season = "Fall")
FR_HSI_IGOM <- fit_HSI(fallRock_trawl_IGOM, Season = "Fall")
FR_HSI_OGOM <- fit_HSI(fallRock_trawl_OGOM, Season = "Fall") 
FR_HSI_CAN <- fit_HSI(fallRock_trawl_CAN, Season = "Fall") 

FallRockHSI <- bind_rows(FR_HSI_ISNE, FR_HSI_OSNE, FR_HSI_IGOM, FR_HSI_OGOM, FR_HSI_CAN, .id = "Stock") |>
  mutate(Stock = factor(Stock, levels = c(5, 4, 3, 2, 1), labels = c("CAN", "IGOM", "OGOM", "OSNE", "ISNE"))) |>
  group_by(row_number()) |>
  mutate(midTemp = get_midpoint(fisher_TEMP),
         midSali = get_midpoint(fisher_SALI), 
         midDept = get_midpoint(fisher_DEPT),
         midLong = get_midpoint(fisher_LONG),
         midLati = get_midpoint(fisher_LATI)) |> 
  ungroup() 

FallRockHSI |> filter(Stock != "ISNE") |>
  ggplot(aes(x = midTemp, y = TEMP_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(2, 25)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Temperature",
       x = "Temperature (°C)", y = "Suitability Index") +
  theme_classic(base_size = 18)  +
  theme(legend.position = "none") -> fr_TSI

FallRockHSI |> filter(Stock != "ISNE") |>
  ggplot(aes(x = midSali, y = SALI_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(24.5, 36)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Salinity",
       x = "Salinity (psu)", y = "") +
  theme_classic(base_size = 18)  +
  theme(legend.position = "none") -> fr_SSI

FallRockHSI |> filter(Stock != "ISNE") |>
  ggplot(aes(x = midDept, y = DEPT_SI, group = Stock, color = Stock)) +
  facet_grid(rows = vars(Stock)) + 
  scale_color_manual(values = c("red2", "steelblue2", "yellow", "pink")) +
  geom_smooth(method = "gam", se = F, linewidth = 2) +
  coord_cartesian(ylim = c(-0.1, 1), xlim = c(0, 500)) +
  geom_hline(yintercept = c(0, 0.2, 0.5, 0.8, 1), lty = 2, color = "gray95") +
  labs(title = "Depth",
       x = "Depth (m)", y = "") +
  theme_classic(base_size = 18)  -> fr_DSI

###### HSI Figure 5 ########
((sj_TSI + sj_SSI + sj_DSI) +
    plot_annotation(tag_levels = "A",
                    title = "Spring Jonah Crab",
                    theme = theme(title = element_text(size = 18))) -> SJ_HSI)

((fj_TSI + fj_SSI + fj_DSI) +
    plot_annotation(tag_levels = list(c("D", "E", "F")),
                    title = "Fall Jonah Crab",
                    theme = theme(title = element_text(size = 20))) -> FJ_HSI)


((sr_TSI + sr_SSI + sr_DSI) + 
    plot_annotation(tag_levels = list(c("G", "H", "I")),
                    title = "Spring Atl. Rock Crab",
                    theme = theme(title = element_text(size = 20))) -> SR_HSI)


((fr_TSI + fr_SSI + fr_DSI) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = list(c("J", "K", "L")),
                    title = "Fall Atl. Rock Crab",
                    theme = theme(title = element_text(size = 18))) -> FR_HSI)

legend2 <- get_legend(FR_HSI)
plot_grid(SJ_HSI, FJ_HSI, SR_HSI, FR_HSI + theme(legend.position = "none"), nrow = 4) -> HSIs

plot_grid(HSIs, legend2,
          rel_widths = c(1, 0.2)) +
  plot_annotation(title = "Seasonal Cancer Crab Suitability Indices", 
                  theme = theme(title = element_text(size = 20))) -> HSI_MSfigure

ggsave("fig5fix_v4_hodgdon.png",
       HSI_MSfigure,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 24,
       width = 14)

############################# HINDCAST HSIs ####################################
# HSI must be added to stock specific rasters, as stock specific data frame, then recombine the HSI rasters across season and species
# Spring evn. raster masking
springTemp_rast_ISNE       <- mask(spngTemp_rast_f, ISNE)
springTemp_rast_OSNE       <- mask(spngTemp_rast_f, OSNE)
springTemp_rast_IGOM       <- mask(spngTemp_rast_f, IGOM)
springTemp_rast_OGOM       <- mask(spngTemp_rast_f, OGOM)
springTemp_rast_CAN       <- mask(spngTemp_rast_f, CAN)

springSalt_rast_ISNE       <- mask(spngSali_rast_f, ISNE)
springSalt_rast_OSNE       <- mask(spngSali_rast_f, OSNE)
springSalt_rast_IGOM       <- mask(spngSali_rast_f, IGOM)
springSalt_rast_OGOM       <- mask(spngSali_rast_f, OGOM)
springSalt_rast_CAN       <- mask(spngSali_rast_f, CAN)

d2 = getNOAA.bathy(lon1 = -77.26, lon2 = -63.32, lat1 = 34.98, lat2 = 44.78, 
                   resolution = 6)
df2 = tibble(fortify.bathy(d2))
depth_raster <- rast(df2, type = "xyz")
summary(values(depth_raster))

crs(depth_raster) <- "EPSG:4269"

depth_raster <- mask(depth_raster, Strata_2)
(depth_raster <- depth_raster*-1)
summary(values(depth_raster))
plot(depth_raster[[1]])
depth_rast_ISNE <- mask(depth_raster, ISNE)
depth_rast_OSNE <- mask(depth_raster, OSNE)
depth_rast_IGOM <- mask(depth_raster, IGOM)
depth_rast_OGOM <- mask(depth_raster, OGOM)
depth_rast_CAN <- mask(depth_raster, CAN)

plot(depth_raster[[1]])
summary(values(depth_rast_CAN))


# Create a function to assign suitability based on ranges
depth_lookup_fun <- function(x, data) {
  sapply(x, function(val) {
    matched <- data[val >= data[,1] & val <= data[,2], 3]
    if (length(matched) > 0) matched[1] else NA  # NA for out-of-range depths
  })
}

## SPRING JONAH
# Can
SJ_HSI_CAN |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SJ_HSI_CAN

sjCAN_mat1 <- as.matrix(SJ_HSI_CAN[, c("mintemp", "maxtemp", "TEMP_SI")])
sjCAN_TSI_rast <- terra::classify(springTemp_rast_CAN, rcl = sjCAN_mat1, include.lowest = T)

sjCAN_mat2 <- as.matrix(SJ_HSI_CAN[, c("minsalt", "maxsalt", "SALI_SI")])
sjCAN_SSI_rast <- terra::classify(springSalt_rast_CAN, rcl = sjCAN_mat2, include.lowest = T)

sjCAN_mat3 <- as.matrix(SJ_HSI_CAN[, c("mindept", "maxdept", "DEPT_SI")])
sjCAN_DSI_rast <- app(depth_rast_CAN, fun = depth_lookup_fun, data = sjCAN_mat3)
# sjCAN_DSI_rast <- terra::classify(depth_rast_CAN, rcl = sjCAN_mat3, include.lowest = T)

# IGOM
SJ_HSI_IGOM |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SJ_HSI_IGOM

sjIGOM_mat1 <- as.matrix(SJ_HSI_IGOM[, c("mintemp", "maxtemp", "TEMP_SI")])
sjIGOM_TSI_rast <- terra::classify(springTemp_rast_IGOM, rcl = sjIGOM_mat1, include.lowest = T)

sjIGOM_mat2 <- as.matrix(SJ_HSI_IGOM[, c("minsalt", "maxsalt", "SALI_SI")])
sjIGOM_SSI_rast <- terra::classify(springSalt_rast_IGOM, rcl = sjIGOM_mat2, include.lowest = T)

sjIGOM_mat3 <- as.matrix(SJ_HSI_IGOM[, c("mindept", "maxdept", "DEPT_SI")])
sjIGOM_DSI_rast <- app(depth_rast_IGOM, fun = depth_lookup_fun, data = sjIGOM_mat3)
# sjIGOM_DSI_rast <- terra::classify(depth_rast_IGOM, rcl = sjIGOM_mat3, include.lowest = T)

# OGOM
SJ_HSI_OGOM |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SJ_HSI_OGOM

sjOGOM_mat1 <- as.matrix(SJ_HSI_OGOM[, c("mintemp", "maxtemp", "TEMP_SI")])
sjOGOM_TSI_rast <- terra::classify(springTemp_rast_OGOM, rcl = sjOGOM_mat1, include.lowest = T)

sjOGOM_mat2 <- as.matrix(SJ_HSI_OGOM[, c("minsalt", "maxsalt", "SALI_SI")])
sjOGOM_SSI_rast <- terra::classify(springSalt_rast_OGOM, rcl = sjOGOM_mat2, include.lowest = T)

sjOGOM_mat3 <- as.matrix(SJ_HSI_OGOM[, c("mindept", "maxdept", "DEPT_SI")])
sjOGOM_DSI_rast <- app(depth_rast_OGOM, fun = depth_lookup_fun, data = sjOGOM_mat3)
# sjOGOM_DSI_rast <- terra::classify(depth_rast_OGOM, rcl = sjOGOM_mat3, include.lowest = T)

# ISNE
SJ_HSI_ISNE |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SJ_HSI_ISNE

sjISNE_mat1 <- as.matrix(SJ_HSI_ISNE[, c("mintemp", "maxtemp", "TEMP_SI")])
sjISNE_TSI_rast <- terra::classify(springTemp_rast_ISNE, rcl = sjISNE_mat1, include.lowest = T)

sjISNE_mat2 <- as.matrix(SJ_HSI_ISNE[, c("minsalt", "maxsalt", "SALI_SI")])
sjISNE_SSI_rast <- terra::classify(springSalt_rast_ISNE, rcl = sjISNE_mat2, include.lowest = T)

sjISNE_mat3 <- as.matrix(SJ_HSI_ISNE[, c("mindept", "maxdept", "DEPT_SI")])
sjISNE_DSI_rast <- app(depth_rast_ISNE, fun = depth_lookup_fun, data = sjISNE_mat3)
# sjISNE_DSI_rast <- terra::classify(depth_rast_ISNE, rcl = sjISNE_mat3, include.lowest = T)

# OSNE
SJ_HSI_OSNE |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SJ_HSI_OSNE

sjOSNE_mat1 <- as.matrix(SJ_HSI_OSNE[, c("mintemp", "maxtemp", "TEMP_SI")])
sjOSNE_TSI_rast <- terra::classify(springTemp_rast_OSNE, rcl = sjOSNE_mat1, include.lowest = T)

sjOSNE_mat2 <- as.matrix(SJ_HSI_OSNE[, c("minsalt", "maxsalt", "SALI_SI")])
sjOSNE_SSI_rast <- terra::classify(springSalt_rast_OSNE, rcl = sjOSNE_mat2, include.lowest = T)

sjOSNE_mat3 <- as.matrix(SJ_HSI_OSNE[, c("mindept", "maxdept", "DEPT_SI")])
sjOSNE_DSI_rast <- app(depth_rast_OSNE, fun = depth_lookup_fun, data = sjOSNE_mat3)
# sjOSNE_DSI_rast <- terra::classify(depth_rast_OSNE, rcl = sjOSNE_mat3, include.lowest = T)

### SPRING ROCK
# Can
SR_HSI_CAN |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SR_HSI_CAN

srCAN_mat1 <- as.matrix(SR_HSI_CAN[, c("mintemp", "maxtemp", "TEMP_SI")])
srCAN_TSI_rast <- terra::classify(springTemp_rast_CAN, rcl = srCAN_mat1, include.lowest = T)

srCAN_mat2 <- as.matrix(SR_HSI_CAN[, c("minsalt", "maxsalt", "SALI_SI")])
srCAN_SSI_rast <- terra::classify(springSalt_rast_CAN, rcl = srCAN_mat2, include.lowest = T)

srCAN_mat3 <- as.matrix(SR_HSI_CAN[, c("mindept", "maxdept", "DEPT_SI")])
srCAN_DSI_rast <- app(depth_rast_CAN, fun = depth_lookup_fun, data = srCAN_mat3)
# srCAN_DSI_rast <- terra::classify(depth_rast_CAN, rcl = srCAN_mat3, include.lowest = T)

# IGOM
SR_HSI_IGOM |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SR_HSI_IGOM

srIGOM_mat1 <- as.matrix(SR_HSI_IGOM[, c("mintemp", "maxtemp", "TEMP_SI")])
srIGOM_TSI_rast <- terra::classify(springTemp_rast_IGOM, rcl = srIGOM_mat1, include.lowest = T)

srIGOM_mat2 <- as.matrix(SR_HSI_IGOM[, c("minsalt", "maxsalt", "SALI_SI")])
srIGOM_SSI_rast <- terra::classify(springSalt_rast_IGOM, rcl = srIGOM_mat2, include.lowest = T)

srIGOM_mat3 <- as.matrix(SR_HSI_IGOM[, c("mindept", "maxdept", "DEPT_SI")])
srIGOM_DSI_rast <- app(depth_rast_IGOM, fun = depth_lookup_fun, data = srIGOM_mat3)
# srIGOM_DSI_rast <- terra::classify(depth_rast_IGOM, rcl = srIGOM_mat3, include.lowest = T)

# OGOM
SR_HSI_OGOM |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SR_HSI_OGOM

srOGOM_mat1 <- as.matrix(SR_HSI_OGOM[, c("mintemp", "maxtemp", "TEMP_SI")])
srOGOM_TSI_rast <- terra::classify(springTemp_rast_OGOM, rcl = srOGOM_mat1, include.lowest = T)

srOGOM_mat2 <- as.matrix(SR_HSI_OGOM[, c("minsalt", "maxsalt", "SALI_SI")])
srOGOM_SSI_rast <- terra::classify(springSalt_rast_OGOM, rcl = srOGOM_mat2, include.lowest = T)

srOGOM_mat3 <- as.matrix(SR_HSI_OGOM[, c("mindept", "maxdept", "DEPT_SI")])
srOGOM_DSI_rast <- app(depth_rast_OGOM, fun = depth_lookup_fun, data = srOGOM_mat3)
# srOGOM_DSI_rast <- terra::classify(depth_rast_OGOM, rcl = srOGOM_mat3, include.lowest = T)

# ISNE
SR_HSI_ISNE |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SR_HSI_ISNE

srISNE_mat1 <- as.matrix(SR_HSI_ISNE[, c("mintemp", "maxtemp", "TEMP_SI")])
srISNE_TSI_rast <- terra::classify(springTemp_rast_ISNE, rcl = srISNE_mat1, include.lowest = T)

srISNE_mat2 <- as.matrix(SR_HSI_ISNE[, c("minsalt", "maxsalt", "SALI_SI")])
srISNE_SSI_rast <- terra::classify(springSalt_rast_ISNE, rcl = srISNE_mat2, include.lowest = T)

srISNE_mat3 <- as.matrix(SR_HSI_ISNE[, c("mindept", "maxdept", "DEPT_SI")])
srISNE_DSI_rast <- app(depth_rast_ISNE, fun = depth_lookup_fun, data = srISNE_mat3)
# srISNE_DSI_rast <- terra::classify(depth_rast_ISNE, rcl = srISNE_mat3, include.lowest = T)

# OSNE
SR_HSI_OSNE |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> SR_HSI_OSNE

srOSNE_mat1 <- as.matrix(SR_HSI_OSNE[, c("mintemp", "maxtemp", "TEMP_SI")])
srOSNE_TSI_rast <- terra::classify(springTemp_rast_OSNE, rcl = srOSNE_mat1, include.lowest = T)

srOSNE_mat2 <- as.matrix(SR_HSI_OSNE[, c("minsalt", "maxsalt", "SALI_SI")])
srOSNE_SSI_rast <- terra::classify(springSalt_rast_OSNE, rcl = srOSNE_mat2, include.lowest = T)

srOSNE_mat3 <- as.matrix(SR_HSI_OSNE[, c("mindept", "maxdept", "DEPT_SI")])
srOSNE_DSI_rast <- app(depth_rast_OSNE, fun = depth_lookup_fun, data = srOSNE_mat3)
# srOSNE_DSI_rast <- terra::classify(depth_rast_OSNE, rcl = srOSNE_mat3, include.lowest = T)

### FALL ENV RASTERS
fallTemp_rast_ISNE       <- mask(fallTemp_rast_f, ISNE)
fallTemp_rast_OSNE       <- mask(fallTemp_rast_f, OSNE)
fallTemp_rast_IGOM       <- mask(fallTemp_rast_f, IGOM)
fallTemp_rast_OGOM       <- mask(fallTemp_rast_f, OGOM)
fallTemp_rast_CAN       <- mask(fallTemp_rast_f, CAN)

fallSalt_rast_ISNE       <- mask(fallSali_rast_f, ISNE)
fallSalt_rast_OSNE       <- mask(fallSali_rast_f, OSNE)
fallSalt_rast_IGOM       <- mask(fallSali_rast_f, IGOM)
fallSalt_rast_OGOM       <- mask(fallSali_rast_f, OGOM)
fallSalt_rast_CAN       <- mask(fallSali_rast_f, CAN)

## FALL JONAH
# Can
FJ_HSI_CAN |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FJ_HSI_CAN

fjCAN_mat1 <- as.matrix(FJ_HSI_CAN[, c("mintemp", "maxtemp", "TEMP_SI")])
fjCAN_TSI_rast <- terra::classify(fallTemp_rast_CAN, rcl = fjCAN_mat1, include.lowest = T)

fjCAN_mat2 <- as.matrix(FJ_HSI_CAN[, c("minsalt", "maxsalt", "SALI_SI")])
fjCAN_SSI_rast <- terra::classify(fallSalt_rast_CAN, rcl = fjCAN_mat2, include.lowest = T)

fjCAN_mat3 <- as.matrix(FJ_HSI_CAN[, c("mindept", "maxdept", "DEPT_SI")])
fjCAN_DSI_rast <- app(depth_rast_CAN, fun = depth_lookup_fun, data = fjCAN_mat3)
# fjCAN_DSI_rast <- terra::classify(depth_rast_CAN, rcl = fjCAN_mat3, include.lowest = T)

# IGOM
FJ_HSI_IGOM |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FJ_HSI_IGOM

fjIGOM_mat1 <- as.matrix(FJ_HSI_IGOM[, c("mintemp", "maxtemp", "TEMP_SI")])
fjIGOM_TSI_rast <- terra::classify(fallTemp_rast_IGOM, rcl = fjIGOM_mat1, include.lowest = T)

fjIGOM_mat2 <- as.matrix(FJ_HSI_IGOM[, c("minsalt", "maxsalt", "SALI_SI")])
fjIGOM_SSI_rast <- terra::classify(fallSalt_rast_IGOM, rcl = fjIGOM_mat2, include.lowest = T)

fjIGOM_mat3 <- as.matrix(FJ_HSI_IGOM[, c("mindept", "maxdept", "DEPT_SI")])
fjIGOM_DSI_rast <- app(depth_rast_IGOM, fun = depth_lookup_fun, data = fjIGOM_mat3)
#fjIGOM_DSI_rast <- terra::classify(depth_rast_IGOM, rcl = fjIGOM_mat3, include.lowest = T)

# OGOM
FJ_HSI_OGOM |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FJ_HSI_OGOM

fjOGOM_mat1 <- as.matrix(FJ_HSI_OGOM[, c("mintemp", "maxtemp", "TEMP_SI")])
fjOGOM_TSI_rast <- terra::classify(fallTemp_rast_OGOM, rcl = fjOGOM_mat1, include.lowest = T)

fjOGOM_mat2 <- as.matrix(FJ_HSI_OGOM[, c("minsalt", "maxsalt", "SALI_SI")])
fjOGOM_SSI_rast <- terra::classify(fallSalt_rast_OGOM, rcl = fjOGOM_mat2, include.lowest = T)

fjOGOM_mat3 <- as.matrix(FJ_HSI_OGOM[, c("mindept", "maxdept", "DEPT_SI")])
fjOGOM_DSI_rast <- app(depth_rast_OGOM, fun = depth_lookup_fun, data = fjOGOM_mat3)
#fjOGOM_DSI_rast <- terra::classify(depth_rast_OGOM, rcl = fjOGOM_mat3, include.lowest = T)

# ISNE
FJ_HSI_ISNE |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FJ_HSI_ISNE

fjISNE_mat1 <- as.matrix(FJ_HSI_ISNE[, c("mintemp", "maxtemp", "TEMP_SI")])
fjISNE_TSI_rast <- terra::classify(fallTemp_rast_ISNE, rcl = fjISNE_mat1, include.lowest = T)

fjISNE_mat2 <- as.matrix(FJ_HSI_ISNE[, c("minsalt", "maxsalt", "SALI_SI")])
fjISNE_SSI_rast <- terra::classify(fallSalt_rast_ISNE, rcl = fjISNE_mat2, include.lowest = T)

fjISNE_mat3 <- as.matrix(FJ_HSI_ISNE[, c("mindept", "maxdept", "DEPT_SI")])
fjISNE_DSI_rast <- app(depth_rast_ISNE, fun = depth_lookup_fun, data = fjISNE_mat3)
# fjISNE_DSI_rast <- terra::classify(depth_rast_ISNE, rcl = fjISNE_mat3, include.lowest = T)

# OSNE
FJ_HSI_OSNE |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FJ_HSI_OSNE

fjOSNE_mat1 <- as.matrix(FJ_HSI_OSNE[, c("mintemp", "maxtemp", "TEMP_SI")])
fjOSNE_TSI_rast <- terra::classify(fallTemp_rast_OSNE, rcl = fjOSNE_mat1, include.lowest = T)

fjOSNE_mat2 <- as.matrix(FJ_HSI_OSNE[, c("minsalt", "maxsalt", "SALI_SI")])
fjOSNE_SSI_rast <- terra::classify(fallSalt_rast_OSNE, rcl = fjOSNE_mat2, include.lowest = T)

fjOSNE_mat3 <- as.matrix(FJ_HSI_OSNE[, c("mindept", "maxdept", "DEPT_SI")])
#fjOSNE_DSI_rast <- terra::classify(depth_rast_OSNE, rcl = fjOSNE_mat3, include.lowest = T)
fjOSNE_DSI_rast <- app(depth_rast_OSNE, fun = depth_lookup_fun, data = fjOSNE_mat3)

### FaLL ROCK
# Can
FR_HSI_CAN |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FR_HSI_CAN

frCAN_mat1 <- as.matrix(FR_HSI_CAN[, c("mintemp", "maxtemp", "TEMP_SI")])
frCAN_TSI_rast <- terra::classify(fallTemp_rast_CAN, rcl = frCAN_mat1, include.lowest = T)

frCAN_mat2 <- as.matrix(FR_HSI_CAN[, c("minsalt", "maxsalt", "SALI_SI")])
frCAN_SSI_rast <- terra::classify(fallSalt_rast_CAN, rcl = frCAN_mat2, include.lowest = T)

frCAN_mat3 <- as.matrix(FR_HSI_CAN[, c("mindept", "maxdept", "DEPT_SI")])
frCAN_DSI_rast <- app(depth_rast_CAN, fun = depth_lookup_fun, data = frCAN_mat3)
# frCAN_DSI_rast <- terra::classify(depth_rast_CAN, rcl = frCAN_mat3, include.lowest = T)

# IGOM
FR_HSI_IGOM |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FR_HSI_IGOM

frIGOM_mat1 <- as.matrix(FR_HSI_IGOM[, c("mintemp", "maxtemp", "TEMP_SI")])
frIGOM_TSI_rast <- terra::classify(fallTemp_rast_IGOM, rcl = frIGOM_mat1, include.lowest = T)

frIGOM_mat2 <- as.matrix(FR_HSI_IGOM[, c("minsalt", "maxsalt", "SALI_SI")])
frIGOM_SSI_rast <- terra::classify(fallSalt_rast_IGOM, rcl = frIGOM_mat2, include.lowest = T)

frIGOM_mat3 <- as.matrix(FR_HSI_IGOM[, c("mindept", "maxdept", "DEPT_SI")])
frIGOM_DSI_rast <- app(depth_rast_IGOM, fun = depth_lookup_fun, data = fjIGOM_mat3)
# frIGOM_DSI_rast <- terra::classify(depth_rast_IGOM, rcl = frIGOM_mat3, include.lowest = T)

# OGOM
FR_HSI_OGOM |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FR_HSI_OGOM

frOGOM_mat1 <- as.matrix(FR_HSI_OGOM[, c("mintemp", "maxtemp", "TEMP_SI")])
frOGOM_TSI_rast <- terra::classify(fallTemp_rast_OGOM, rcl = frOGOM_mat1, include.lowest = T)

frOGOM_mat2 <- as.matrix(FR_HSI_OGOM[, c("minsalt", "maxsalt", "SALI_SI")])
frOGOM_SSI_rast <- terra::classify(fallSalt_rast_OGOM, rcl = frOGOM_mat2, include.lowest = T)

frOGOM_mat3 <- as.matrix(FR_HSI_OGOM[, c("mindept", "maxdept", "DEPT_SI")])
frOGOM_DSI_rast <- app(depth_rast_OGOM, fun = depth_lookup_fun, data = fjOGOM_mat3)
# frOGOM_DSI_rast <- terra::classify(depth_rast_OGOM, rcl = frOGOM_mat3, include.lowest = T)

# ISNE
FR_HSI_ISNE |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FR_HSI_ISNE

frISNE_mat1 <- as.matrix(FR_HSI_ISNE[, c("mintemp", "maxtemp", "TEMP_SI")])
frISNE_TSI_rast <- terra::classify(fallTemp_rast_ISNE, rcl = frISNE_mat1, include.lowest = T)

frISNE_mat2 <- as.matrix(FR_HSI_ISNE[, c("minsalt", "maxsalt", "SALI_SI")])
frISNE_SSI_rast <- terra::classify(fallSalt_rast_ISNE, rcl = frISNE_mat2, include.lowest = T)

frISNE_mat3 <- as.matrix(FR_HSI_ISNE[, c("mindept", "maxdept", "DEPT_SI")])
frISNE_DSI_rast <- app(depth_rast_ISNE, fun = depth_lookup_fun, data = fjISNE_mat3)
# frISNE_DSI_rast <- terra::classify(depth_rast_ISNE, rcl = frISNE_mat3, include.lowest = T)

# OSNE
FR_HSI_OSNE |> 
  group_by(row_number()) |>
  mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
         minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
         mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> FR_HSI_OSNE

frOSNE_mat1 <- as.matrix(FR_HSI_OSNE[, c("mintemp", "maxtemp", "TEMP_SI")])
frOSNE_TSI_rast <- terra::classify(fallTemp_rast_OSNE, rcl = frOSNE_mat1, include.lowest = T)

frOSNE_mat2 <- as.matrix(FR_HSI_OSNE[, c("minsalt", "maxsalt", "SALI_SI")])
frOSNE_SSI_rast <- terra::classify(fallSalt_rast_OSNE, rcl = frOSNE_mat2, include.lowest = T)

frOSNE_mat3 <- as.matrix(FR_HSI_OSNE[, c("mindept", "maxdept", "DEPT_SI")])
frOSNE_DSI_rast <- app(depth_rast_OSNE, fun = depth_lookup_fun, data = fjOSNE_mat3)
# frOSNE_DSI_rast <- terra::classify(depth_rast_OSNE, rcl = frOSNE_mat3, include.lowest = T)

### RECOMBINE RASTERS
# Spring Jonah
SJ_TSI <- mosaic(sjCAN_TSI_rast, sjIGOM_TSI_rast, sjOGOM_TSI_rast, sjISNE_TSI_rast, sjOSNE_TSI_rast)
SJ_SSI <- mosaic(sjCAN_SSI_rast, sjIGOM_SSI_rast, sjOGOM_SSI_rast, sjISNE_SSI_rast, sjOSNE_SSI_rast)
SJ_DSI <- mosaic(sjCAN_DSI_rast, sjIGOM_DSI_rast, sjOGOM_DSI_rast, sjISNE_DSI_rast, sjOSNE_DSI_rast)
SJ_DSI2 <- resample(SJ_DSI, SJ_TSI)

u_SJ_TSI_9211 <- terra::app(SJ_TSI[[1:19]], fun = mean, na.rm = F) # 1992 - 2011
u_SJ_TSI_1224 <- terra::app(SJ_TSI[[20:32]], fun = mean, na.rm = F) # 2012 - 2024

u_SJ_SSI_9211 <- terra::app(SJ_SSI[[1:19]], fun = mean, na.rm = F) # 1992 - 2011
u_SJ_SSI_1224 <- terra::app(SJ_SSI[[20:32]], fun = mean, na.rm = F) # 2012 - 2024

SJSI_stack_9211 <- c(u_SJ_TSI_9211, u_SJ_SSI_9211, SJ_DSI2)
SJSI_stack_1224 <- c(u_SJ_TSI_1224, u_SJ_SSI_1224, SJ_DSI2)

SJ_HSI_9211 <- terra::app(SJSI_stack_9211, fun = mean, na.rm = F) # 1992 - 2011
SJ_HSI_1224 <- terra::app(SJSI_stack_1224, fun = mean, na.rm = F) # 1992 - 2011

# Spring Rock
SR_TSI <- mosaic(srCAN_TSI_rast, srIGOM_TSI_rast, srOGOM_TSI_rast, srISNE_TSI_rast, srOSNE_TSI_rast)
SR_SSI <- mosaic(srCAN_SSI_rast, srIGOM_SSI_rast, srOGOM_SSI_rast, srISNE_SSI_rast, srOSNE_SSI_rast)
SR_DSI <- mosaic(srCAN_DSI_rast, srIGOM_DSI_rast, srOGOM_DSI_rast, srISNE_DSI_rast, srOSNE_DSI_rast)
SR_DSI2 <- resample(SR_DSI, SR_TSI)

u_SR_TSI_9211 <- terra::app(SR_TSI[[1:19]], fun = mean, na.rm = F) # 1992 - 2011
u_SR_TSI_1224 <- terra::app(SR_TSI[[20:32]], fun = mean, na.rm = F) # 2012 - 2024

u_SR_SSI_9211 <- terra::app(SR_SSI[[1:19]], fun = mean, na.rm = F) # 1992 - 2011
u_SR_SSI_1224 <- terra::app(SR_SSI[[20:32]], fun = mean, na.rm = F) # 2012 - 2024

SRSI_stack_9211 <- c(u_SR_TSI_9211, u_SR_SSI_9211, SR_DSI2)
SRSI_stack_1224 <- c(u_SR_TSI_1224, u_SR_SSI_1224, SR_DSI2)

SR_HSI_9211 <- terra::app(SRSI_stack_9211, fun = mean, na.rm = F) # 1992 - 2011
SR_HSI_1224 <- terra::app(SRSI_stack_1224, fun = mean, na.rm = F) # 2012 - 2024

# Fall Jonah
FJ_TSI <- mosaic(fjCAN_TSI_rast, fjIGOM_TSI_rast, fjOGOM_TSI_rast, fjISNE_TSI_rast, fjOSNE_TSI_rast)
FJ_SSI <- mosaic(fjCAN_SSI_rast, fjIGOM_SSI_rast, fjOGOM_SSI_rast, fjISNE_SSI_rast, fjOSNE_SSI_rast)
FJ_DSI <- mosaic(fjCAN_DSI_rast, fjIGOM_DSI_rast, fjOGOM_DSI_rast, fjISNE_DSI_rast, fjOSNE_DSI_rast)
FJ_DSI2 <- resample(FJ_DSI, FJ_TSI)

u_FJ_TSI_9211 <- terra::app(FJ_TSI[[1:19]], fun = mean, na.rm = F) # 1992 - 2011
u_FJ_TSI_1224 <- terra::app(FJ_TSI[[20:32]], fun = mean, na.rm = F) # 2012 - 2024

u_FJ_SSI_9211 <- terra::app(FJ_SSI[[1:19]], fun = mean, na.rm = F) # 1992 - 2011
u_FJ_SSI_1224 <- terra::app(FJ_SSI[[20:32]], fun = mean, na.rm = F) # 2012 - 2024

FJSI_stack_9211 <- c(u_FJ_TSI_9211, u_FJ_SSI_9211, FJ_DSI2)
FJSI_stack_1224 <- c(u_FJ_TSI_1224, u_FJ_SSI_1224, FJ_DSI2)

FJ_HSI_9211 <- terra::app(FJSI_stack_9211, fun = mean, na.rm = F) # 1992 - 2011
FJ_HSI_1224 <- terra::app(FJSI_stack_1224, fun = mean, na.rm = F) # 2012 - 2024

# Fall Rock 
FR_TSI <- mosaic(frCAN_TSI_rast, frIGOM_TSI_rast, frOGOM_TSI_rast, frISNE_TSI_rast, frOSNE_TSI_rast)
FR_SSI <- mosaic(frCAN_SSI_rast, frIGOM_SSI_rast, frOGOM_SSI_rast, frISNE_SSI_rast, frOSNE_SSI_rast)
FR_DSI <- mosaic(frCAN_DSI_rast, frIGOM_DSI_rast, frOGOM_DSI_rast, frISNE_DSI_rast, frOSNE_DSI_rast)
FR_DSI2 <- resample(FR_DSI, FR_TSI)

u_FR_TSI_9211 <- terra::app(FR_TSI[[1:19]], fun = mean, na.rm = F) # 1992 - 2011
u_FR_TSI_1224 <- terra::app(FR_TSI[[20:32]], fun = mean, na.rm = F) # 2012 - 2024

u_FR_SSI_9211 <- terra::app(FR_SSI[[1:19]], fun = mean, na.rm = F) # 1992 - 2011
u_FR_SSI_1224 <- terra::app(FR_SSI[[20:32]], fun = mean, na.rm = F) # 2012 - 2024

FRSI_stack_9211 <- c(u_FR_TSI_9211, u_FR_SSI_9211, FR_DSI2)
FRSI_stack_1224 <- c(u_FR_TSI_1224, u_FR_SSI_1224, FR_DSI2)

FR_HSI_9211 <- terra::app(FRSI_stack_9211, fun = mean, na.rm = F) # 1992 - 2011
FR_HSI_1224 <- terra::app(FRSI_stack_1224, fun = mean, na.rm = F) # 2012 - 2024

################ MAPPING ##################
## jonah spring temp SI map
ggplot() +
  geom_raster(data = u_SJ_TSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Temperature",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sjTSImap

ggplot() +
  geom_raster(data = u_SJ_SSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Salinity",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sjSSImap

ggplot() +
  geom_raster(data = SJ_DSI2, aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Depth",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sjDSImap

((sjTSImap + sjSSImap + sjDSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Spring Jonah Crab Suitability Indices 1992 - 2011",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> sjSI_9211

ggsave(filename = "springJonah_SI9211.png",
       plot = sjSI_9211,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 9,
       width = 22)


# spring jonah 2012 
# jonah spring temp SI map
ggplot() +
  geom_raster(data = u_SJ_TSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Temperature",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sjTSImap

ggplot() +
  geom_raster(data = u_SJ_SSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Salinity",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sjSSImap

ggplot() +
  geom_raster(data = SJ_DSI2, aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Depth",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sjDSImap


((sjTSImap + sjSSImap + sjDSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Spring Jonah Crab Suitability Indices 2012 - 2024",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> sjSI_1224

ggsave(filename = "springJonah_1224.png",
       plot = sjSI_1224,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 9,
       width = 22)

########## SPRING ROCK CRAB SI PLOTS #########################
####### Time window 1 ##############
ggplot() +
  geom_raster(data = u_SR_TSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Temperature",
       y = "Latitude",
       x = "Longtiude",
       fill = "SI") +
  theme_classic(base_size = 17) -> srTSImap

ggplot() +
  geom_raster(data = u_SR_SSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Salinity",
       y = "",
       x = "Longtiude",
       fill = "SI") +
  theme_classic(base_size = 17) -> srSSImap

ggplot() +
  geom_raster(data = SR_DSI2, aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Depth",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> srDSImap


((srTSImap + srSSImap + srDSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Spring Atl. Rock Crab Suitability Indices 1992 - 2012",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> srSI_9211

ggsave(filename = "springRock_9211.png",
       plot = srSI_9211,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 9,
       width = 22)


####### Time window 2 ##############
ggplot() +
  geom_raster(data = u_SR_TSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Temperature",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> srTSImap

ggplot() +
  geom_raster(data = u_SR_SSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Salinity",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> srSSImap

ggplot() +
  geom_raster(data = SR_DSI2, aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Depth",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> srDSImap


((srTSImap + srSSImap + srDSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Spring Atl. Rock Crab Suitability Indices 2012 - 2024",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> srHSI_1224

ggsave(filename = "springRock_1224.png",
       plot = srHSI_1224,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 9,
       width = 22)


############## FALL JONAH ###############
####### TIME WINDOW 1 ###################
ggplot() +
  geom_raster(data = u_FJ_TSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Temperature",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> fjTSImap

ggplot() +
  geom_raster(data = u_FJ_SSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Salinity",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> fjSSImap

ggplot() +
  geom_raster(data = FJ_DSI2, aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Depth",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> fjDSImap


((fjTSImap + fjSSImap + fjDSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Fall Jonah Crab Suitability Indices 1992 - 2011",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> fjHSI_9211

ggsave(filename = "fallJonah_9211.png",
       plot = fjHSI_9211,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 9,
       width = 22)


# fall jonah 2012 
ggplot() +
  geom_raster(data = u_FJ_TSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Temperature",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> fjTSImap_1224

ggplot() +
  geom_raster(data = u_FJ_SSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Salinity",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> fjSSImap_1224

ggplot() +
  geom_raster(data = FJ_DSI2, aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Depth",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> fjDSImap

((fjTSImap_1224 + fjSSImap_1224 + fjDSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Fall Jonah Crab Suitability Indices 2012 - 2024",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> fjHSI_1224

ggsave(filename = "fallJonah_1224.png",
       plot = fjHSI_1224,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 9,
       width = 22)

########## Fall ROCK CRAB SI PLOTS #########################
####### Time window 1 ##############
ggplot() +
  geom_raster(data = u_FR_TSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Temperature",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> frTSImap

ggplot() +
  geom_raster(data = u_FR_SSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Salinity",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> frSSImap

ggplot() +
  geom_raster(data = FR_DSI2, aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Depth",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> frDSImap

((frTSImap + frSSImap + frDSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Fall Atl. Rock Crab Suitability Indices 1992 - 2012",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> frHSI_9211

ggsave(filename = "fallRock_9211.png",
       plot = frHSI_9211,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 9,
       width = 22)

####### Time window 2 ##############
ggplot() +
  geom_raster(data = u_FR_TSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Temperature",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> frTSImap

ggplot() +
  geom_raster(data = u_FR_SSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Salinity",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> frSSImap

ggplot() +
  geom_raster(data = FR_DSI2, aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Depth",
       y = "",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> frDSImap

((frTSImap + frSSImap + frDSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Fall Atl. Rock Crab Suitability Indices 2012 - 2024",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> frHSI_1224

ggsave(filename = "fallRock_1224.png",
       plot = frHSI_1224,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 9,
       width = 22)

############## HSI MAPPING ######

## JONAH CRAB HSI ####

ggplot() +
  geom_raster(data = SJ_HSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring, 1992-2011",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sjHSImap_9211

ggplot() +
  geom_raster(data = SJ_HSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring, 2012-2024",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sjHSImap_1224

ggplot() +
  geom_raster(data = FJ_HSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Fall, 1992-2011",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> fjHSImap_9211

ggplot() +
  geom_raster(data = SJ_HSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Fall, 2012-2024",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> fjHSImap_1224


((sjHSImap_9211 + sjHSImap_1224) / (fjHSImap_9211 + fjHSImap_1224)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Jonah Crab Habitat Suitability",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> jonahHSI

ggsave(filename = "JonahHSImaps.png",
       plot = jonahHSI,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 14,
       width = 16)

######### ROCK CRAB HSI ###########

ggplot() +
  geom_raster(data = SR_HSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring, 1992-2011",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> srHSImap_9211

ggplot() +
  geom_raster(data = SR_HSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring, 2012-2024",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> srHSImap_1224

ggplot() +
  geom_raster(data = FR_HSI_9211, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Fall, 1992-2011",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> frHSImap_9211

ggplot() +
  geom_raster(data = FR_HSI_1224, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) + 
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Fall, 2012-2024",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> frHSImap_1224


((srHSImap_9211 + srHSImap_1224) / (frHSImap_9211 + frHSImap_1224)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Atl. Rock Crab Habitat Suitability",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> rockHSI

ggsave(filename = "RockHSImaps.png",
       plot = rockHSI,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 14,
       width = 16)
