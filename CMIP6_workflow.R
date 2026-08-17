# Load required libraries
library(ncdf4)
library(tidyterra)
library(raster)
library(fields)   # for image.plot()

# data downloaded from https://esgf-metagrid.cloud.dkrz.de/search

# Jonah crab stocks
NOAA_stat_areas <-st_read("~/Documents/R Projects/Shapefiles/NEFSC_GIS/Statistical_Areas_2010.shp") |>
  dplyr::filter(Id > 460 & Id < 641) 

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

ISNE <- st_as_sf(NOAA_stat_areas |> filter(Stock == "ISNE"))
OSNE <- st_as_sf(NOAA_stat_areas |> filter(Stock == "OSNE"))
IGOM <- st_as_sf(NOAA_stat_areas |> filter(Stock == "IGOM"))
OGOM <- st_as_sf(NOAA_stat_areas |> filter(Stock == "OGOM"))
CAN <- st_as_sf(NOAA_stat_areas |> filter(Stock == "CAN"))


# FUNCTIONS
# Create a function to assign suitability based on ranges
depth_lookup_fun <- function(x, data) {
  sapply(x, function(val) {
    matched <- data[val >= data[,1] & val <= data[,2], 3]
    if (length(matched) > 0) matched[1] else NA  # NA for out-of-range depths
  })
}

rast_project_HSI <- function(HSI, 
                             rast1,  # Temperature
                             rast2,  # Salinity
                             rast3) # Depth
{ 
  HSI |> 
    group_by(row_number()) |>
    mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
           minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
           mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> HSI
  
  # Temperature 
  TSI.mat <- as.matrix(HSI[, c("mintemp", "maxtemp", "TEMP_SI")])
  rast1a <- terra::classify(rast1, rcl = TSI.mat, include.lowest = T)
  rast1a[ rast1 < 0 | rast1 > 1 ] <- NA  
  
  # Salinity 
  SSI.mat <- as.matrix(HSI[, c("minsalt", "maxsalt", "SALI_SI")])
  rast2a <- terra::classify(rast2, rcl = SSI.mat, include.lowest = T)
  rast2a[ rast2 < 0 | rast2 > 1 ] <- NA
  
  # Depth
  DSI <- as.matrix(HSI[, c("mindept", "maxdept", "DEPT_SI")])
  rast3a <- app(rast3, fun = depth_lookup_fun, data = DSI)
  
  # Stack
  rast_stack <- c(rast1a, rast2a, rast3a)
  rast_out <- terra::app(rast_stack, fun = mean, na.rm = T)
  
  return(rast_out)
}

rasterize_2d_downscale_CMIP6 <- function(ncfile, # string with path to netcdf
                                         project_crs, # string describing the CRS used by project files
                                         year1, # filter value, "2075"
                                         year2, # filter value 2, "2094"
                                         lon_call = "longitude",
                                         lat_call = "latitude",
                                         months, # filter value 3 c("03", "04", "05")
                                         mask_poly, # polygon to mask raster to area of interest
                                         coords, # coordinate grid for downscaling
                                         name, # name for coords downscale 
                                         t1, 
                                         t2,
                                         shift_x,
                                         shift_y,
                                         flip_it
){
  ncfile1 <- ncfile
  nc <- nc_open(ncfile1)
  lon <- ncvar_get(nc, lon_call)
  lat <- ncvar_get(nc, lat_call)
  nc_close(nc)
  
  # rasterize data
  r <- rast(ncfile)
  
  if (flip_it == T) {
    r <- terra::flip(r, direction = "vertical")
  }
  
  
  # reproject with new extent 
  ext(r) <- ext(min(lon), max(lon), min(lat), max(lat)) 
  
  # assign project CRS
  crs(r) <- project_crs
  
  # extract values 
  yrs <- format(time(r), "%Y")
  mos <- format(time(r), "%m")
  
  sub_r <- r[[ yrs >= year1 & yrs <= year2 & mos %in% months ]]
  
  # mean spring 
  mean_sub_r <- terra::app(sub_r, fun = mean, na.rm = F)
  # plot(mean_sub_r)
  
  # put in correct orientation if necessary
  # if (max(lon) > 180) {
  #   mean_sub_r <- rotate(mean_sub_r)
  # }
  mean_sub_r <- shift(mean_sub_r, dx = shift_x, shift_y)
  plot(mean_sub_r)
  plot(Strata, add = TRUE, col = "red")
  abline(v = -67)
  
  plot(mean_sub_r, ylim = c(30,50), xlim = c(-60, -75))
  plot(Strata, add = TRUE, col = "red")
  abline(v = -67)
  
  # plot(mean_sub_r)
  # abline(v = 0)
  # plot(Strata, add = TRUE, col = "red")
  
  # mask
  mean_sub_r <- terra::mask(mean_sub_r, mask_poly)
  
  # make data point for downscaling 
  mean_sub_r.pt <- as.data.frame(mean_sub_r, xy = TRUE, na.rm = TRUE)
  
  # downscale model
  mean_sub_r_m1 <-             gam(mean ~ s(x, y, bs = "tp"),
                                   data = mean_sub_r.pt,
                                   method = "REML")
  
  # downscale predictions
  coords$name <- predict(mean_sub_r_m1, newdata = coords)
  coords <- coords[!is.na(coords$name), ]
  
  coords_rast <- terra::rast(coords, type = "xyz", crs = "EPSG:4269")
  coords_rast <- terra::mask(coords_rast, Strata_2)
  
  
  
  
  (ggplot() +
      geom_raster(data = mean_sub_r, aes(x = x, y = y, fill = mean)) +
      scale_fill_distiller(palette = "Spectral", na.value = "grey") +
      geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey", alpha = 0.1) +
      geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey", alpha = 0.1) +
      coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
      labs(#title = t1,
        y = "Latitude",
        x = "Longitude",
        fill = "") +
      theme_classic(base_size = 17) -> p1)
  
  
  (ggplot() +
      geom_raster(data = coords_rast, aes(x = x, y = y, fill = name)) +
      scale_fill_distiller(palette = "Spectral", na.value = "grey") +
      geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
      geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
      coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
      labs(#title = t2,
        y = "Latitude",
        x = "Longitude",
        fill = "") +
      theme_classic(base_size = 17) -> p2)
  
  
  return(list(p1 + p2,
              coords_rast))
  
}

extract_bottom_salinity <- function(
    salinity_nc, # string with path to netcdf
    depth_nc, # string to depth netcdf --- from the same model grid 
    project_crs, # string describing the CRS used by project files
    lon_call = "longitude",
    lat_call = "latitude",
    year1, # filter value, "2075"
    year2, # filter value 2, "2094"
    months, # filter value 3 c("03", "04", "05")
    mask_poly, # polygon to mask raster to area of interest
    coords, # coordinate grid for downscaling
    name, # name for coords downscale 
    t1, 
    t2,
    shift_x,
    shift_y,
    flip_it
){
  # load salinity data to extract values and levels
  ncfile1 <- salinity_nc
  nc <- nc_open(ncfile1)
  lon <- ncvar_get(nc, lon_call)
  lat <- ncvar_get(nc, lat_call)
  lev <- ncvar_get(nc, "lev")   # numeric depth levels (m)
  nc_close(nc)
  
  # salinity raster to filter and overlay
  so_r <- rast(salinity_nc, subds = "so")
  plot(so_r[[1]])
  so_r
  
  # depth raster to extract bottom values 
  depth_o <- rast(depth_nc, subds = "deptho")
  plot(depth_o[[1]])
  depth_o
  
  # filter salinity raster to relevant dates for meaning 
  yrs <- format(time(so_r), "%Y")
  mos <- format(time(so_r), "%m")
  
  sub_r <- so_r[[ yrs >= year1 & yrs <= year2 & mos %in% months ]]
  
  # compute bottom level index per cell (1..length(lev))
  lev <- sort(lev)  # make sure levels are ascending
  
  k <- app(depth_o, fun = function(d) {
    out <- rep(NA_integer_, length(d))
    ok  <- !is.na(d)
    
    out[ok] <- findInterval(d[ok], lev)
    
    # findInterval returns 0 if d < min(lev)
    out[out == 0] <- NA_integer_
    out
  })
  
  plot(k)
  
  # then for each time, pick the layer corresponding to (time, bottom lev)
  # You need to know how your layers are ordered: usually lev varies fastest within time.
  nlev  <- length(lev)
  nlay  <- nlyr(so_r)
  ntime <- nlay / nlev
  
  # benthic <- rast(so_r[[1]])  # template
  # benthic
  # plot(benthic[[1]])
  benthic_list <- vector("list", ntime)
  
  for (t in 1:ntime) {
    
    lev_layers_t <- so_r[[ ((t - 1) * nlev + 1) : (t * nlev) ]]
    
    # combine salinity levels + index raster
    stack_t <- c(lev_layers_t, k)
    
    benthic_list[[t]] <- app(
      stack_t,
      fun = function(v) {
        kk <- v[nlev + 1]
        if (is.na(kk)) return(NA_real_)
        v[kk]
      }
    )
  }
  
  rm(so_r)
  
  mean_benthic <- terra::app(rast(benthic_list), fun = mean, na.rm = F)
  
  if (flip_it == T) {
    mean_benthic <- terra::flip(mean_benthic, direction = "vertical")
  }
  
  
  # reproject with new extent 
  ext(mean_benthic) <- ext(min(lon), max(lon), min(lat), max(lat)) 
  
  # assign project CRS
  crs(mean_benthic) <- project_crs
  
  r3 <- shift(mean_benthic, dx = shift_x, shift_y)
  plot(r3)
  plot(Strata, add = TRUE, col = "red")
  abline(v = -67)
  
  plot(r3, ylim = c(30,50), xlim = c(-60, -75))
  plot(Strata, add = TRUE, col = "red")
  abline(v = -67)
  
  # mask
  masked_rast <- terra::mask(r3, mask_poly)
  
  # make data point for downscaling 
  masked_rast.pt <- as.data.frame(masked_rast, xy = TRUE, na.rm = TRUE)
  
  # downscale model
  masked_rast_m1 <-             gam(mean ~ s(x, y, bs = "tp"),
                                    data = masked_rast.pt,
                                    method = "REML")
  
  # downscale predictions
  coords$name <- predict(masked_rast_m1, newdata = coords)
  coords <- coords[!is.na(coords$name), ]
  
  coords_rast <- terra::rast(coords, type = "xyz", crs = "EPSG:4269")
  coords_rast <- terra::mask(coords_rast, Strata_2)
  
  
  (ggplot() +
      geom_raster(data = masked_rast, aes(x = x, y = y, fill = mean)) +
      scale_fill_distiller(palette = "Spectral", na.value = "grey") +
      geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey", alpha = 0.1) +
      geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey", alpha = 0.1) +
      coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
      labs(#title = t1,
        y = "Latitude",
        x = "Longitude",
        fill = "") +
      theme_classic(base_size = 17) -> p1)
  
  
  (ggplot() +
      geom_raster(data = coords_rast, aes(x = x, y = y, fill = name)) +
      scale_fill_distiller(palette = "Spectral", na.value = "grey") +
      geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
      geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
      coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
      labs(#title = t2,
        y = "Latitude",
        x = "Longitude",
        fill = "") +
      theme_classic(base_size = 17) -> p2)
  
  
  return(list(p1 + p2,
              coords_rast))
}


#### sanity check code ####
# CanESM5
ncfile <- "~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_CanESM5_ssp585_r10i1p1f1_gn_201501-210012.nc"
nc <- nc_open(ncfile)
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
nc_close(nc)

c(min(lon), max(lon), min(lat), max(lat))

r <- rast("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_CanESM5_ssp585_r10i1p1f1_gn_201501-210012.nc")
plot(r[[100]])
r
plot(Strata, add = TRUE, col = "red")

if (lat[1] > lat[length(lat)]) {
  r2 <- terra::flip(r, direction = "vertical")
}

plot(r2[[100]])
r2
plot(Strata, add = TRUE, col = "red")
abline(v = 0)



r3 <- shift(r2, dx = -286)
plot(r3[[100]])
plot(Strata, add = TRUE, col = "red")
abline(v = -67)

# UKESM1
ncfile <- "~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_UKESM1-0-LL_ssp585_r2i1p1f2_gn_205001-210012.nc"
nc <- nc_open(ncfile)
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
nc_close(nc)

c(min(lon), max(lon), min(lat), max(lat))

r <- rast("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_UKESM1-0-LL_ssp585_r2i1p1f2_gn_205001-210012.nc")
plot(r[[100]])
r
plot(Strata, add = TRUE, col = "red")

if (lat[1] > lat[length(lat)]) {
  r2 <- terra::flip(r, direction = "vertical")
}

ext(r2) <- ext(min(lon), max(lon), min(lat), max(lat))
plot(r2[[100]])


r3 <- shift(r2, dx = -107, dy = -4)
plot(r3[[100]])
plot(Strata, add = TRUE, col = "red")
abline(v = -67)

plot(r3[[100]], ylim = c(30,50), xlim = c(-60, -75))
plot(Strata, add = TRUE, col = "red")
abline(v = -67)


# MPI-ESM1-2-LR
ncfile <- "~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_MPI-ESM1-2-LR_ssp585_r28i1p1f1_gn_207501-209412.nc"
nc <- nc_open(ncfile)
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
nc_close(nc)

c(min(lon), max(lon), min(lat), max(lat))

r <- rast("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_MPI-ESM1-2-LR_ssp585_r28i1p1f1_gn_207501-209412.nc")
plot(r[[100]])
r
plot(Strata, add = TRUE, col = "red")

ext(r) <- ext(min(lon), max(lon), min(lat), max(lat)) 
plot(r[[100]])
# r2 <- rotate(r)
# plot(r2[[100]])

r4 <- shift(r, dx = -185, dy = 16)
plot(r3[[100]])
plot(Strata, add = TRUE, col = "red")
abline(v = -67)

plot(r4[[100]], ylim = c(30,50), xlim = c(-60, -75))
plot(Strata, add = TRUE, col = "red")
abline(v = -67)

# CNRM-CM6-1
ncfile <- "~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_CNRM-CM6-1_ssp585_r1i1p1f2_gn_201501-210012.nc"
nc <- nc_open(ncfile)
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

c(min(lon), max(lon), min(lat), max(lat))

r <- rast("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_CNRM-CM6-1_ssp585_r1i1p1f2_gn_201501-210012.nc")
plot(r[[100]])
r
plot(Strata, add = TRUE, col = "red")

ext(r) <- ext(min(lon), max(lon), min(lat), max(lat)) 
plot(r[[100]])


r2 <- terra::flip(r, direction = "vertical")
plot(r2[[100]])

r3 <- shift(r2, dx = -107, dy = 0)
plot(r3[[100]])
plot(Strata, add = TRUE, col = "red")
abline(v = -67)

plot(r3[[100]], ylim = c(30,50), xlim = c(-60, -75))
plot(Strata, add = TRUE, col = "red")
abline(v = -67)


############# TEMPERATURE PROJECTIONS / DOWNSCALING ###########################

countries <- map_data("world")
states <- map_data("state")

# prediction grid 
coords <- xyFromCell(clim_salt_S, 1:ncell(clim_salt_S)) |> as.data.frame()
 

# TEMPERATURE MODELING
# MODEL CNRM-CM6-1, SPRING, TEMPERATURE 
Temp_S_CNRM_CM6_1_down <- rasterize_2d_downscale_CMIP6(ncfile = "~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_CNRM-CM6-1_ssp585_r1i1p1f2_gn_201501-210012.nc",
                                                       project_crs = "EPSG:4269",
                                                       year1 = "2075",
                                                       year2 = "2094",
                                                       lon_call = "lon",
                                                       lat_call = "lat",
                                                       months = c("03", "04", "05"),
                                                       mask_poly = Strata_2,
                                                       coords = coords,
                                                       name = name,
                                                       t1 = "Temp_S_CNRM-CM6-1",
                                                       t2 = "Temp_S_CNRM-CM6-1 downscale",
                                                       shift_x = -106,
                                                       shift_y = 0,
                                                    flip_it = T)

Temp_S_CNRM_CM6_1_down

# MODEL CNRM-CM6-1, FALL, TEMPERATURE
Temp_F_CNRM_CM6_1_down <- rasterize_2d_downscale_CMIP6( "~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_CNRM-CM6-1_ssp585_r1i1p1f2_gn_201501-210012.nc",
                                                    "EPSG:4269",
                                                    "2075",
                                                    "2094",
                                                    lon_call = "lon",
                                                    lat_call = "lat",
                                                    c("09", "10", "11"),
                                                    Strata_2,
                                                    coords,
                                                    Temp_F_CNRM_CM6_1,
                                                    "Temp_F_CNRM-CM6-1",
                                                    "Temp_F_CNRM-CM6-1 downscale",
                                                    shift_x = -106,
                                                    shift_y = 0,
                                                    flip_it = T)

Temp_F_CNRM_CM6_1_down

# MODEL CanESM5, SPRING, TEMPERATURE 
Temp_S_CanESM5_down <- rasterize_2d_downscale_CMIP6(ncfile = "~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_CanESM5_ssp585_r10i1p1f1_gn_201501-210012.nc",
                                                    project_crs = "EPSG:4269",
                                                    year1 = "2075",
                                                    year2 = "2094",
                                                    months = c("03", "04", "05"),
                                                    mask_poly = Strata_2,
                                                    coords = coords,
                                                    name = name,
                                                    t1 = "Temp_S_CanESM5",
                                                    t2 = "Temp_S_CanESM5 downscale",
                                                    shift_x = -286,
                                                    shift_y = 0,
                                                    flip_it = T)

Temp_S_CanESM5_down

# MODEL CanESM5, FALL, TEMPERATURE 
Temp_F_CanESM5_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_CanESM5_ssp585_r10i1p1f1_gn_201501-210012.nc",
                                                    "EPSG:4269",
                                                    lon_call = "longitude",
                                                    lat_call = "latitude",
                                                    "2075",
                                                    "2094",
                                                    c("09", "10", "11"),
                                                    Strata_2,
                                                    coords,
                                                    Temp_F_CanESM5,
                                                    "Temp_F_CanESM5",
                                                    "Temp_F_CanESM5 downscale",
                                                    shift_x = -286,
                                                    shift_y = 0,
                                                    flip_it = T)

Temp_F_CanESM5_down


# MODEL MPI_ESM1_2, SPRING, TEMPERATURE
Temp_S_MPI_ESM1_2_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_MPI-ESM1-2-LR_ssp585_r28i1p1f1_gn_207501-209412.nc",
                                                       "EPSG:4269",
                                                       lon_call = "longitude",
                                                       lat_call = "latitude",
                                                       "2075",
                                                       "2094",
                                                       c("03", "04", "05"),
                                                       Strata_2,
                                                       coords,
                                                       Temp_S_MPI_ESM1_2,
                                                       "Temp_S_MPI_ESM1_2",
                                                       "Temp_S_MPI_ESM1_2 downscale",
                                                       shift_x = -185,
                                                       shift_y = 16,
                                                       flip_it = F)

Temp_S_MPI_ESM1_2_down

# MODEL MPI_ESM1_2, FALL, TEMPERATURE
Temp_F_MPI_ESM1_2_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_MPI-ESM1-2-LR_ssp585_r28i1p1f1_gn_207501-209412.nc",
                                                       "EPSG:4269",
                                                       lon_call = "longitude",
                                                       lat_call = "latitude",
                                                       "2075",
                                                       "2094",
                                                       c("09", "10", "11"),
                                                       Strata_2,
                                                       coords,
                                                       Temp_F_MPI_ESM1_2,
                                                       "Temp_F_MPI_ESM1_2",
                                                       "Temp_F_MPI_ESM1_2 downscale",
                                                       shift_x = -185,
                                                       shift_y = 16,
                                                       flip_it = F)

Temp_F_MPI_ESM1_2_down

# MODEL UKESM1, SPRING, TEMPERATURE
Temp_S_UKESM1_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_UKESM1-0-LL_ssp585_r2i1p1f2_gn_205001-210012.nc",
                                                       "EPSG:4269",
                                                   lon_call = "longitude",
                                                   lat_call = "latitude",
                                                       "2075",
                                                       "2094",
                                                       c("03", "04", "05"),
                                                       Strata_2,
                                                       coords,
                                                       name,
                                                       "Temp_S_UKESM1",
                                                       "Temp_S_UKESM1 downscale",
                                                       shift_x = -107,
                                                       shift_y = -4,
                                                       flip_it = T)

Temp_S_UKESM1_down

# MODEL UKESM1, SPRING, TEMPERATURE
Temp_F_UKESM1_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_tob_output/tob_Omon_UKESM1-0-LL_ssp585_r2i1p1f2_gn_205001-210012.nc",
                                                   "EPSG:4269",
                                                   lon_call = "longitude",
                                                   lat_call = "latitude",
                                                   "2075",
                                                   "2094",
                                                   c("09", "10", "11"),
                                                   Strata_2,
                                                   coords,
                                                   Temp_F_UKESM1,
                                                   "Temp_F_UKESM1",
                                                   "Temp_F_UKESM1 downscale",
                                                   shift_x = -107,
                                                   shift_y = -4,
                                                   flip_it = T)

Temp_F_UKESM1_down

# STACK RASTERS FROM DOWNSCALED MODELS 
CM6_Spring_Ensemble  <- c(Temp_S_CanESM5_down[[2]], Temp_S_CNRM_CM6_1_down[[2]], Temp_S_MPI_ESM1_2_down[[2]], Temp_S_UKESM1_down[[2]])
CM6_Spring_Ensemble_mean <- terra::app(CM6_Spring_Ensemble, fun = mean, na.rm = T)

ggplot() +
    geom_raster(data = CM6_Spring_Ensemble_mean, aes(x = x, y = y, fill = mean)) +
    scale_fill_distiller(palette = "YlOrRd",  limits = c(5, 30), na.value = "grey", direction = 1) +
    geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
    coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
    labs(title = "Spring Temperature",
         y = "Latitude",
      x = "Longitude",
      fill = "") +
    theme_classic(base_size = 17) -> Spng_Temp_Ensemble

CM6_Fall_Ensemble <- c(Temp_F_CanESM5_down[[2]], Temp_F_CNRM_CM6_1_down[[2]], Temp_F_MPI_ESM1_2_down[[2]], Temp_F_UKESM1_down[[2]])
CM6_Fall_Ensemble_mean <- terra::app(CM6_Fall_Ensemble, fun = mean, na.rm = T)

ggplot() +
    geom_raster(data = CM6_Fall_Ensemble_mean, aes(x = x, y = y, fill = mean)) +
    scale_fill_distiller(palette = "YlOrRd",  limits = c(5, 30), na.value = "grey", direction = 1) +
    geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
    coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
    labs(title = "Fall Temperature",
      y = "Latitude",
      x = "Longitude",
      fill = "") +
    theme_classic(base_size = 17) -> Fall_Temp_Ensemble


########################### SALINITY PROJECTIONS / DOWNSCALING #################
# MODEL CNRM-CM6-1, SPRING, SALINITY 
Salt_S_CNRM_CM6_1_down <- rasterize_2d_downscale_CMIP6(ncfile = "~/Documents/R Projects/ClimateProjections/CMIP6_so_output/sob_Omon_CNRM-CM6-1_ssp585_r1i1p1f2_gn_201501-210012.nc",
                                                       project_crs = "EPSG:4269",
                                                       year1 = "2075",
                                                       year2 = "2094",
                                                       lon_call = "lon",
                                                       lat_call = "lat",
                                                       months = c("03", "04", "05"),
                                                       mask_poly = Strata_2,
                                                       coords = coords,
                                                       name = name,
                                                       t1 = "Salt_S_CNRM-CM6-1",
                                                       t2 = "Salt_S_CNRM-CM6-1 downscale",
                                                       shift_x = -106,
                                                       shift_y = 0,
                                                       flip_it = T)

Salt_S_CNRM_CM6_1_down

# MODEL CNRM-CM6-1, FALL, SALINITY 
Salt_F_CNRM_CM6_1_down <- rasterize_2d_downscale_CMIP6(ncfile = "~/Documents/R Projects/ClimateProjections/CMIP6_so_output/sob_Omon_CNRM-CM6-1_ssp585_r1i1p1f2_gn_201501-210012.nc",
                                                       project_crs = "EPSG:4269",
                                                       year1 = "2075",
                                                       year2 = "2094",
                                                       lon_call = "lon",
                                                       lat_call = "lat",
                                                       months = c("09", "10", "11"),
                                                       mask_poly = Strata_2,
                                                       coords = coords,
                                                       name = name,
                                                       t1 = "Salt_F_CNRM-CM6-1",
                                                       t2 = "Salt_F_CNRM-CM6-1 downscale",
                                                       shift_x = -106,
                                                       shift_y = 0,
                                                       flip_it = T)


Salt_F_CNRM_CM6_1_down

# MODEL CanESM5, SPRING, SALINITY 
Salt_S_CanESM5_down <- rasterize_2d_downscale_CMIP6(ncfile = "~/Documents/R Projects/ClimateProjections/CMIP6_so_output/sob_Omon_CanESM5_ssp585_r10i1p2f1_gn_201501-210012.nc",
                                                    project_crs = "EPSG:4269",
                                                    year1 = "2075",
                                                    year2 = "2094",
                                                    months = c("03", "04", "05"),
                                                    mask_poly = Strata_2,
                                                    coords = coords,
                                                    name = name,
                                                    t1 = "Salt_S_CanESM5",
                                                    t2 = "Salt_S_CanESM5 downscale",
                                                    shift_x = -286,
                                                    shift_y = 0,
                                                    flip_it = T)

Salt_S_CanESM5_down

# MODEL CanESM5, FALL, SALINITY 
Salt_F_CanESM5_down <- rasterize_2d_downscale_CMIP6(ncfile = "~/Documents/R Projects/ClimateProjections/CMIP6_so_output/sob_Omon_CanESM5_ssp585_r10i1p2f1_gn_201501-210012.nc",
                                                    project_crs = "EPSG:4269",
                                                    year1 = "2075",
                                                    year2 = "2094",
                                                    months = c("09", "10", "11"),
                                                    mask_poly = Strata_2,
                                                    coords = coords,
                                                    name = name,
                                                    t1 = "Salt_F_CanESM5",
                                                    t2 = "Salt_F_CanESM5 downscale",
                                                    shift_x = -286,
                                                    shift_y = 0,
                                                    flip_it = T)

Salt_F_CanESM5_down

# MODEL MPI_ESM1_2, SPRING, SALINITY
Salt_S_MPI_ESM1_2_down <- extract_bottom_salinity(salinity_nc = "~/Documents/R Projects/ClimateProjections/CMIP6_so_output/so_Omon_MPI-ESM1-2-LR_ssp585_r28i1p1f1_gn_207501-209412.nc",
                                                  depth_nc = "~/Documents/R Projects/ClimateProjections/CMIP6_deptho_output/deptho_Ofx_MPI-ESM1-2-LR_ssp585_r2i1p1f1_gn.nc",
                                                       "EPSG:4269",
                                                       lon_call = "longitude",
                                                       lat_call = "latitude",
                                                       "2075",
                                                       "2094",
                                                       c("03", "04", "05"),
                                                       Strata_2,
                                                       coords,
                                                       Sali_S_MPI_ESM1_2,
                                                       "Salt_S_MPI_ESM1_2",
                                                       "Salt_S_MPI_ESM1_2 downscale",
                                                       shift_x = -185,
                                                       shift_y = 16,
                                                       flip_it = F)

Salt_S_MPI_ESM1_2_down

Salt_S_MPI_ESM1_2_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_so_output/sob_Omon_MPI-ESM1-2-LR_ssp585_r28i1p1f1_gn_207501-209412.nc",
                                                       "EPSG:4269",
                                                       lon_call = "longitude",
                                                       lat_call = "latitude",
                                                       "2075",
                                                       "2094",
                                                       c("03", "04", "05"),
                                                       Strata_2,
                                                       coords,
                                                       Salt_S_MPI_ESM1_2,
                                                       "Salt_S_MPI_ESM1_2",
                                                       "Salt_S_MPI_ESM1_2 downscale",
                                                       shift_x = -185,
                                                       shift_y = 16,
                                                       flip_it = F)

Salt_S_MPI_ESM1_2_down

# MODEL MPI_ESM1_2, FALL, SALINITY
Salt_F_MPI_ESM1_2_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_so_output/sob_Omon_MPI-ESM1-2-LR_ssp585_r28i1p1f1_gn_207501-209412.nc",
                                                   "EPSG:4269",
                                                   lon_call = "longitude",
                                                   lat_call = "latitude",
                                                   "2075",
                                                   "2094",
                                                   c("09", "10", "11"),
                                                   Strata_2,
                                                   coords,
                                                   Sali_F_MPI_ESM1_2,
                                                   "Salt_F_MPI_ESM1_2",
                                                   "Salt_F_MPI_ESM1_2 downscale",
                                                   shift_x = -185,
                                                   shift_y = 16,
                                                   flip_it = F)
  
Salt_F_MPI_ESM1_2_down

# MODEL UKESM1, SPRING, TEMPERATURE
Salt_S_UKESM1_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_so_output/sob_Omon_UKESM1-0-LL_ssp585_r2i1p1f2_gn_205001-210012.nc",
                                                   "EPSG:4269",
                                                   lon_call = "longitude",
                                                   lat_call = "latitude",
                                                   "2075",
                                                   "2094",
                                                   c("03", "04", "05"),
                                                   Strata_2,
                                                   coords,
                                                   Salt_S_UKESM1,
                                                   "Salt_S_UKESM1",
                                                   "Salt_S_UKESM1 downscale",
                                                   shift_x = -107,
                                                   shift_y = -4,
                                                   flip_it = T)

Salt_S_UKESM1_down

# MODEL UKESM1, FALL, TEMPERATURE
Salt_F_UKESM1_down <- rasterize_2d_downscale_CMIP6("~/Documents/R Projects/ClimateProjections/CMIP6_so_output/sob_Omon_UKESM1-0-LL_ssp585_r2i1p1f2_gn_205001-210012.nc",
                                                   "EPSG:4269",
                                                   lon_call = "longitude",
                                                   lat_call = "latitude",
                                                   "2075",
                                                   "2094",
                                                   c("09", "10", "11"),
                                                   Strata_2,
                                                   coords,
                                                   Salt_F_UKESM1,
                                                   "Salt_F_UKESM1",
                                                   "Salt_F_UKESM1 downscale",
                                                   shift_x = -107,
                                                   shift_y = -4,
                                                   flip_it = T)



Salt_F_UKESM1_down


# Ensemble Salinity 
CM6_Spring_Salt_Ensemble <- c(Salt_S_CanESM5_down[[2]], Salt_S_CNRM_CM6_1_down[[2]], Salt_S_MPI_ESM1_2_down[[2]], Salt_S_UKESM1_down[[2]])
CM6_Spring_Salt_Ensemble_mean <- terra::app(CM6_Spring_Salt_Ensemble, fun = mean, na.rm = T)

ggplot() +
    geom_raster(data = CM6_Spring_Salt_Ensemble_mean, aes(x = x, y = y, fill = mean)) +
    scale_fill_distiller(palette = "Blues",  limits = c(31, 37), na.value = "grey", direction = 1) +
    geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
    coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
    labs(title = "Spring Salinity ",
      y = "Latitude",
      x = "Longitude",
      fill = "") +
    theme_classic(base_size = 17) -> Spng_Salt_Ensemble

CM6_Fall_Salt_Ensemble <- c(Salt_F_CanESM5_down[[2]], Salt_F_CNRM_CM6_1_down[[2]], Salt_F_MPI_ESM1_2_down[[2]], Salt_F_UKESM1_down[[2]])
CM6_Fall_Salt_Ensemble_mean <- terra::app(CM6_Fall_Salt_Ensemble, fun = mean, na.rm = T)

ggplot() +
    geom_raster(data = CM6_Fall_Salt_Ensemble_mean, aes(x = x, y = y, fill = mean)) +
    scale_fill_distiller(palette = "Blues",  limits = c(31, 37), na.value = "grey", direction = 1) +
    geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
    coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
    labs(title = "Fall Salinity",
      y = "Latitude",
      x = "Longitude",
      fill = "") +
    theme_classic(base_size = 17) -> Fall_Salt_Ensemble


(Spng_Temp_Ensemble + Fall_Temp_Ensemble) / (Spng_Salt_Ensemble + Fall_Salt_Ensemble) + plot_layout(guides = "collect") +
  plot_annotation(title = "CMIP 6 Ensemble Mean",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> ensemble_plots

ggsave(filename = "figure_CMIP6_output.png",
       plot = ensemble_plots,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 14,
       width = 20)

################################### DEPTH Raster ###############################
d3 = getNOAA.bathy(lon1 = -78, lon2 = -63, lat1 = 34, lat2 = 45, 
                   resolution = 6) # 1 degree x 1 se
df3 = tibble(fortify.bathy(d3))

depth_raster <- rast(df3, type = "xyz")
crs(depth_raster) <- "EPSG:4269"
depth_raster <- mask(depth_raster, Strata_2)
(depth_raster <- depth_raster*-1)
depth_raster2 <- resample(depth_raster, CM6_Fall_Ensemble_mean)
plot(depth_raster2)



################ Generate HSI information from projections
sal_proj_F_OSNE       <- mask(CM6_Fall_Salt_Ensemble_mean, OSNE)
sal_proj_F_IGOM       <- mask(CM6_Fall_Salt_Ensemble_mean, IGOM)
sal_proj_F_OGOM       <- mask(CM6_Fall_Salt_Ensemble_mean, OGOM)
sal_proj_F_CAN        <- mask(CM6_Fall_Salt_Ensemble_mean, CAN)

temp_proj_F_OSNE       <- mask(CM6_Fall_Ensemble_mean, OSNE)
temp_proj_F_IGOM       <- mask(CM6_Fall_Ensemble_mean, IGOM)
temp_proj_F_OGOM       <- mask(CM6_Fall_Ensemble_mean, OGOM)
temp_proj_F_CAN        <- mask(CM6_Fall_Ensemble_mean, CAN)

sal_proj_S_OSNE       <- mask(CM6_Spring_Salt_Ensemble_mean, OSNE)
sal_proj_S_IGOM       <- mask(CM6_Spring_Salt_Ensemble_mean, IGOM)
sal_proj_S_OGOM       <- mask(CM6_Spring_Salt_Ensemble_mean, OGOM)
sal_proj_S_CAN        <- mask(CM6_Spring_Salt_Ensemble_mean, CAN)

temp_proj_S_OSNE       <- mask(CM6_Spring_Ensemble_mean, OSNE)
temp_proj_S_IGOM       <- mask(CM6_Spring_Ensemble_mean, IGOM)
temp_proj_S_OGOM       <- mask(CM6_Spring_Ensemble_mean, OGOM)
temp_proj_S_CAN        <- mask(CM6_Spring_Ensemble_mean, CAN)

depth_proj_OSNE       <- mask(depth_raster2, OSNE)
depth_proj_IGOM       <- mask(depth_raster2, IGOM)
depth_proj_OGOM       <- mask(depth_raster2, OGOM)
depth_proj_CAN        <- mask(depth_raster2, CAN)



################################ HSIs #########################################
# RUN FIT_HSIs.R to generate
SJ_CAN_fHSI  <- rast_project_HSI(SJ_HSI_CAN, temp_proj_S_CAN, sal_proj_S_CAN, depth_proj_CAN) 
SJ_OGOM_fHSI <- rast_project_HSI(SJ_HSI_OGOM, temp_proj_S_OGOM, sal_proj_S_OGOM, depth_proj_OGOM) 
SJ_IGOM_fHSI <- rast_project_HSI(SJ_HSI_IGOM, temp_proj_S_IGOM, sal_proj_S_IGOM, depth_proj_IGOM) 
SJ_OSNE_fHSI <- rast_project_HSI(SJ_HSI_OSNE, temp_proj_S_OSNE, sal_proj_S_OSNE, depth_proj_OSNE) 
SJ_fHSI <- terra::mosaic(SJ_CAN_fHSI, SJ_IGOM_fHSI, SJ_OGOM_fHSI, SJ_OSNE_fHSI)
plot(SJ_fHSI)

FJ_CAN_fHSI  <- rast_project_HSI(FJ_HSI_CAN, temp_proj_F_CAN, sal_proj_F_CAN, depth_proj_CAN) 
FJ_OGOM_fHSI <- rast_project_HSI(FJ_HSI_OGOM, temp_proj_F_OGOM, sal_proj_F_OGOM, depth_proj_OGOM) 
FJ_IGOM_fHSI <- rast_project_HSI(FJ_HSI_IGOM, temp_proj_F_IGOM, sal_proj_F_IGOM, depth_proj_IGOM) 
FJ_OSNE_fHSI <- rast_project_HSI(FJ_HSI_OSNE, temp_proj_F_OSNE, sal_proj_F_OSNE, depth_proj_OSNE) 
FJ_fHSI <- terra::mosaic(FJ_CAN_fHSI, FJ_IGOM_fHSI, FJ_OGOM_fHSI, FJ_OSNE_fHSI)
plot(FJ_fHSI)


SR_CAN_fHSI  <- rast_project_HSI(SR_HSI_CAN, temp_proj_S_CAN, sal_proj_S_CAN, depth_proj_CAN) 
SR_OGOM_fHSI <- rast_project_HSI(SR_HSI_OGOM, temp_proj_S_OGOM, sal_proj_S_OGOM, depth_proj_OGOM) 
SR_IGOM_fHSI <- rast_project_HSI(SR_HSI_IGOM, temp_proj_S_IGOM, sal_proj_S_IGOM, depth_proj_IGOM) 
SR_OSNE_fHSI <- rast_project_HSI(SR_HSI_OSNE, temp_proj_S_OSNE, sal_proj_S_OSNE, depth_proj_OSNE) 
SR_fHSI <- terra::mosaic(SR_CAN_fHSI, SR_IGOM_fHSI, SR_OGOM_fHSI, SR_OSNE_fHSI)
plot(SR_fHSI)

FR_CAN_fHSI  <- rast_project_HSI(FR_HSI_CAN, temp_proj_F_CAN, sal_proj_F_CAN, depth_proj_CAN) 
FR_OGOM_fHSI <- rast_project_HSI(FR_HSI_OGOM, temp_proj_F_OGOM, sal_proj_F_OGOM, depth_proj_OGOM) 
FR_IGOM_fHSI <- rast_project_HSI(FR_HSI_IGOM, temp_proj_F_IGOM, sal_proj_F_IGOM, depth_proj_IGOM) 
FR_OSNE_fHSI <- rast_project_HSI(FR_HSI_OSNE, temp_proj_F_OSNE, sal_proj_F_OSNE, depth_proj_OSNE) 
FR_fHSI <- terra::mosaic(FR_CAN_fHSI, FR_IGOM_fHSI, FR_OGOM_fHSI, FR_OSNE_fHSI)
plot(FR_fHSI)


#################### Jonah crab projection HSI mas ############################
Stocks2 <- st_read("~/Documents/R Projects/Cancer_HSI2/NOAA_Edit_Shape/editShape.shp")

ggplot() +
  geom_raster(data = SJ_fHSI, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring, 2075-2094",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sj_Spring_HSImap



ggplot() +
    geom_raster(data = FJ_fHSI, aes(x = x, y = y, fill = mean)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
    geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
                 color = "gray65", linewidth = 0.1) +
    geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
    coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
    labs(title = "Fall, 2075-2094",
         y = "Latitude",
         x = "Longitude",
         fill = "SI") +
    theme_classic(base_size = 17) -> sj_Fall_HSImap


((sj_Spring_HSImap) + (sj_Fall_HSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Jonah Crab Habitat Suitability",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> jonahHSI_proj

ggsave(filename = "JonahHSImaps_proj.png",
       plot = jonahHSI_proj,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 8,
       width = 16)

######################### ROCK crab projection HSI maps #######################
ggplot() +
  geom_raster(data = SR_fHSI, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring, 2075-2094",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> sr_HSImap



ggplot() +
    geom_raster(data = FR_fHSI, aes(x = x, y = y, fill = mean)) +
    scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
    geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
    geom_contour(aes(x = x, y = y, z = z),
                 data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
                 color = "gray65", linewidth = 0.1) +
    geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
    geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
    coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
    labs(title = "Fall, 2075-2094",
         y = "Latitude",
         x = "Longitude",
         fill = "SI") +
    theme_classic(base_size = 17) -> fr_HSImap



((sr_HSImap) + (fr_HSImap)) + plot_layout(guides = "collect") +
  plot_annotation(title = "Jonah Crab Habitat Suitability",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> rockHSI_proj

ggsave(filename = "RockHSImaps_proj.png",
       plot = rockHSI_proj,
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 8,
       width = 16)

##### HSIs seem too good to be true and should be nvestigated with individual layers ########

# Use line 1317 down in FIT_HSIs.R 
## Spring Jonah 
# CAN HSI
generate_suitability_indces <- function(HSI, temp_rast, salt_rast, depth_rast){
  HSI |> 
    group_by(row_number()) |>
    mutate(mintemp = get_minpoint(fisher_TEMP), maxtemp = get_maxpoint(fisher_TEMP),
           minsalt = get_minpoint(fisher_SALI), maxsalt = get_maxpoint(fisher_SALI),
           mindept = get_minpoint(fisher_DEPT), maxdept = get_maxpoint(fisher_DEPT)) -> HSI
  
  mat1 <- as.matrix(HSI[, c("mintemp", "maxtemp", "TEMP_SI")])
  df_temp1 = data.frame(mintemp = max(HSI$maxtemp), maxtemp = minmax(temp_rast)[2], TEMP_SI = HSI$TEMP_SI[20] / 2)
  
  if (df_temp1$mintemp > df_temp1$maxtemp){
    df_temp1$maxtemp = df_temp1$mintemp
  }
  
  mat1 <- bind_rows(as.data.frame(mat1), df_temp1)
  TSI_rast <- terra::classify(temp_rast, rcl = mat1, include.lowest = T)
  
  mat2 <- as.matrix(HSI[, c("minsalt", "maxsalt", "SALI_SI")])
  df_temp2 = data.frame(minsalt = max(HSI$minsalt), maxsalt = minmax(salt_rast)[2] + 3.5, SALI_SI = HSI$SALI_SI[20] / 2)
  
  mat2 <- bind_rows(as.data.frame(mat2), df_temp2)
  SSI_rast <- terra::classify(salt_rast, rcl = mat2, include.lowest = T)
  
  mat3 <- as.matrix(HSI[, c("mindept", "maxdept", "DEPT_SI")])
  DSI_rast <- app(depth_rast, fun = depth_lookup_fun, data = mat3)
  DSI_rast <- resample(DSI_rast, TSI_rast)
  
  return(list(
    TSI_rast,
    SSI_rast,
    DSI_rast
  ))
}

#####
Spring_Jonah_CAN_Proj <- generate_suitability_indces(SJ_HSI_CAN, temp_proj_S_CAN, sal_proj_S_CAN, depth_proj_F_CAN)
Spring_Jonah_CAN_Proj <- rast(Spring_Jonah_CAN_Proj)
Spring_Jonah_CAN_HSI_Proj <- terra::app(Spring_Jonah_CAN_Proj, fun = mean, na.rm = T)

Spring_Jonah_OGOM_Proj <- generate_suitability_indces(SJ_HSI_OGOM, temp_proj_S_OGOM, sal_proj_S_OGOM, depth_proj_OGOM)
Spring_Jonah_OGOM_Proj <- rast(Spring_Jonah_OGOM_Proj)
Spring_Jonah_OGOM_HSI_Proj <- terra::app(Spring_Jonah_OGOM_Proj, fun = mean, na.rm = T)

Spring_Jonah_IGOM_Proj <- generate_suitability_indces(SJ_HSI_IGOM, temp_proj_S_IGOM, sal_proj_S_IGOM, depth_proj_IGOM)
Spring_Jonah_IGOM_Proj <- rast(Spring_Jonah_IGOM_Proj)
Spring_Jonah_IGOM_HSI_Proj <- terra::app(Spring_Jonah_IGOM_Proj, fun = mean, na.rm = T)

Spring_Jonah_OSNE_Proj <- generate_suitability_indces(SJ_HSI_OSNE, temp_proj_S_OSNE, sal_proj_S_OSNE, depth_proj_OSNE)
Spring_Jonah_OSNE_Proj <- rast(Spring_Jonah_OSNE_Proj)
Spring_Jonah_OSNE_HSI_Proj <- terra::app(Spring_Jonah_OSNE_Proj, fun = mean, na.rm = T)

SJ_HSI_Proj <- mosaic(Spring_Jonah_CAN_HSI_Proj, Spring_Jonah_OGOM_HSI_Proj, Spring_Jonah_IGOM_HSI_Proj, Spring_Jonah_OSNE_HSI_Proj)

Fall_Jonah_CAN_Proj <- generate_suitability_indces(FJ_HSI_CAN, temp_proj_F_CAN, sal_proj_F_CAN, depth_proj_F_CAN)
Fall_Jonah_CAN_Proj <- rast(Fall_Jonah_CAN_Proj)
Fall_Jonah_CAN_HSI_Proj <- terra::app(Fall_Jonah_CAN_Proj, fun = mean, na.rm = T)

Fall_Jonah_OGOM_Proj <- generate_suitability_indces(FJ_HSI_OGOM, temp_proj_F_OGOM, sal_proj_F_OGOM, depth_proj_OGOM)
Fall_Jonah_OGOM_Proj <- rast(Fall_Jonah_OGOM_Proj)
Fall_Jonah_OGOM_HSI_Proj <- terra::app(Fall_Jonah_OGOM_Proj, fun = mean, na.rm = T)

Fall_Jonah_IGOM_Proj <- generate_suitability_indces(FJ_HSI_IGOM, temp_proj_F_IGOM, sal_proj_F_IGOM, depth_proj_IGOM)
Fall_Jonah_IGOM_Proj <- rast(Fall_Jonah_IGOM_Proj)
Fall_Jonah_IGOM_HSI_Proj <- terra::app(Fall_Jonah_IGOM_Proj, fun = mean, na.rm = T)

Fall_Jonah_OSNE_Proj <- generate_suitability_indces(FJ_HSI_OSNE, temp_proj_F_OSNE, sal_proj_F_OSNE, depth_proj_OSNE)
Fall_Jonah_OSNE_Proj <- rast(Fall_Jonah_OSNE_Proj)
Fall_Jonah_OSNE_HSI_Proj <- terra::app(Fall_Jonah_OSNE_Proj, fun = mean, na.rm = T)

FJ_HSI_Proj <- mosaic(Fall_Jonah_CAN_HSI_Proj, Fall_Jonah_OGOM_HSI_Proj, Fall_Jonah_IGOM_HSI_Proj, Fall_Jonah_OSNE_HSI_Proj)

###### Rock crab 
Spring_Rock_CAN_Proj <- generate_suitability_indces(SR_HSI_CAN, temp_proj_S_CAN, sal_proj_S_CAN, depth_proj_F_CAN)
Spring_Rock_CAN_Proj <- rast(Spring_Rock_CAN_Proj)
Spring_Rock_CAN_HSI_Proj <- terra::app(Spring_Rock_CAN_Proj, fun = mean, na.rm = T)

Spring_Rock_OGOM_Proj <- generate_suitability_indces(SR_HSI_OGOM, temp_proj_S_OGOM, sal_proj_S_OGOM, depth_proj_OGOM)
Spring_Rock_OGOM_Proj <- rast(Spring_Rock_OGOM_Proj)
Spring_Rock_OGOM_HSI_Proj <- terra::app(Spring_Rock_OGOM_Proj, fun = mean, na.rm = T)

Spring_Rock_IGOM_Proj <- generate_suitability_indces(SR_HSI_IGOM, temp_proj_S_IGOM, sal_proj_S_IGOM, depth_proj_IGOM)
Spring_Rock_IGOM_Proj <- rast(Spring_Rock_IGOM_Proj)
Spring_Rock_IGOM_HSI_Proj <- terra::app(Spring_Rock_IGOM_Proj, fun = mean, na.rm = T)

Spring_Rock_OSNE_Proj <- generate_suitability_indces(SR_HSI_OSNE, temp_proj_S_OSNE, sal_proj_S_OSNE, depth_proj_OSNE)
Spring_Rock_OSNE_Proj <- rast(Spring_Rock_OSNE_Proj)
Spring_Rock_OSNE_HSI_Proj <- terra::app(Spring_Rock_OSNE_Proj, fun = mean, na.rm = T)

SR_HSI_Proj <- mosaic(Spring_Rock_CAN_HSI_Proj, Spring_Rock_OGOM_HSI_Proj, Spring_Rock_IGOM_HSI_Proj, Spring_Rock_OSNE_HSI_Proj)

Fall_Rock_CAN_Proj <- generate_suitability_indces(FR_HSI_CAN, temp_proj_F_CAN, sal_proj_F_CAN, depth_proj_F_CAN)
Fall_Rock_CAN_Proj <- rast(Fall_Rock_CAN_Proj)
Fall_Rock_CAN_HSI_Proj <- terra::app(Fall_Rock_CAN_Proj, fun = mean, na.rm = T)

Fall_Rock_OGOM_Proj <- generate_suitability_indces(FR_HSI_OGOM, temp_proj_F_OGOM, sal_proj_F_OGOM, depth_proj_OGOM)
Fall_Rock_OGOM_Proj <- rast(Fall_Rock_OGOM_Proj)
Fall_Rock_OGOM_HSI_Proj <- terra::app(Fall_Rock_OGOM_Proj, fun = mean, na.rm = T)

Fall_Rock_IGOM_Proj <- generate_suitability_indces(FR_HSI_IGOM, temp_proj_F_IGOM, sal_proj_F_IGOM, depth_proj_IGOM)
Fall_Rock_IGOM_Proj <- rast(Fall_Rock_IGOM_Proj)
Fall_Rock_IGOM_HSI_Proj <- terra::app(Fall_Rock_IGOM_Proj, fun = mean, na.rm = T)

Fall_Rock_OSNE_Proj <- generate_suitability_indces(FR_HSI_OSNE, temp_proj_F_OSNE, sal_proj_F_OSNE, depth_proj_OSNE)
Fall_Rock_OSNE_Proj <- rast(Fall_Rock_OSNE_Proj)
Fall_Rock_OSNE_HSI_Proj <- terra::app(Fall_Rock_OSNE_Proj, fun = mean, na.rm = T)

FR_HSI_Proj <- mosaic(Fall_Rock_CAN_HSI_Proj, Fall_Rock_OGOM_HSI_Proj, Fall_Rock_IGOM_HSI_Proj, Fall_Rock_OSNE_HSI_Proj)

################# TSI, SSI, DSI Maps ###########################################
ggplot() +
  geom_raster(data = SJ_HSI_Proj, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> SJ_Proj

ggplot() +
  geom_raster(data = FJ_HSI_Proj, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Fall",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> FJ_Proj

SJ_Proj + FJ_Proj + plot_layout(guides = "collect") +
  plot_annotation(title = "Jonah Crab HSI, RCP 8.5 2075-2094",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> jonah_crab_proj

ggsave(filename = "Jonah_HSI_projectionv2.png",
       plot = jonah_crab_proj, 
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 8,
       width = 16)

ggplot() +
  geom_raster(data = SR_HSI_Proj, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> SR_Proj

ggplot() +
  geom_raster(data = FR_HSI_Proj, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "Spectral",  limits = c(0, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Fall",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> FR_Proj

SR_Proj + FR_Proj + plot_layout(guides = "collect") +
  plot_annotation(title = "Atl. Rock Crab HSI, RCP 8.5 2075-2094",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> rock_crab_proj

ggsave(filename = "Rock_HSI_projectionv2.png",
       plot = rock_crab_proj, 
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 8,
       width = 16)

########### Compare to hindcast 
# env rasters
spngTemp_rast_f <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringTemp.tif")
fallTemp_rast_f <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallTemp.tif")
spngSali_rast_f <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullSpringSali.tif")
fallSali_rast_f <- rast("~/Documents/R Projects/Cancer_HSI2/Data_Raster_Stack/fullFallSali.tif")

spngTemp_rast_f <- terra::mask(spngTemp_rast_f, Strata_2)
fallTemp_rast_f <- terra::mask(fallTemp_rast_f, Strata_2)
spngSali_rast_f <- terra::mask(spngSali_rast_f, Strata_2)
fallSali_rast_f <- terra::mask(fallSali_rast_f, Strata_2)

spngTemp_rast_f <- terra::app(spngTemp_rast_f, fun = mean, na.rm = T)
fallTemp_rast_f <- terra::app(fallTemp_rast_f, fun = mean, na.rm = T)
spngSali_rast_f <- terra::app(spngSali_rast_f, fun = mean, na.rm = T)
fallSali_rast_f <- terra::app(fallSali_rast_f, fun = mean, na.rm = T)

ext(spngTemp_rast_f) <- ext(depth_raster)
ext(fallTemp_rast_f) <- ext(depth_raster)
ext(spngSali_rast_f) <- ext(depth_raster)
ext(fallSali_rast_f) <- ext(depth_raster)

depth_raster <- resample(depth_raster, spngTemp_rast_f)


# SPRING ENV RASTERS
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

depth_rast_ISNE <- mask(depth_raster, ISNE)
depth_rast_OSNE <- mask(depth_raster, OSNE)
depth_rast_IGOM <- mask(depth_raster, IGOM)
depth_rast_OGOM <- mask(depth_raster, OGOM)
depth_rast_CAN <- mask(depth_raster, CAN)

######### 
#####
Spring_Jonah_CAN <- generate_suitability_indces(SJ_HSI_CAN, springTemp_rast_CAN, springSalt_rast_CAN, depth_rast_CAN)
Spring_Jonah_CAN <- rast(Spring_Jonah_CAN)
Spring_Jonah_CAN_HSI <- terra::app(Spring_Jonah_CAN, fun = mean, na.rm = T)

Spring_Jonah_OGOM <- generate_suitability_indces(SJ_HSI_OGOM, springTemp_rast_OGOM, springSalt_rast_OGOM, depth_rast_OGOM)
Spring_Jonah_OGOM <- rast(Spring_Jonah_OGOM)
Spring_Jonah_OGOM_HSI <- terra::app(Spring_Jonah_OGOM, fun = mean, na.rm = T)

Spring_Jonah_IGOM <- generate_suitability_indces(SJ_HSI_IGOM, springTemp_rast_IGOM, springSalt_rast_IGOM, depth_rast_IGOM)
Spring_Jonah_IGOM <- rast(Spring_Jonah_IGOM)
Spring_Jonah_IGOM_HSI <- terra::app(Spring_Jonah_IGOM, fun = mean, na.rm = T)

Spring_Jonah_OSNE <- generate_suitability_indces(SJ_HSI_OSNE, springTemp_rast_OSNE, springSalt_rast_OSNE, depth_rast_OSNE)
Spring_Jonah_OSNE <- rast(Spring_Jonah_OSNE)
Spring_Jonah_OSNE_HSI <- terra::app(Spring_Jonah_OSNE, fun = mean, na.rm = T)

SJ_HSI <- mosaic(Spring_Jonah_CAN_HSI, Spring_Jonah_OGOM_HSI, Spring_Jonah_IGOM_HSI, Spring_Jonah_OSNE_HSI)

Fall_Jonah_CAN <- generate_suitability_indces(FJ_HSI_CAN, fallTemp_rast_CAN, fallSalt_rast_CAN, depth_rast_CAN)
Fall_Jonah_CAN <- rast(Fall_Jonah_CAN)
Fall_Jonah_CAN_HSI <- terra::app(Fall_Jonah_CAN, fun = mean, na.rm = T)

Fall_Jonah_OGOM <- generate_suitability_indces(FJ_HSI_OGOM, fallTemp_rast_OGOM, fallSalt_rast_OGOM, depth_rast_OGOM)
Fall_Jonah_OGOM <- rast(Fall_Jonah_OGOM)
Fall_Jonah_OGOM_HSI <- terra::app(Fall_Jonah_OGOM, fun = mean, na.rm = T)

Fall_Jonah_IGOM <- generate_suitability_indces(FJ_HSI_IGOM, fallTemp_rast_IGOM, fallSalt_rast_IGOM, depth_rast_IGOM)
Fall_Jonah_IGOM <- rast(Fall_Jonah_IGOM)
Fall_Jonah_IGOM_HSI <- terra::app(Fall_Jonah_IGOM, fun = mean, na.rm = T)

Fall_Jonah_OSNE <- generate_suitability_indces(FJ_HSI_OSNE, fallTemp_rast_OSNE, fallSalt_rast_OSNE, depth_rast_OSNE)
Fall_Jonah_OSNE <- rast(Fall_Jonah_OSNE)
Fall_Jonah_OSNE_HSI <- terra::app(Fall_Jonah_OSNE, fun = mean, na.rm = T)

FJ_HSI <- mosaic(Fall_Jonah_CAN_HSI, Fall_Jonah_OGOM_HSI, Fall_Jonah_IGOM_HSI, Fall_Jonah_OSNE_HSI)

###### Rock crab 
Spring_Rock_CAN <- generate_suitability_indces(SR_HSI_CAN, springTemp_rast_CAN, springSalt_rast_CAN, depth_rast_CAN)
Spring_Rock_CAN <- rast(Spring_Rock_CAN)
Spring_Rock_CAN_HSI <- terra::app(Spring_Rock_CAN, fun = mean, na.rm = T)

Spring_Rock_OGOM <- generate_suitability_indces(SR_HSI_OGOM, springTemp_rast_OGOM, springSalt_rast_OGOM, depth_rast_OGOM)
Spring_Rock_OGOM <- rast(Spring_Rock_OGOM)
Spring_Rock_OGOM_HSI <- terra::app(Spring_Rock_OGOM, fun = mean, na.rm = T)

Spring_Rock_IGOM <- generate_suitability_indces(SR_HSI_IGOM, springTemp_rast_IGOM, springSalt_rast_IGOM, depth_rast_IGOM)
Spring_Rock_IGOM <- rast(Spring_Rock_IGOM)
Spring_Rock_IGOM_HSI <- terra::app(Spring_Rock_IGOM, fun = mean, na.rm = T)

Spring_Rock_OSNE <- generate_suitability_indces(SR_HSI_OSNE, springTemp_rast_OSNE, springSalt_rast_OSNE, depth_rast_OSNE)
Spring_Rock_OSNE <- rast(Spring_Rock_OSNE)
Spring_Rock_OSNE_HSI <- terra::app(Spring_Rock_OSNE, fun = mean, na.rm = T)

SR_HSI <- mosaic(Spring_Rock_CAN_HSI, Spring_Rock_OGOM_HSI, Spring_Rock_IGOM_HSI, Spring_Rock_OSNE_HSI)

Fall_Rock_CAN <- generate_suitability_indces(FR_HSI_CAN, fallTemp_rast_CAN, fallSalt_rast_CAN, depth_rast_CAN)
Fall_Rock_CAN <- rast(Fall_Rock_CAN)
Fall_Rock_CAN_HSI <- terra::app(Fall_Rock_CAN, fun = mean, na.rm = T)

Fall_Rock_OGOM <- generate_suitability_indces(FR_HSI_OGOM, fallTemp_rast_OGOM, fallSalt_rast_OGOM, depth_rast_OGOM)
Fall_Rock_OGOM <- rast(Fall_Rock_OGOM)
Fall_Rock_OGOM_HSI <- terra::app(Fall_Rock_OGOM, fun = mean, na.rm = T)

Fall_Rock_IGOM <- generate_suitability_indces(FR_HSI_IGOM, fallTemp_rast_IGOM, fallSalt_rast_IGOM, depth_rast_IGOM)
Fall_Rock_IGOM <- rast(Fall_Rock_IGOM)
Fall_Rock_IGOM_HSI <- terra::app(Fall_Rock_IGOM, fun = mean, na.rm = T)

Fall_Rock_OSNE <- generate_suitability_indces(FR_HSI_OSNE, fallTemp_rast_OSNE, fallSalt_rast_OSNE, depth_rast_OSNE)
Fall_Rock_OSNE <- rast(Fall_Rock_OSNE)
Fall_Rock_OSNE_HSI <- terra::app(Fall_Rock_OSNE, fun = mean, na.rm = T)

FR_HSI <- mosaic(Fall_Rock_CAN_HSI, Fall_Rock_OGOM_HSI, Fall_Rock_IGOM_HSI, Fall_Rock_OSNE_HSI)

################### JONAH CRAB #####################
ext(SJ_HSI) <- ext(SJ_HSI_Proj)
SJ_HSI <- resample(SJ_HSI, SJ_HSI_Proj)

diff_SJ_HSI <- SJ_HSI_Proj - SJ_HSI

ggplot() +
  geom_raster(data = diff_SJ_HSI, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "RdBu",  limits = c(-1, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> SJ_HSI_diff_map

ext(FJ_HSI) <- ext(FJ_HSI_Proj)
FJ_HSI <- resample(FJ_HSI, FJ_HSI_Proj)

diff_FJ_HSI <- FJ_HSI_Proj - FJ_HSI

ggplot() +
  geom_raster(data = diff_FJ_HSI, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "RdBu",  limits = c(-1, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Fall",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> FJ_HSI_diff_map

SJ_HSI_diff_map + FJ_HSI_diff_map + plot_layout(guides = "collect") +
  plot_annotation(title = "Jonah Crab HSI Differences 1992-2024 vs RCP 8.5 2075-2094",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> jonah_crab_diff

ggsave(filename = "Jonah_HSI_difference.png",
       plot = jonah_crab_diff, 
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 8,
       width = 16)

################### ROCK CRAB #####################
ext(SR_HSI) <- ext(SR_HSI_Proj)
SR_HSI <- resample(SR_HSI, SR_HSI_Proj)

diff_SR_HSI <- SR_HSI_Proj - SR_HSI

ggplot() +
  geom_raster(data = diff_SR_HSI, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "RdBu",  limits = c(-1, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Spring",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> SR_HSI_diff_map

ext(FR_HSI) <- ext(FR_HSI_Proj)
FR_HSI <- resample(FR_HSI, FR_HSI_Proj)

diff_FR_HSI <- FR_HSI_Proj - FR_HSI

ggplot() +
  geom_raster(data = diff_FR_HSI, aes(x = x, y = y, fill = mean)) +
  scale_fill_distiller(palette = "RdBu",  limits = c(-1, 1), na.value = "grey") +
  geom_map(data = countries, map = countries, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_map(data = states, map = states, aes(long, lat, map_id = region), color = "grey99", fill = "grey") +
  geom_contour(aes(x = x, y = y, z = z),
               data = bf2, breaks = c(0, -50, -100, -200, -500, -1000),
               color = "gray65", linewidth = 0.1) +
  geom_sf(data = Stocks2, aes(group = Stock, color = Stock), alpha = 0.1, linewidth = 0.5) +
  geom_label(data = Stock_Label, aes(x = lon, y = lat, label = Stock), alpha = 0.1, label.size = 0.) +
  coord_sf(xlim = c(-76, -65.7), ylim= c(36, 44.5)) +
  labs(title = "Fall",
       y = "Latitude",
       x = "Longitude",
       fill = "SI") +
  theme_classic(base_size = 17) -> FR_HSI_diff_map


SR_HSI_diff_map + FR_HSI_diff_map + plot_layout(guides = "collect") +
  plot_annotation(title = "Atl. Rock Crab HSI Differences 1992-2024 vs RCP 8.5 2075-2094",
                  tag_levels = "A",
                  theme = theme(title = element_text(size = 20))) -> rock_crab_diff

ggsave(filename = "Rock_HSI_difference.png",
       plot = rock_crab_diff, 
       path = "~/Documents/R Projects/Cancer_HSI2/Figures",
       dpi = 750,
       height = 8,
       width = 16)

## PERCENT HABITAT ABOVE THRESHOLDS
percent_above_threshold <- function(r, threshold = 0.5) {
  total <- global(!is.na(r), fun = "sum", na.rm = FALSE)
  above <- global(r > threshold, fun = "sum", na.rm = TRUE)
  (percent <- (above / total) * 100)
}

# Jonah crab
# Spring 
# Moderate habitat
percent_above_threshold(mask(SJ_HSI_Proj, CAN)) - percent_above_threshold(mask(SJ_HSI_1224, CAN))
percent_above_threshold(mask(SJ_HSI_Proj, OGOM)) - percent_above_threshold(mask(SJ_HSI_1224, OGOM))
percent_above_threshold(mask(SJ_HSI_Proj, IGOM)) - percent_above_threshold(mask(SJ_HSI_1224, IGOM))
percent_above_threshold(mask(SJ_HSI_Proj, OSNE)) - percent_above_threshold(mask(SJ_HSI_1224, OSNE))

# Optimal habitat
percent_above_threshold(mask(SJ_HSI_Proj, CAN), threshold = 0.8) - percent_above_threshold(mask(SJ_HSI_1224, CAN), threshold = 0.8)
percent_above_threshold(mask(SJ_HSI_Proj, OGOM), threshold = 0.8) - percent_above_threshold(mask(SJ_HSI_1224, OGOM), threshold = 0.8)
percent_above_threshold(mask(SJ_HSI_Proj, IGOM), threshold = 0.8) - percent_above_threshold(mask(SJ_HSI_1224, IGOM), threshold = 0.8)
percent_above_threshold(mask(SJ_HSI_Proj, OSNE), threshold = 0.8) - percent_above_threshold(mask(SJ_HSI_1224, OSNE), threshold = 0.8)

# Fall 
# Moderate habitat
percent_above_threshold(mask(FJ_HSI_Proj, CAN)) - percent_above_threshold(mask(FJ_HSI_1224, CAN))
percent_above_threshold(mask(FJ_HSI_Proj, OGOM)) - percent_above_threshold(mask(FJ_HSI_1224, OGOM))
percent_above_threshold(mask(FJ_HSI_Proj, IGOM)) - percent_above_threshold(mask(FJ_HSI_1224, IGOM))
percent_above_threshold(mask(FJ_HSI_Proj, OSNE)) - percent_above_threshold(mask(FJ_HSI_1224, OSNE))

# Optimal habitat
percent_above_threshold(mask(FJ_HSI_Proj, CAN), threshold = 0.8)  - percent_above_threshold(mask(FJ_HSI_1224, CAN), threshold = 0.8)
percent_above_threshold(mask(FJ_HSI_Proj, OGOM), threshold = 0.8) - percent_above_threshold(mask(FJ_HSI_1224, OGOM), threshold = 0.8)
percent_above_threshold(mask(FJ_HSI_Proj, IGOM), threshold = 0.8) - percent_above_threshold(mask(FJ_HSI_1224, IGOM), threshold = 0.8)
percent_above_threshold(mask(FJ_HSI_Proj, OSNE), threshold = 0.8) - percent_above_threshold(mask(FJ_HSI_1224, OSNE), threshold = 0.8)

# Atl. rock crab
# Spring 
# Moderate habitat
percent_above_threshold(mask(SR_HSI_Proj, CAN)) - percent_above_threshold(mask(SR_HSI_1224, CAN))
percent_above_threshold(mask(SR_HSI_Proj, OGOM)) - percent_above_threshold(mask(SR_HSI_1224, OGOM))
percent_above_threshold(mask(SR_HSI_Proj, IGOM)) - percent_above_threshold(mask(SR_HSI_1224, IGOM))
percent_above_threshold(mask(SR_HSI_Proj, OSNE)) - percent_above_threshold(mask(SR_HSI_1224, OSNE))

# Optimal habitat
percent_above_threshold(mask(SR_HSI_Proj, CAN), threshold = 0.8) - percent_above_threshold(mask(SR_HSI_1224, CAN), threshold = 0.8)
percent_above_threshold(mask(SR_HSI_Proj, OGOM), threshold = 0.8) - percent_above_threshold(mask(SR_HSI_1224, OGOM), threshold = 0.8)
percent_above_threshold(mask(SR_HSI_Proj, IGOM), threshold = 0.8) - percent_above_threshold(mask(SR_HSI_1224, IGOM), threshold = 0.8)
percent_above_threshold(mask(SR_HSI_Proj, OSNE), threshold = 0.8) - percent_above_threshold(mask(SR_HSI_1224, OSNE), threshold = 0.8)

# Fall 
# Moderate habitat
percent_above_threshold(mask(FR_HSI_Proj, CAN)) - percent_above_threshold(mask(FR_HSI_1224, CAN))
percent_above_threshold(mask(FR_HSI_Proj, OGOM)) - percent_above_threshold(mask(FR_HSI_1224, OGOM))
percent_above_threshold(mask(FR_HSI_Proj, IGOM)) - percent_above_threshold(mask(FR_HSI_1224, IGOM))
percent_above_threshold(mask(FR_HSI_Proj, OSNE)) - percent_above_threshold(mask(FR_HSI_1224, OSNE))

# Optimal habitat
percent_above_threshold(mask(FR_HSI_Proj, CAN), threshold = 0.8)  - percent_above_threshold(mask(FR_HSI_1224, CAN), threshold = 0.8)
percent_above_threshold(mask(FR_HSI_Proj, OGOM), threshold = 0.8) - percent_above_threshold(mask(FR_HSI_1224, OGOM), threshold = 0.8)
percent_above_threshold(mask(FR_HSI_Proj, IGOM), threshold = 0.8) - percent_above_threshold(mask(FR_HSI_1224, IGOM), threshold = 0.8)
percent_above_threshold(mask(FR_HSI_Proj, OSNE), threshold = 0.8) - percent_above_threshold(mask(FR_HSI_1224, OSNE), threshold = 0.8)
