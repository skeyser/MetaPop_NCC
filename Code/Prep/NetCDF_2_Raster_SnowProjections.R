## ----------------------------------------------------------
##
## Script name: Future Snow Depth Projections
##
## Script purpose: Formatting and working with future
## snow depth projections from downscaled RCM
##
## Author: Spencer R Keyser
##
## Date Created: 2025-06-04
##
## Email: skeyser@wisc.edu
##
## Github: https://github.com/skeyser
##
## -----------------------------------------------------------
##
## Notes:
##
##
## -----------------------------------------------------------

## Defaults
options(scipen = 6, digits = 4)

## -----------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(data.table)
library(here)
library(sf)
library(terra)
library(raster)
library(ncdf4)
## -----------------------------------------------------------

## Load in the extent for the marten study
cost <- rast(here("Inputs/Cost_Water_1km_2.asc"))
plot(cost)
roi <- ext(cost)
crs.cost <- crs(cost)

## Path to netCDF files
snow.path <- list.files("R:/Users/skeyser/PhD/NASA Project/Marten_Data/SnowDepthProjs/", full.names = T)
snow.path <- list.files("R:/Users/skeyser/PhD/NASA Project/Marten_Data/Expanded_SnowDepth_Future/", full.names = T)
# # Read in first snow depth file to check
# snow.nc <- ncdf4::nc_open(snow.path[7])
# print(snow.nc)
#
# attributes(snow.nc$var)
# attributes(snow.nc$dim)
#
# snow.nc$var
# snow.nc
#
# lat <- ncvar_get(snow.nc, "xlat")
# nlat <- dim(lat)
#
# lon <- ncvar_get(snow.nc, "xlon")
# nlon <- dim(lon)
#
# time <- ncvar_get(snow.nc, "time")
# head(time)
# tunits <- ncatt_get(snow.nc, "time", "units")
#
# ## Create separate lists for each RCM
# pr_list <- list(
#   sp.access <- snow.path[stringr::str_detect(snow.path, "access")],
#   sp.cnrm <- snow.path[stringr::str_detect(snow.path, "cnrm")],
#   sp.gfdl <- snow.path[stringr::str_detect(snow.path, "gfdl")],
#   sp.ipsl <- snow.path[stringr::str_detect(snow.path, "ipsl")],
#   sp.miroc <- snow.path[stringr::str_detect(snow.path, "miroc5")],
#   sp.mri <- snow.path[stringr::str_detect(snow.path, "mri")]
# )

## Projection is set to LAMCON per Notaro
prj <- paste("+proj=lcc",  # Lambert Conformal Conic
             "+lat_0=45",    # Latitude of projection origin
             "+lon_0=-97",   # Longitude of projection origin
             "+lat_1=36",    # First standard parallel
             "+lat_2=52",    # Second standard parallel
             "+x_0=0",       # False easting
             "+y_0=0",       # False northing
             "+datum=WGS84", # Datum
             "+units=m",     # Units in meters
             sep=" ")


## Load all the rasters in
snow.ras <- lapply(snow.path, function(x) {rast(x, drivers = "NETCDF") + 0})

## Set PCS for all layers
snow.ras <- lapply(snow.ras, function(x) {
  crs(x) <- prj
  return(x)
  })

## Find the maximum seasonal snow values for each layer
snow.ras.max <- lapply(snow.ras, function(x) {
  x[[grep("max", names(x))]]})

snow.ras.first <- lapply(snow.ras, function(x) {
  x[[grep("first", names(x))]]})

snow.ras.last <- lapply(snow.ras, function(x) {
  x[[grep("last", names(x))]]})

## Rename the layers
rcm.names <- c("access_late20", "access_late21", "access_mid21",
               "cnrm_late20", "cnrm_late21", "cnrm_mid21",
               "gfdl_late20", "gfdl_late21", "gfdl_mid21",
               "ipsl_late20", "ipsl_late21", "ipsl_mid21",
               "miroc5_late20", "miroc5_late21", "miroc5_mid21",
               "mri_late20", "mri_late21", "mri_mid21")

names(snow.ras.max) <- rcm.names
names(snow.ras.first) <- rcm.names
names(snow.ras.last) <- rcm.names

## Project to the CRS for the project
snow.rp.max <- lapply(snow.ras, function(x) project(x, crs.cost))
snow.rp.first <- lapply(snow.ras.first, function(x) project(x, crs.cost))
snow.rp.last <- lapply(snow.ras.last, function(x) project(x, crs.cost))

## Handle the weird orientation issue from NCDF
snow.c.max <- lapply(snow.rp.max, function(x) flip(x))
snow.c.first <- lapply(snow.rp.first, function(x) flip(x))
snow.c.last <- lapply(snow.rp.last, function(x) flip(x))

## Crop to the ROI with a buffer
snow.c.max <- lapply(snow.c.max, function(x) crop(x, roi + 250000))
snow.c.first <- lapply(snow.c.first, function(x) crop(x, roi + 250000))
snow.c.last <- lapply(snow.c.last, function(x) crop(x, roi + 250000))

## Interpolating into the lakes to deal with edge effects
snow.c.max <- lapply(snow.c.max, function(x){
  focal(x, 7, "mean", na.policy = "only", na.rm = T)
})
snow.c.first <- lapply(snow.c.first, function(x){
  focal(x, 7, "mean", na.policy = "only", na.rm = T)
})
snow.c.last <- lapply(snow.c.last, function(x){
  focal(x, 7, "mean", na.policy = "only", na.rm = T)
})

## Check
plot(snow.c.first$access_late20)
plot(snow.c.last$access_late20)

snow.dur <- snow.c.last - snow.c.first
snow.dur <- mapply(function(x, y) x - y, snow.c.last, snow.c.first, SIMPLIFY = FALSE)
plot(snow.dur[[1]])
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Calculate the mean per model for each time
## bin - Baseline, Mid, and Future snow cover
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Calculate % difference raster between historical and mid/late century
late20 <- snow.c[grep("late20", names(snow.c))]
mid21 <- snow.c[grep("mid21", names(snow.c))]
late21 <- snow.c[grep("late21", names(snow.c))]

late20.dur <- snow.dur[grep("late20", names(snow.dur))]
mid21.dur <- snow.dur[grep("mid21", names(snow.dur))]
late21.dur <- snow.dur[grep("late21", names(snow.dur))]

## Calculate average max for the baseline
baseline <- lapply(late20, function(x) app(x, mean))
mid21_mean <- lapply(mid21, function(x) app(x, mean))
late21_mean <- lapply(late21, function(x) app(x, mean))

baseline_dur <- lapply(late20.dur, function(x) app(x, mean))
mid21_mean_dur <- lapply(mid21.dur, function(x) app(x, mean))
late21_mean_dur <- lapply(late21.dur, function(x) app(x, mean))

## Verify the output
baseline[1]
baseline_dur[1]

## Estimate the % change between baseline to mid AND baseline to late
pchange_fun <- function(r1, r2){

  diff.rast <- r2

  for(i in 1:length(diff.rast)){
    for(j in 1:nlyr(diff.rast[[1]])){

      r1.tmp <- r1[[i]]
      r2.tmp <- r2[[i]][[j]]

      ## Percent difference
      p.change <- ((r2.tmp - r1.tmp) / r1.tmp) * 100
      diff.rast[[i]][[j]] <- p.change

    }
  }
  return(diff.rast)
}

## Execute the function across the means
mid.change.mean <- pchange_fun(r1 = baseline, r2 = mid21_mean)
late.change.mean <- pchange_fun(r1 = baseline, r2 = late21_mean)

mid.change.mean.dur <- pchange_fun(r1 = baseline_dur, r2 = mid21_mean_dur)
late.change.mean.dur <- pchange_fun(r1 = baseline_dur, r2 = late21_mean_dur)


## Visualizations
nr <- length(mid.change.mean)/2
nc <- (length(mid.change.mean)/2)-1

par(mfrow = c(nr, nc))
{
  plot(mid.change.mean$access_mid21)
  plot(mid.change.mean$cnrm_mid21)
  plot(mid.change.mean$gfdl_mid21)
  plot(mid.change.mean$ipsl_mid21)
  plot(mid.change.mean$miroc5_mid21)
  plot(mid.change.mean$mri_mid21)
}
{
  plot(late.change.mean$access_late21)
  plot(late.change.mean$cnrm_late21)
  plot(late.change.mean$gfdl_late21)
  plot(late.change.mean$ipsl_late21)
  plot(late.change.mean$miroc5_late21)
  plot(late.change.mean$mri_late21)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Disaggregate to match the spatial resolution
## for marten data
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Resample to 1 km
template <- rast(ext(mid.change.mean$access_mid21), res=1000) # res in meters
crs(template) <- crs(mid.change.mean$access_mid21)

mid.change1_km <- lapply(mid.change.mean, function(x){ resample(x, template, method = "bilinear") })
late.change1_km <- lapply(late.change.mean, function(x){ resample(x, template, method = "bilinear") })

dev.off()
plot(mid.change1_km$access_mid21)
plot(late.change1_km$access_late21)

## Make the lists a stacked raster
mid.change.ras <- rast(mid.change1_km)
late.change.ras <- rast(late.change1_km)

## Write the rasters to file
writeRaster(mid.change.ras, filename = here("./Data/MidCentury_PctChange.tif"), overwrite = T)
writeRaster(late.change.ras, filename = here("./Data/LateCentury_PctChange.tif"), overwrite = T)

# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ##
# ## Subsection: Interpolating the yearly values
# ##
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# late20df <- as.data.frame(late20$access_late20, xy = T)
# late20df <- late20df |>
#   tidyr::pivot_longer(cols = contains("max"), names_to = "Year", values_to = "Snow") |>
#   mutate(
#     year_num = as.numeric(stringr::str_extract(Year, "\\d+")),
#     Year = year_num + 1979
#   ) |>
#   dplyr::select(-year_num)
#
# mid21df <- as.data.frame(mid21$access_mid21, xy = T)
# mid21df <- mid21df |>
#   tidyr::pivot_longer(cols = contains("max"), names_to = "Year", values_to = "Snow") |>
#   mutate(
#     year_num = as.numeric(stringr::str_extract(Year, "\\d+")),
#     Year = year_num + 2039
#   ) |>
#   dplyr::select(-year_num)
#
# late21df <- as.data.frame(late21$access_late21, xy = T)
# late21df <- late21df |>
#   tidyr::pivot_longer(cols = contains("max"), names_to = "Year", values_to = "Snow") |>
#   mutate(
#     year_num = as.numeric(stringr::str_extract(Year, "\\d+")),
#     Year = year_num + 2079
#   ) |>
#   dplyr::select(-year_num)
#
#
# ## Full ts
# snow.ts <- bind_rows(late20df, mid21df, late21df)
# snow.ts <- snow.ts |>
#   group_by(x, y) |>
#   tidyr::complete(Year = 1980:2099) |>
#   ungroup()
#
# ## Interpolating the data
# str(snow.ts)
#
# library(mgcv)
# library(dplyr)
#
# # First, create spatial factor for random effect
# df <- snow.ts %>%
#   group_by(x, y) |>
#   mutate(
#     site_id = group_indices(),
#     space_id = factor(paste(x, y)),
#     Year_c = scale(Year, scale = FALSE)
#   ) |>
#   ungroup()
#
# # Fit GAM with spatial random effect and AR temporal structure
# m <- bam(Snow ~ s(Year_c, k = 20) + # smooth term for year
#            s(space_id, bs = 're')     # spatial random effect
#          data = df,
#          method = 'REML',
#          AR.start = df$Year == min(df$Year, na.rm = TRUE), # Start of AR process
#          rho = 0.7) # Initial AR1 correlation parameter
#
# # Predict missing values
# df$Snow_pred <- predict(m)
#
# # Fill in NAs with predictions
# df$Snow_filled <- ifelse(is.na(df$Snow), df$Snow_pred, df$Snow)

