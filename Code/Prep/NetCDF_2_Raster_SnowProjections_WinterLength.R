## ----------------------------------------------------------
##
## Script name: Future Snow Season Length Projections for Ben
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

## Path to netCDF files
snow.path <- list.files("R:/Users/skeyser/PhD/NASA Project/Marten_Data/Expanded_SnowDepth_Future/", full.names = T)

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

names(snow.ras.first) <- rcm.names
names(snow.ras.last) <- rcm.names

## Project to the CRS for the project
snow.rp.first <- lapply(snow.ras.first, function(x) project(x, crs.cost))
snow.rp.last <- lapply(snow.ras.last, function(x) project(x, crs.cost))

## Handle the weird orientation issue from NCDF
snow.c.first <- lapply(snow.ras.first, function(x) flip(x))
snow.c.last <- lapply(snow.ras.last, function(x) flip(x))

## Check
plot(snow.c.first$access_late20)
plot(snow.c.last$access_late20)

snow.dur <- mapply(function(x, y) x - y, snow.c.last, snow.c.first, SIMPLIFY = FALSE)
plot(snow.dur[[1]])
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Calculate the mean per model for each time
## bin - Baseline, Mid, and Future snow cover
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Calculate % difference raster between historical and mid/late century
late20.dur <- snow.dur[grep("late20", names(snow.dur))]
mid21.dur <- snow.dur[grep("mid21", names(snow.dur))]
late21.dur <- snow.dur[grep("late21", names(snow.dur))]

## Calculate average max for the baseline
baseline_dur <- lapply(late20.dur, function(x) app(x, mean))
mid21_mean_dur <- lapply(mid21.dur, function(x) app(x, mean))
late21_mean_dur <- lapply(late21.dur, function(x) app(x, mean))

## Verify the output
plot(baseline_dur[1]$access_late20)
plot(mid21_mean_dur[1]$access_mid21)
plot(late21_mean_dur[1]$access_late21)

## Write these files
str(baseline_dur_r)
baseline_dur_r <- rast(baseline_dur)
mid21_mean_dur_r <- rast(mid21_mean_dur)
late21_mean_dur_r <- rast(late21_mean_dur)
writeRaster(baseline_dur_r, here("./Data/BaselineWinterLength.tif"))
writeRaster(mid21_mean_dur_r, here("./Data/Mid21WinterLength.tif"))
writeRaster(late21_mean_dur_r, here("./Data/Late21WinterLength.tif"))


## Estimate the % change between baseline to mid AND baseline to late
pchange_fun <- function(r1, r2){

  diff.rast <- r2
  loss.rast <- r2

  for(i in 1:length(diff.rast)){
    for(j in 1:nlyr(diff.rast[[1]])){

      r1.tmp <- r1[[i]]
      r2.tmp <- r2[[i]][[j]]

      ## Percent difference
      p.change <- ((r2.tmp - r1.tmp) / r1.tmp) * 100
      diff.rast[[i]][[j]] <- p.change

      ## Track complete loss (where r1 had snow but r2 has none)
      # Consider both NA and 0 as no snow
      loss <- ((!is.na(r1.tmp) & r1.tmp > 0) & (is.na(r2.tmp) | r2.tmp == 0))
      loss.rast[[i]][[j]] <- loss

    }
  }
  return(list(percent_change = diff.rast,
              complete_loss = loss.rast))
}

## Execute the function across the means
mid.change.mean.dur <- pchange_fun(r1 = baseline_dur, r2 = mid21_mean_dur)
late.change.mean.dur <- pchange_fun(r1 = baseline_dur, r2 = late21_mean_dur)


## Visualizations
nr <- length(mid.change.mean.dur)/2
nc <- (length(mid.change.mean.dur)/2)-1

par(mfrow = c(nr, nc))
{
  plot(mid.change.mean.dur$access_mid21)
  plot(mid.change.mean.dur$cnrm_mid21)
  plot(mid.change.mean.dur$gfdl_mid21)
  plot(mid.change.mean.dur$ipsl_mid21)
  plot(mid.change.mean.dur$miroc5_mid21)
  plot(mid.change.mean.dur$mri_mid21)
}
{
  plot(late.change.mean.dur$access_late21)
  plot(late.change.mean.dur$cnrm_late21)
  plot(late.change.mean.dur$gfdl_late21)
  plot(late.change.mean.dur$ipsl_late21)
  plot(late.change.mean.dur$miroc5_late21)
  plot(late.change.mean.dur$mri_late21)
}

## Just take the access data
mid.change <- mid.change.mean.dur$complete_loss$access_mid21
late.change <- late.change.mean.dur$complete_loss$access_late21
mid.change <- mid.change.mean.dur$percent_change$access_mid21
late.change <- late.change.mean.dur$percent_change$access_late21

## % change rasters
writeRaster(mid.change, filename = "R:/Users/skeyser/PhD/NASA Project/Marten_Data/MidCentury_WinterLength_PctChange.tif")
writeRaster(late.change, filename = "R:/Users/skeyser/PhD/NASA Project/Marten_Data/LateCentury_WinterLength_PctChange.tif")

mid.late.change <- c(mid.change, late.change)
names(mid.late.change) <- c("PctChange_Mid_WL", "PctChange_Late_WL")

writeRaster(mid.late.change, filename = "R:/Users/skeyser/PhD/NASA Project/Marten_Data/MidLateCentury_WinterLength_PctChange.tif", overwrite = T)


plot(mid.change)

mid.change.df <- as.data.frame(mid.change, xy = T)
colnames(mid.change.df) <- c("x", "y", "pct_change")
late.change.df <- as.data.frame(late.change, xy = T)
colnames(late.change.df) <- c("x", "y", "pct_change")

## Add in a polygon for reference
region <- st_read("G:/eBird/US_Map/Canada_USA_map.shp") |>
  filter(!PROVINCE_S %in% c("HAWAII")) |>
  st_transform(crs = crs(mid.change))

ggplot() +
  geom_raster(data = mid.change.df, aes(x=x, y=y, fill = pct_change)) +
  scale_fill_viridis_c(option = "G", direction = -1, limits = c(-70, 0)) +
  geom_sf(data = region, size = 1, fill = NA, color = "black") +
  theme_void()

ggplot() +
  geom_raster(data = late.change.df, aes(x=x, y=y, fill = pct_change)) +
  scale_fill_viridis_c(option = "G", direction = -1, limits = c(-70, 0)) +
  geom_sf(data = region, size = 1, fill = NA, color = "black") +
  theme_void()


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
