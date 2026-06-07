## ----------------------------------------------------------
##
## Script name: Marten Survival Projections
##
## Script purpose: Use snow depth and entropy to predict
## marten survival rates across the identified patches.
##
## Author: Spencer R Keyser
##
## Date Created: 2025-04-19
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
library(here)
library(dplyr)
library(ggplot2)
library(data.table)
library(sf)
library(terra)
library(mapview)
library(tidyverse)

## -----------------------------------------------------------

######################################################################################################################################################
################################## Read in data ######################################################################################################
######################################################################################################################################################
patches <- st_read(here("./MartenPatchPolygons_NoAreaThresh.shp"))
mapview(patches)
mapview(patches[patches$patches == 1,])

## Calculate the area in m2 and hectare for book keeping
patches <- patches |>
  mutate(Area_m2 = st_area(geometry),
         Area_ha = units::set_units(st_area(geometry), "ha"))

## Snow-free entropy raster
filelist.ent <- list.files("D:/My Drive/Entropy/", full.names = TRUE)
dates <- seq(2013, 2022, by = 1) # -- Ran in 5 year groups because of ram issues -- ##

ent <- rast()
for(i in 1:length(dates)){
  temp.list <- as.list(filelist.ent[str_detect(filelist.ent, as.character(dates[i]))])
  temp.list <- lapply(temp.list, rast)
  rsrc <- sprc(temp.list)
  temp.merge <- merge(rsrc)
  print(paste("Merge Finished for year", dates[i]))
  temp.merge[temp.merge < 3.24] <- NA # Mask values outside observed
  names(temp.merge) <- dates[i]
  ent <- c(ent, temp.merge)
  print(paste("File finished for year", dates[i]))
  gc()
}
ent
names(ent) <- dates

## Annual snow depth rasters
filelist.snow <- list.files("D:/My Drive/SnowDepth/", full.names = TRUE)
snow <- rast(filelist.snow)
snow <- snow * 1000 # convert to mm
names(snow) <- seq(2003,2022, by = 1)

snow
#mapview(snow$`2003`) + mapview(patches[patches$patches == 109,])

## Future snowdepth rasters
A1mid <- rast(here("./Data/MidCentury_PctChange.tif"))
A1mid <- A1mid/100
A1late <- rast(here("./Data/LateCentury_PctChange.tif"))
A1late <- A1late/100

## Future Winter Lenght rasters
wl.baseline <- rast(here("./Data/BaselineWinterLength.tif"))
wl.A1mid <- rast(here("./Data/Mid21WinterLength.tif"))
wl.A1late <- rast(here("./Data/Late21WinterLength.tif"))

## Scale these rasters to be a proportion of the full year with snow
## These values need to be supplied the scale the demographic rates
## with the weighted mean for annual survival as
## (surv_so * p) + (surv_sf * (1-p)) where surv_so = survival snow-on
## surv_sf = survival snow-off, and p = proportion of days with snow
wl.baseline <- wl.baseline/365
wl.A1mid <- wl.A1mid/365
wl.A1late <- wl.A1late/365

## Reproject into the snow layer
A1mid <- project(A1mid, snow)
A1mid <- resample(A1mid, snow)
A1late <- project(A1late, snow)
A1late <- resample(A1late, snow)

snow.c <- crop(snow, A1mid)

## Scale the snow cover
# Create a list to store results
future_snow <- list()

# For each percent change scenario
rcm_pct_change_fun <- function(future,
                               current){
  future_snow <- list()

  for(i in 1:nlyr(future)) {
    # Get the percent change for this scenario
    pct <- future[[i]]

    # Apply the change to the snow values
    # When pct is negative, it will decrease the snow values
    # When pct is positive, it will increase the snow values
    adjusted_snow <- current * (1 + pct)

    # Store in list
    names(adjusted_snow) <- paste0(names(adjusted_snow), "_", names(pct))
    future_snow[[i]] <- adjusted_snow
  }
  return(future_snow)
}

# Combine all results
mid_list <- rcm_pct_change_fun(future = A1mid,
                          current = snow.c)

late_list <- rcm_pct_change_fun(future = A1late,
                                current = snow.c)

########################################################################################################################################################
################################## Estimate survival ###################################################################################################
########################################################################################################################################################
### Estimate snow-free survival

## transform values based on observed and test survival
ent.std <- terra::app(ent, function(x) (x-4.14)/0.18) ## standardize raster to observed data
gc()

## Let's try to reduce the extent for faster processing
mart.ext <- patches |>
  st_transform(crs = crs(ent.std)) |>
  st_bbox() |>
  st_as_sfc() |>
  st_sf() |>
  vect()

ent.std.c <- crop(ent.std, mart.ext)
gc()
plot(ent.std.c[[1]])
#writeRaster(ent.std.c, filename = "G:/Entropy2013_2022_Standardized_Cropped.tif")
ent.std.c <- rast("G:/Entropy2013_2022_Standardized_Cropped.tif")
ent.surv <- rast()

ent.std.c1 <- ent.std.c[[1:5]]
ent.surv <- terra::app(ent.std.c1, function(x) 0.87^(exp(x*(-1.8984))))

for(i in 1:nlyr(ent.std.c)){
  print(paste("Ras", i))
  tmp.r <- ent.std.c[[i]]
  tmp.surv <- terra::app(tmp.r, function(x) 0.87^(exp(x*(-1.8984)))) ## standardized betas
  ent.surv <- c(ent.surv, tmp.surv)
  gc()
}

names(ent.surv) <- names(ent.std.c)

ent.surv.patches <- terra::extract(ent.surv, patches, mean, na.rm = TRUE) ## extract survival by patch
patches.ent <- patches %>%
  st_drop_geometry %>%
  bind_cols(ent.surv.patches) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_ent")

#write.csv(patches.ent, here("./Data/marten_patch_survival/ent_surv_2013_2022_long_AllPatches.csv"), row.names = FALSE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Snow-on
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mid_list
late_list

### Estimate snow-on survival

surv_proj_fun <- function(current,
                          future){
for(i in 1:nlyr(snow)){
  current[[i]][current[[i]] > 1172] <- NA # Mask values outside observed
  future[[i]][future[[i]] > 1172] <- NA # Mask values outside observed
}

snow.std <- terra::app(current, function(x) (x-662)/301) ## standardize raster to observed data
snow.surv <- app(snow.std, function(x) 0.86^(exp(x*(-2.509)))) # estimate survival from standardized betas
snowA1mid.std <- terra::app(future, function(x) (x-662)/301) ## standardize raster to observed data
snowA1mid.surv <- app(snowA1mid.std, function(x) 0.86^(exp(x*(-2.509)))) # estimate survival from standardized betas

snow.out <- list(
  current.surv = snow.surv,
  future.surv = snowA1mid.surv)

return(snow.out)

}

## Apply this function for each of the
## RCMS
snow_mid_cent <- lapply(mid_list, function(x) surv_proj_fun(current = snow.c, future = x))
snow_late_cent <- lapply(late_list, function(x) surv_proj_fun(current = snow.c, future = x))

plot(snow_mid_cent[[1]]$current.surv)
plot(snow_late_cent[[1]]$current.surv)
plot(snow_mid_cent[[1]]$future.surv)
plot(snow_late_cent[[1]]$future.surv)

## Extract the values
snow.surv.patches <- terra::extract(snow_mid_cent[[1]]$current.surv, patches, mean, na.rm = TRUE)

snow.surv.patches.midAccess <- terra::extract(snow_mid_cent[[1]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.midCnrm <- terra::extract(snow_mid_cent[[2]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.midGfdl <- terra::extract(snow_mid_cent[[3]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.midIpsl <- terra::extract(snow_mid_cent[[4]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.midMiroc5 <- terra::extract(snow_mid_cent[[5]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.midMri <- terra::extract(snow_mid_cent[[6]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch

snow.surv.patches.lateAccess <- terra::extract(snow_late_cent[[1]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.lateCnrm <- terra::extract(snow_late_cent[[2]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.lateGfdl <- terra::extract(snow_late_cent[[3]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.lateIpsl <- terra::extract(snow_late_cent[[4]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.lateMiroc5 <- terra::extract(snow_late_cent[[5]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.lateMri <- terra::extract(snow_late_cent[[6]]$future.surv, patches, mean, na.rm = TRUE) # extract survival by patch

# Combine all datasets
patches.snow <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(scenario = "Current")

patches.snow.midAccess <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.midAccess) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Mid_Access")

patches.snow.midCnrm <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.midCnrm) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Mid_CNRM")

patches.snow.midGfdl <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.midGfdl) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Mid_GFDL")

patches.snow.midIpsl <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.midIpsl) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Mid_IPSL")

patches.snow.midMiroc5 <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.midMiroc5) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Mid_MIROC5")

patches.snow.midMri <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.midMri) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Mid_MRI")

## Late century
patches.snow.lateAccess <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.lateAccess) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Late_Access")

patches.snow.lateCnrm <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.lateCnrm) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Late_CNRM")

patches.snow.lateGfdl <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.lateGfdl) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Late_GFDL")

patches.snow.lateIpsl <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.lateIpsl) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Late_IPSL")

patches.snow.lateMiroc5 <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.lateMiroc5) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Late_MIROC5")

patches.snow.lateMri <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.lateMri) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(year = stringr::str_extract(year, "^\\d{4}")) |>
  mutate(scenario = "Late_MRI")

patches.snow <- rbind(patches.snow,
                      patches.snow.midAccess,
                      patches.snow.midCnrm,
                      patches.snow.midGfdl,
                      patches.snow.midIpsl,
                      patches.snow.midMiroc5,
                      patches.snow.midMri,
                      patches.snow.lateAccess,
                      patches.snow.lateCnrm,
                      patches.snow.lateGfdl,
                      patches.snow.lateIpsl,
                      patches.snow.lateMiroc5,
                      patches.snow.lateMri)
head(patches.snow)
tail(patches.snow)
#write.csv(snow.surv.patches, "~/Desktop/snow_surv_2003-2022.csv", row.names = FALSE)
write.csv(patches.snow, here("./Data/marten_patch_survival/snow_surv_2003_2022_long_CurrentFuture_AllPatches_Projections.csv"), row.names = FALSE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Winter Length Extraction
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Rasters of duration
## Future Winter Lenght rasters
wl.baseline <- rast(here("./Data/BaselineWinterLength.tif"))
wl.A1mid <- rast(here("./Data/Mid21WinterLength.tif"))
wl.A1late <- rast(here("./Data/Late21WinterLength.tif"))

## Scale these rasters to be a proportion of the full year with snow
## These values need to be supplied the scale the demographic rates
## with the weighted mean for annual survival as
## (surv_so * p) + (surv_sf * (1-p)) where surv_so = survival snow-on
## surv_sf = survival snow-off, and p = proportion of days with snow
wl.baseline <- wl.baseline/365
wl.A1mid <- wl.A1mid/365
wl.A1late <- wl.A1late/365

## Interpolating into the lakes to deal with edge effects
wl.baseline_focal <- lapply(wl.baseline, function(x){
  focal(x, 7, "mean", na.policy = "only", na.rm = T)
})
wl.A1mid_focal <- lapply(wl.A1mid, function(x){
  focal(x, 7, "mean", na.policy = "only", na.rm = T)
})
wl.A1late_focal <- lapply(wl.A1late, function(x){
  focal(x, 7, "mean", na.policy = "only", na.rm = T)
})

## Resample to 1 km
template <- rast(ext(wl.baseline_focal[[1]]), res=1000) # res in meters
crs(template) <- crs(wl.baseline_focal[[1]])

wl.baseline_km <- lapply(wl.baseline_focal, function(x){ resample(x, template, method = "bilinear") })
wl.A1mid_km <- lapply(wl.A1mid_focal, function(x){ resample(x, template, method = "bilinear") })
wl.A1late_km <- lapply(wl.A1late_focal, function(x){ resample(x, template, method = "bilinear") })


## Reproject into the snow layer
names(wl.baseline_km) <- names(wl.baseline)
wl.baseline_km <- project(rast(wl.baseline_km), snow)
wl.baseline_km <- resample(wl.baseline_km, snow)
names(wl.A1mid_km) <- names(wl.A1mid)
wl.A1mid_km <- project(rast(wl.A1mid_km), snow)
wl.A1mid_km <- resample(wl.A1mid_km, snow)
names(wl.A1late_km) <- names(wl.A1late)
wl.A1late_km <- project(rast(wl.A1late_km), snow)
wl.A1late_km <- resample(wl.A1late_km, snow)


wlength.patches.baseAccess <- terra::extract(wl.baseline_km[[1]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.baseCnrm <- terra::extract(wl.baseline_km[[2]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.baseGfdl <- terra::extract(wl.baseline_km[[3]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.baseIpsl <- terra::extract(wl.baseline_km[[4]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.baseMiroc5 <- terra::extract(wl.baseline_km[[5]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.baseMri <- terra::extract(wl.baseline_km[[6]], patches, mean, na.rm = TRUE) # extract survival by patch

## Calculate the average for the baseline since we don't want GCM specific estimates but want to
## respect scaling in the contemporary period
wlength.patches.base <- wlength.patches.baseAccess |>
  left_join(wlength.patches.baseCnrm, by = "ID") |>
  left_join(wlength.patches.baseGfdl, by = "ID") |>
  left_join(wlength.patches.baseIpsl, by = "ID") |>
  left_join(wlength.patches.baseMiroc5, by = "ID") |>
  left_join(wlength.patches.baseMri, by = "ID") |>
  rowwise() |>
  mutate(Mean_WL = mean(c_across(-ID))) |>
  ungroup() |>
  dplyr::select(ID, Mean_WL)


wlength.patches.midAccess <- terra::extract(wl.A1mid_km[[1]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.midCnrm <- terra::extract(wl.A1mid_km[[2]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.midGfdl <- terra::extract(wl.A1mid_km[[3]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.midIpsl <- terra::extract(wl.A1mid_km[[4]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.midMiroc5 <- terra::extract(wl.A1mid_km[[5]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.midMri <- terra::extract(wl.A1mid_km[[6]], patches, mean, na.rm = TRUE) # extract survival by patch

wlength.patches.lateAccess <- terra::extract(wl.A1late_km[[1]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.lateCnrm <- terra::extract(wl.A1late_km[[2]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.lateGfdl <- terra::extract(wl.A1late_km[[3]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.lateIpsl <- terra::extract(wl.A1late_km[[4]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.lateMiroc5 <- terra::extract(wl.A1late_km[[5]], patches, mean, na.rm = TRUE) # extract survival by patch
wlength.patches.lateMri <- terra::extract(wl.A1late_km[[6]], patches, mean, na.rm = TRUE) # extract survival by patch

# Combine all datasets
patches.wl <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.base) %>%
  mutate(scenario = "Current") |>
  rename(WL = Mean_WL)

patches.wl.midAccess <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.midAccess) %>%
  mutate(scenario = "Mid_Access") |>
  rename(WL = access_mid21)

patches.wl.midCnrm <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.midCnrm) %>%
  mutate(scenario = "Mid_CNRM") |>
  rename(WL = cnrm_mid21)

patches.wl.midGfdl <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.midGfdl) %>%
  mutate(scenario = "Mid_GFDL") |>
  rename(WL = gfdl_mid21)


patches.wl.midIpsl <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.midIpsl) %>%
  mutate(scenario = "Mid_IPSL") |>
  rename(WL = ipsl_mid21)

patches.wl.midMiroc5 <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.midMiroc5) %>%
  mutate(scenario = "Mid_MIROC5") |>
  rename(WL = miroc5_mid21)


patches.wl.midMri <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.midMri) %>%
  mutate(scenario = "Mid_MRI") |>
  rename(WL = mri_mid21)


## Late century
patches.wl.lateAccess <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.lateAccess) %>%
  mutate(scenario = "Late_Access") |>
  rename(WL = access_late21)


patches.wl.lateCnrm <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.lateCnrm) %>%
  mutate(scenario = "Late_CNRM") |>
  rename(WL = cnrm_late21)


patches.wl.lateGfdl <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.lateGfdl) %>%
  mutate(scenario = "Late_GFDL") |>
  rename(WL = gfdl_late21)


patches.wl.lateIpsl <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.lateIpsl) %>%
  mutate(scenario = "Late_IPSL") |>
  rename(WL = ipsl_late21)


patches.wl.lateMiroc5 <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.lateMiroc5) %>%
  mutate(scenario = "Late_MIROC5") |>
  rename(WL = miroc5_late21)


patches.wl.lateMri <- patches %>%
  st_drop_geometry %>%
  bind_cols(wlength.patches.lateMri) %>%
  mutate(scenario = "Late_MRI") |>
  rename(WL = mri_late21)


patches.wl <- rbind(patches.wl,
                    patches.wl.midAccess,
                    patches.wl.midCnrm,
                    patches.wl.midGfdl,
                    patches.wl.midIpsl,
                    patches.wl.midMiroc5,
                    patches.wl.midMri,
                    patches.wl.lateAccess,
                    patches.wl.lateCnrm,
                    patches.wl.lateGfdl,
                    patches.wl.lateIpsl,
                    patches.wl.lateMiroc5,
                    patches.wl.lateMri) |>
  dplyr::select(patches, ID, scenario, WL)

#########################################################################################################################################################
################################## Summary ##############################################################################################################
#########################################################################################################################################################
## Removed SpatRasters and read in csv files
## I had to split up the entropy rasters into 5 year groups for processing. Combine into single file and join with snow data
patches.snow <- read.csv(here("./Data/marten_patch_survival/snow_surv_2003-2022_long_CurrentFuture.csv")) %>% mutate(year = as.character(year))
patches.snow <- read.csv(here("./Data/marten_patch_survival/snow_surv_2003_2022_long_CurrentFuture_AllPatches_Projections.csv")) |>
  mutate(year = as.character(year))

patches.ent <- read.csv("./Data/marten_patch_survival/ent_surv_2013_2022_long_AllPatches.csv") %>% mutate(year = as.character(year))

## Add some columns to the ent DF so we can make the naming consistent
patches.ent <- bind_rows(
  patches.ent %>% mutate(scenario = "Current"),
  patches.ent %>% mutate(scenario = "Mid_Access"),
  patches.ent %>% mutate(scenario = "Mid_CNRM"),
  patches.ent %>% mutate(scenario = "Mid_GFDL"),
  patches.ent %>% mutate(scenario = "Mid_IPSL"),
  patches.ent %>% mutate(scenario = "Mid_MIROC5"),
  patches.ent %>% mutate(scenario = "Mid_MRI"),
  patches.ent %>% mutate(scenario = "Late_Access"),
  patches.ent %>% mutate(scenario = "Late_CNRM"),
  patches.ent %>% mutate(scenario = "Late_GFDL"),
  patches.ent %>% mutate(scenario = "Late_IPSL"),
  patches.ent %>% mutate(scenario = "Late_MIROC5"),
  patches.ent %>% mutate(scenario = "Late_MRI")
)

## Join with snow survival to estiamte annual average
patches.join <- left_join(patches.snow, patches.ent, by = c("patches", "year", "scenario", "ID")) %>%
  filter(as.numeric(year) >= 2013) %>%
  group_by(scenario) %>%  # Add grouping if you want means by scenario
  mutate(
    surv_ent = ifelse(is.nan(surv_ent), mean(surv_ent, na.rm = TRUE), surv_ent),
    surv_snow = ifelse(is.nan(surv_snow), mean(surv_snow, na.rm = TRUE), surv_snow)  # Fixed typo: - to =
  ) %>%
  ungroup() %>%
  left_join(patches.wl, by = c("patches", "ID", "scenario")) |>
  mutate(surv_annu = (surv_ent + surv_snow)/2) |>
  mutate(surv_annu_wladj = (surv_snow * WL) + (surv_ent * (1-WL)))


## CI function
calculate_confidence_interval <- function(mean, sd, n, confidence = 0.95) {
  error_margin <- qt((1 + confidence) / 2, df = n - 1) * (sd / sqrt(n))
  lower_bound <- round(mean - error_margin,2)
  upper_bound <- round(mean + error_margin,2)
  return(c(lower = lower_bound, upper = upper_bound))
}

## Summarize by patch over 20 years 2003-2022 ## esimated median from time series
patches.surv.sum <- patches.join %>%
  group_by(patches, scenario) %>%
  summarise(surv_ent_20y = round(median(surv_ent, na.rm = T),2),
            surv_snow_20y = round(median(surv_snow, na.rm = T),2),
            surv_annu_20y = round(median(surv_annu, na.rm = T),2),
            surv_annu_20y_adj = round(median(surv_annu_wladj, na.rm = T),2),
            sd_annu_20y = round(sd(surv_annu, na.rm = T),2), n = n(),
            lower_CI_annu = calculate_confidence_interval(surv_annu_20y, sd_annu_20y, n)[1],
            upper_CI_annu = calculate_confidence_interval(surv_annu_20y, sd_annu_20y, n)[2],
            sd_annu_20y_adj = round(sd(surv_annu_wladj, na.rm = T),2), n = n(),
            lower_CI_annu_adj = calculate_confidence_interval(surv_annu_20y_adj, sd_annu_20y_adj, n)[1],
            upper_CI_annu_adj = calculate_confidence_interval(surv_annu_20y_adj, sd_annu_20y_adj, n)[2],
            sd_ent_20y = round(sd(surv_ent, na.rm = T),2),
            sd_snow_20y = round(sd(surv_snow, na.rm = T),2),
            lower_CI_ent = calculate_confidence_interval(surv_ent_20y, sd_ent_20y, n)[1],
            upper_CI_ent = calculate_confidence_interval(surv_ent_20y, sd_ent_20y, n)[2],
            lower_CI_snow = calculate_confidence_interval(surv_snow_20y, sd_ent_20y, n)[1],
            upper_CI_snow = calculate_confidence_interval(surv_snow_20y, sd_ent_20y, n)[2])

# View the summary statistics with confidence interval
head(patches.surv.sum)
tail(patches.surv.sum)
write.csv(patches.surv.sum[,1:9], here("./Data/marten_patch_survival/all_patches_annual_survival_current_future_projections_2013_2022_wl_adjusted.csv"), row.names = FALSE)

#########################################################################################################################################################
################################## Plot #################################################################################################################
#########################################################################################################################################################
patches.join.long <- patches.join %>% pivot_longer(cols = starts_with("surv"), names_to = "Season", values_to = "Survival")

ggplot() + geom_boxplot(data = patches.join.long, aes(x = factor(patches), y = Survival, color = Season)) +
  theme_bw(base_size = 30) + theme(legend.position = "top")
ggplot() + geom_boxplot(data = patches.join.long[patches.join.long$Season == "surv_annu",], aes(x = factor(patches), y = Survival, color = Season)) +
  theme_bw(base_size = 30) + theme(legend.position = "top")

