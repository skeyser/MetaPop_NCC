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
A1mid <- (1 - 0.3)
A1late <- (1 - 0.5)

snowA1mid <- snow * A1mid
snowA1late <- snow * A1late

plot(snow[[1]])
plot(snowA1mid[[1]])
plot(snowA1late[[1]])

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
### Estimate snow-on survival
for(i in 1:nlyr(snow)){
  snow[[i]][snow[[i]] > 1172] <- NA # Mask values outside observed
  snowA1mid[[i]][snowA1mid[[i]] > 1172] <- NA # Mask values outside observed
  snowA1late[[i]][snowA1late[[i]] > 1172] <- NA # Mask values outside observed
}

snow.std <- terra::app(snow, function(x) (x-662)/301) ## standardize raster to observed data
snow.surv <- app(snow.std, function(x) 0.86^(exp(x*(-2.509)))) # estimate survival from standardized betas
snowA1mid.std <- terra::app(snowA1mid, function(x) (x-662)/301) ## standardize raster to observed data
snowA1mid.surv <- app(snowA1mid.std, function(x) 0.86^(exp(x*(-2.509)))) # estimate survival from standardized betas
snowA1late.std <- terra::app(snowA1late, function(x) (x-662)/301) ## standardize raster to observed data
snowA1late.surv <- app(snowA1late.std, function(x) 0.86^(exp(x*(-2.509)))) # estimate survival from standardized betas

snow.surv.patches <- terra::extract(snow.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.a1mid <- terra::extract(snowA1mid.surv, patches, mean, na.rm = TRUE) # extract survival by patch
snow.surv.patches.a1late <- terra::extract(snowA1late.surv, patches, mean, na.rm = TRUE) # extract survival by patch

# Combine all datasets
patches.snow <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(scenario = "Current")

patches.snow.a1mid <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.a1mid) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(scenario = "A1mid")

patches.snow.a1late <- patches %>%
  st_drop_geometry %>%
  bind_cols(snow.surv.patches.a1late) %>%
  pivot_longer(cols = starts_with("20"), names_to = "year", values_to = "surv_snow") |>
  mutate(scenario = "A1late")

patches.snow <- rbind(patches.snow, patches.snow.a1mid, patches.snow.a1late)

#write.csv(snow.surv.patches, "~/Desktop/snow_surv_2003-2022.csv", row.names = FALSE)
write.csv(patches.snow, here("./Data/marten_patch_survival/snow_surv_2003_2022_long_CurrentFuture_AllPatches.csv"), row.names = FALSE)

#########################################################################################################################################################
################################## Summary ##############################################################################################################
#########################################################################################################################################################
## Removed SpatRasters and read in csv files
## I had to split up the entropy rasters into 5 year groups for processing. Combine into single file and join with snow data
patches.snow <- read.csv(here("./Data/marten_patch_survival/snow_surv_2003-2022_long_CurrentFuture.csv")) %>% mutate(year = as.character(year))
patches.ent <- read.csv("./Data/marten_patch_survival/ent_surv_2003-2022_long.csv") %>% mutate(year = as.character(year))

## Add some columns to the ent DF so we can make the naming consistent
patches.ent <- bind_rows(
  patches.ent %>% mutate(scenario = "Current"),
  patches.ent %>% mutate(scenario = "A1mid"),
  patches.ent %>% mutate(scenario = "A1late")
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
  mutate(surv_annu = (surv_ent + surv_snow)/2)

## CI function
calculate_confidence_interval <- function(mean, sd, n, confidence = 0.95) {
  error_margin <- qt((1 + confidence) / 2, df = n - 1) * (sd / sqrt(n))
  lower_bound <- round(mean - error_margin,2)
  upper_bound <- round(mean + error_margin,2)
  return(c(lower = lower_bound, upper = upper_bound))
}

## Summarize by patch over 20 years 2003-2022 ## esimated median from time series
patches.surv.sum <- patches.join %>% group_by(patches, scenario) %>%
  summarise(surv_ent_20y = round(median(surv_ent),2), surv_snow_20y = round(median(surv_snow),2),
            surv_annu_20y = round(median(surv_annu),2), sd_annu_20y = round(sd(surv_annu),2), n = n(),
            lower_CI_annu = calculate_confidence_interval(surv_annu_20y, sd_annu_20y, n)[1],
            upper_CI_annu = calculate_confidence_interval(surv_annu_20y, sd_annu_20y, n)[2],
            sd_ent_20y = round(sd(surv_ent),2), sd_snow_20y = round(sd(surv_snow),2),
            lower_CI_ent = calculate_confidence_interval(surv_ent_20y, sd_ent_20y, n)[1],
            upper_CI_ent = calculate_confidence_interval(surv_ent_20y, sd_ent_20y, n)[2],
            lower_CI_snow = calculate_confidence_interval(surv_snow_20y, sd_ent_20y, n)[1],
            upper_CI_snow = calculate_confidence_interval(surv_snow_20y, sd_ent_20y, n)[2])

# View the summary statistics with confidence interval
print(patches.surv.sum)
write.csv(patches.surv.sum[,1:9], here("./Data/marten_patch_survival/all_patches_annual_survival_current_future_2013_2022.csv"), row.names = FALSE)

#########################################################################################################################################################
################################## Plot #################################################################################################################
#########################################################################################################################################################
patches.join.long <- patches.join %>% pivot_longer(cols = starts_with("surv"), names_to = "Season", values_to = "Survival")

ggplot() + geom_boxplot(data = patches.join.long, aes(x = factor(patches), y = Survival, color = Season)) +
  theme_bw(base_size = 30) + theme(legend.position = "top")
ggplot() + geom_boxplot(data = patches.join.long[patches.join.long$Season == "surv_annu",], aes(x = factor(patches), y = Survival, color = Season)) +
  theme_bw(base_size = 30) + theme(legend.position = "top")
