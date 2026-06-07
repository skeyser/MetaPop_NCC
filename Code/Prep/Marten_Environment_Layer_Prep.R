## ----------------------------------------------------------
##
## Script name: Marten SDM Environmental Variable Prep Script
##
## Script purpose: Prep all environmental variables for SDMs @ 1 km resolution
## for mammal SDMs.
##
## Author: Spencer R Keyser
##
## Date Created: 2023-09-29
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
library(terra)
library(sf)
library(rasterVis)

## -----------------------------------------------------------

## -----------------------------------------------------------
##
## Begin Section: Defining the region for predictor creation
##
## -----------------------------------------------------------

## ***********************************************************
##
## Section Notes: Things to consider
## 1. Extent of all the known marten samples
## 2. The cost layer constrains the minimum size of the domain
##  i. Modelling domain > model projection
##
## ***********************************************************

## Load the cost surface for bounding box
cost <- rast("R:/Users/skeyser/PhD/NASA Project/Spencer connectivity data/Circuitscape/Input/predcosts.txt")

## Take box for cost layer (buffer a little but for edges)
bbox <- st_bbox(cost)
cost.box <- st_polygon(list(
    cbind(
      c(floor(bbox$xmin),
      ceiling(bbox$xmax),
      ceiling(bbox$xmax),
      floor(bbox$xmin),
      floor(bbox$xmin))
    ,
    c(floor(bbox$ymin),
      floor(bbox$ymin),
      ceiling(bbox$ymax),
      ceiling(bbox$ymax),
      floor(bbox$ymin))
  ))) |>
  st_sfc(crs = crs(cost))

plot(cost.box)

dom <- vect(cost.box)

st_write(cost.box, "R:/Users/skeyser/PhD/NASA Project/Marten_Data/CostLayerBBox.shp")

## Region for cropping
# dom.og <- st_read("G:/eBird/US_Map/states.shp") |>
#   filter(STATE_NAME %in% c("Minnesota", "Wisconsin", "Michigan")) |>
#   st_union() |>
#   st_transform(crs = 4326)
#
# dom <- dom.og |>
#   st_transform(crs = "EPSG:5070") |>
#   st_buffer(dist = 10000) |>
#   st_transform(crs = 4326)
#
# wi <- st_read("G:/eBird/US_Map/states.shp") |>
#   filter(STATE_NAME == "Wisconsin") |>
#   st_buffer(dist = 10000) |>
#   st_transform(crs = 4326)

## Dom
#dom <- st_read("C:/Users/skeyser/Documents/ArcGIS/Projects/Marten_Habitat_Layers/GreatLakesBBox.shp")
#dom <- vect(dom)

dom.tmp <- project(dom, "EPSG:5070")

ggplot() +
  geom_sf(data = dom.og, color = "red", fill = NA) +
  geom_sf(data = dom, fill = NA)

## SpatVector Transforms
dom.og <- vect(dom.og)

dom <- vect(dom)

## Raster Template
r.tmp <- rast(xmin = xmin(dom.tmp),
              xmax = xmax(dom.tmp),
              ymin = ymin(dom.tmp),
              ymax = ymax(dom.tmp),
              resolution = 1000,
              crs = crs(dom.tmp))


## -----------------------------------------------------------
##
## Begin Section: Relevant Environmental Predictors
##
## -----------------------------------------------------------

## ***********************************************************
##
## Section Notes: Environmental Predictor List
## 1. Bioclim variables
## 2. WHIs
## 3. Snow depth
## 4. Elevation
## 5. Human Modification Index
## 6. Landcover
## 7. Tree height
##
##
## ***********************************************************

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Bioclim variables
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ***********************************************************
##
## Section Notes:
## BIO1: Annual Mean Temperature
## *BIO2: Mean Diurnal Range (Mean of monthly (max T - min T))
## BIO3: Isothermality (BIO2/BIO7)(*100)
## *BIO4: Temperature Seasonality (SD * 100)
## *BIO5: Max Temp of Warmest Month
## *BIO6: Min Temp of Coldest Month
## BIO7: Temperature Annual Range (BIO5 - BIO6)
## BIO8: Mean Temperature of Wettest Quarter
## BIO9: Mean Temperature of Driest Quarter
## BIO10: Mean Temperature of Warmest Quarter
## BIO11: Mean Temperature of Coldest Quarter
## BIO12: Annual Precipitation
## BIO13: Precipitation of Wettest Month
## BIO14: Precipitation of Driest Month
## BIO15: Precipitation Seasonality (CV)
## BIO16: Precipitation of Wettest Quarter
## BIO17: Precipitation of Driest Quarter
## *BIO18: Precipitation of Warmest Quarter
## *BIO19: Precipitation of Coldest Quarter
##
## *denotes the Bioclim vars used for Bobcat and Lynx Peers et al., 2013
## ***********************************************************

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Water masking
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fill in the missing NAs
water.mask <- rast("E:/My Drive/GEE_Output/Water_Mask500.tif")
plot(water.mask)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: CHELSA
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dom.wgs84 <- project(dom, "EPSG:4326")

chel <- stringr::str_sort(list.files(path = "G:/CHELSA_BioClim/",
                                     pattern = "*.tif$",
                                     full.names = TRUE),
                          numeric = T)
chel <- rast(chel)
chel <- crop(chel, dom.wgs84, mask = T)

names(chel) <- c("MAT",
                 "MDR",
                 "Iso",
                 "TSea",
                 "TMaxWM",
                 "TMinCM",
                 "TR",
                 "MATWetQ",
                 "MATDryQ",
                 "MATWarmQ",
                 "MATColdQ",
                 "PcpAnnual",
                 "PcpWetM",
                 "PcpDryM",
                 "PcpSea",
                 "MeanPcpWetQ",
                 "MeanPcpDryQ",
                 "MeanPcpWarmQ",
                 "MeanPcpDryQ",
                 "SCD")

plot(chel[["MATColdQ"]])

## Take only a couple variables from CHELSA Climate Normals (1980-2010)
chel <- chel[[c("MATColdQ", "MATWarmQ", "MAT", "PcpAnnual")]]

## Resample
chel <- project(chel, crs(r.tmp))
chel <- resample(chel, r.tmp, method = "bilinear")
plot(chel)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Elevation
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ele <- rast("V:/GudexCross/GMTED_GlobalElev/GMTED2010_Elevation_Global.tif")
ele <- crop(ele, dom.wgs84, mask = T)
ele <- project(ele, crs(r.tmp))
ele <- resample(ele, r.tmp, method = "bilinear")
plot(ele)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: TPI & TR
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tpi <- rast("G:/eBird/DEM/tpi_1KMmd_GMTEDmd.tif")
tpi <- crop(tpi, dom.wgs84, mask = T)
tpi <- project(tpi, crs(r.tmp))
tpi <- resample(tpi, r.tmp, method = "bilinear")
plot(tpi)

tr <- rast("G:/eBird/DEM/roughness_1KMmd_GMTEDmd.tif")
tr <- crop(tr, dom.wgs84, mask = T)
tr <- project(tr, crs(r.tmp))
tr <- resample(tr, r.tmp, method = "bilinear")
plot(tr)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Human Modification Index
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hum <- rast("G:/HumanMod/gHM/gHM/gHM.tif")
dom.hum <- project(dom, crs(hum))
hum <- crop(hum, dom.hum, mask = T)
hum <- project(hum, crs(r.tmp))
hum <- resample(hum, r.tmp, method = "bilinear")

plot(hum)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Forest Height
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fh <- rast("G:/ForestHeight/Forest_height_2019_NAM.tif")

fh <- crop(fh, dom.wgs84)
fh <- mask(fh, dom.wgs84)
## Reclassify the fh product to deal with alternative classes
## 0-60 FH
## 101 water
## 102 snow/ice
## 103 No DATA
fh <- classify(fh, cbind(61, 103, NA))
fh <- project(fh, crs(r.tmp))
fh <- resample(fh, r.tmp, method = "bilinear")
names(fh) <- "For_Hgt"
plot(fh)

## Match extents
# ext(water.mask) <- ext(fh)
# fh.c <- crop(fh, water.mask)
# wm.fh <- resample(water.mask, fh.c, "near")
# fh.c <- mask(fh.c, wm.fh, maskvalues = 1)
# plot(fh.c)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Canopy Cover
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
can <- rast("E:/My Drive/GEE_Output/CostLayerCanopyCover.tif")
dom.can <- project(dom.wgs84, crs(can))
can <- crop(can, dom.can, mask = F)
can <- mask(can, dom.can)
can <- project(can, crs(r.tmp))
can <- resample(can, r.tmp, method = "bilinear")
plot(can)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: WHI
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load snow season length backfilled data for 2003-2020
ssl <- rast("G:/WHI/Backfilled/SSL/SSL_BF_20032020.tif")
plot(ssl)

ssl <- crop(ssl, dom.wgs84, mask = T)
plot(ssl)

## Remove the NAs caused by the backfilling
#ssl.f <- focal(ssl, w = 9, fun = "mean", na.policy = "only", fillvalue = NA)
ssl.f <- project(ssl, crs(r.tmp))
ssl.f <- resample(ssl.f, r.tmp, method = "bilinear")
plot(ssl.f)

# ## Match extents
# ext(water.mask) <- ext(ssl.f)
# water.mask <- resample(water.mask, ssl.f, method = "near")
#
# ## Mask by water pixels
# ssl.f <- mask(ssl.f, water.mask, maskvalues = 1)
# plot(ssl.f)
names(ssl.f) <- "SSL"

## Repeat for all WHIs
fwos <- rast("G:/WHI/Backfilled/FWOS/FWOS_BF_20032020.tif")
fwos <- crop(fwos, dom.wgs84, mask = T)
#fwos.f <- focal(fwos, w = 9, fun = "mean", na.policy = "only", fillvalue = NA)
fwos.f <- project(fwos, crs(r.tmp))
fwos.f <- resample(fwos.f, r.tmp, method = "bilinear")
plot(fwos.f)

# fwos.f <- mask(fwos.f, water.mask, maskvalues = 1)
names(fwos.f) <- "FWOS"

## SCV
scv <- rast("G:/WHI/Backfilled/SCv/SCV_BF_20032020.tif")
scv <- crop(scv, dom.wgs84, mask = T)
#scv.f <- focal(scv, w = 9, fun = "mean", na.policy = "only", fillvalue = NA)
scv.f <- project(scv, crs(r.tmp))
scv.f <- resample(scv.f, r.tmp, method = "bilinear")
plot(scv.f)

# scv.f <- mask(scv.f, water.mask, maskvalues = 1)
# plot(scv.f)
names(scv.f) <- "SCV"

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: LULC %
## % Evergreen, Mixed, Deciduous Forest, Wetland, Open, Urban, and Cropland
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dom <- cost.box |>
  st_transform(crs = crs(lc)) |>
  vect()

if(!file.exists("R:/Users/skeyser/PhD/NASA Project/Marten_Data/PctLandCover_NALCMS/")){
  lc <- rast("G:/Land_cover_2020_30m_TIF/NA_NALCMS_landcover_2020_30m/data/NA_NALCMS_landcover_2020_30m.tif")

  dom.proj <- project(dom, crs(lc))
  lc <- crop(lc, dom.proj, mask = F)
  lc <- mask(lc, dom.proj)
  plot(lc)

  pclass.needle <- function(x, y=c(1,2)) {
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.dec <- function(x, y=c(4,5)) {
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.mix <- function(x, y=c(6)) {
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.wet <- function(x, y=c(14)){
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.grass <- function(x, y=c(9, 10, 12)){
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.shrub <- function(x, y=c(7, 8, 11)){
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.barren <- function(x, y=c(13, 16)){
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.urban <- function(x, y=c(17)){
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.crop <- function(x, y=c(15)){
    return( length(which(x %in% y)) / length(x) )
  }

  ## Broad cats
  pclass.for <- function(x, y=c(1,2,4,5,6)){
    return( length(which(x %in% y)) / length(x) )
  }
  pclass.open <- function(x, y=c(8,10,11,12,13,16)){
    return( length(which(x %in% y)) / length(x) )
  }

  ## Specific cats
  pf.nl <- terra::focal(lc, w=matrix(1, 9, 9), pclass.needle)
  plot(pf.nl)
  pf.dec <- terra::focal(lc, w=matrix(1, 9, 9), pclass.dec)
  plot(pf.dec)
  pf.mix <- terra::focal(lc, w=matrix(1,9,9), pclass.mix)
  plot(pf.mix)

  ## LC layers
  plulc <- c(pf.nl, pf.dec, pf.mix)#, pf.wet, pf.grass, pf.shrub, pf.barren, pf.urban, pf.crop, pf.for, pf.open)
  names(plulc) <- c("Conifer", "Deciduous", "Mixed")#, "Wetland", "Grassland", "Shrubland", "Barren", "Urban", "Crop", "Forest", "Open")
  #pf.rs <- resample(plulc[[1]], bioclim[[1]], method = "bilinear")
  #plot(pf.rs)
  writeRaster(plulc, filename = "R:/Users/skeyser/PhD/NASA Project/Marten_Data/Pct_ForestTypes_30mMODIS_NALCMS_APIS.tif", overwrite = T)
} else {
  plulc <- rast("R:/Users/skeyser/PhD/NASA Project/Marten_Data/Pct_ForestTypes_30mMODIS_NALCMS_APIS.tif")
}

plulc3k <- resample(plulc, bioclim3k[[1]], method = "bilinear")

plulc <- project(plulc, crs(r.tmp))
plulc <- resample(plulc, r.tmp, method = "bilinear")
plot(plulc)

## LC-only
plulc <- project(lc, crs(r.tmp))
plulc <- resample(plulc, r.tmp, method = "near")
plot(plulc)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Productivity Data (DHIs)
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dhi.min <- rast("G:/DHI/MODIS_Aqua_NDVI_500m_2005_2020/MinDHI_scaled_MODIS_Aqua_NDVI_500m_2005-2020.tif")
dhi.cum <- rast("G:/DHI/MODIS_Aqua_NDVI_500m_2005_2020/CumDHI_scaled_MODIS_Aqua_NDVI_500m_2005-2020.tif")
dhi.var <- rast("G:/DHI/MODIS_Aqua_NDVI_500m_2005_2020/VarDHI_scaled_MODIS_Aqua_NDVI_500m_2005_2020.tif")
dhi.st <- c(dhi.cum, dhi.min, dhi.var)
dhi.st <- crop(dhi.st, dom.wgs84, mask = T)

dhi.st <- project(dhi.st, crs(r.tmp))
dhi.st <- resample(dhi.st, r.tmp, method = "bilinear")
names(dhi.st) <- c("DHI_Min", "DHI_Cum", "DHI_Var")

plot(dhi.st)


## Crop the DHIs to the correct extent
# bb <- terra::ext(chel[[1]])
# dhi.st <- crop(dhi.st, dom.wgs84, mask = T)
# plot(dhi.st[[1]])
#
# dhi.st <- resample(dhi.st, chel[[1]])

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Snow Water Equivalent
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

swe <- rast("E:/My Drive/GEE_Output/SWE_1980_2022_1km_CostLayerExt-0000000000-0000000000.tif")
swe <- crop(swe, dom.wgs84, mask = T)
swe <- project(swe, crs(r.tmp))
swe <- resample(swe, r.tmp, method = "bilinear")
names(swe) <- "SWE"

plot(swe)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Forest Structure ABG
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# abg <- rast("G:/GEDI_L4B_Gridded_Biomass_V2_1_2299/data/GEDI04_B_MW019MW223_02_002_02_R01000M_MU.tif")
# bbox <- st_read("R:/Users/skeyser/PhD/NASA Project/Marten_Data/CostLayerBBox.shp")
# bbox <- bbox |>
#   st_transform(crs = crs(abg)) |>
#   vect()
#
# abg.crop <- crop(abg, bbox)
# plot(abg.crop)

## -----------------------------------------------------------
##
## End Section: Rasters Loaded
##
## -----------------------------------------------------------

## Mask by SSL to remove lakes
env.st <- c(chel,
            plulc,
            can,
            dhi.st,
            ele,
            tpi,
            tr,
            hum,
            fh,
            #ssl.f,
            #scv.f,
            #fwos.f,
            ssl.f,
            scv.f,
            fwos.f,
            swe)

env.mask <- mask(env.st, ssl.f)
plot(env.mask)

## Write the rasterbrick for quick loading
writeRaster(env.st, filename = "R:/Users/skeyser/PhD/NASA Project/Marten_Data/EnvLayers_1km_MartenSuitability_APIS.tif", overwrite = T)
writeRaster(env.mask, filename = "R:/Users/skeyser/PhD/NASA Project/Marten_Data/EnvLayers_1km_MartenSuitability_APIS_Masked_Test.tif", overwrite = T)

## -----------------------------------------------------------
##
## Begin Section: Resample and stack all rasters together
##
## -----------------------------------------------------------

## Make a list of all the rasters
r.list <- list(chel,
               plulc,
               can,
               dhi.st,
               ele,
               tpi,
               tr,
               hum,
               fh,
               ssl.f,
               scv.f,
               fwos.f,
               swe)


rast.std.fun <- function(r.list, coord_ref = "EPSG:5070", dom.fun){
  ## We want everything to the 1km resolution so we will need to aggregate some layers
  ## Holder list
  std.list <- vector(mode = "list", length = length(r.list))
  std.list[[1]] <- r.list[[1]]

  ## Raster Template
  r.tmp <- rast(xmin = xmin(dom.fun),
                xmax = xmax(dom.fun),
                ymin = ymin(dom.fun),
                ymax = ymax(dom.fun),
                resolution = 1000,
                crs = crs(coord_ref),
                vals = NA)

  ## Loop for standardizing geospatial products
  ## Predicated on matching to bioclim specs cropped to the US
  for(i in 1:length(r.list)){
    print(paste("Processing:", names(r.list[[i]])))
    if(crs(r.tmp) != crs(r.list[[i]])){
      crs.tmp <- crs(r.list[[i]])
      dom.tmp <- project(dom.fun, crs.tmp)
      r.prj <- crop(r.list[[i]], dom.tmp, mask = TRUE)
      r.prj <- terra::project(r.prj, crs(r.tmp))
      r.prj <- terra::resample(r.prj, r.tmp, method = "bilinear")
      std.list[[i]] <- r.prj
    } else if (any(res(r.tmp) != res(r.list[[i]]))){
      r.crop <- crop(r.list[[i]], dom.fun, mask = TRUE)
      std.list[[i]] <- terra::resample(r.crop, r.tmp, method = "bilinear")
    } else {
      std.list[[i]] <- r.list[[i]]
    }}

}


env.st <- rast(std.list)
names(env.st) <- c("MATColdQ",
                   "MATWarmQ",
                   "MAT",
                   "PcpAnnual",
                   "Pct_Con",
                   "Pct_Dec",
                   "Pct_Mix",
                   "Pct_Canopy",
                   "DHI_Cum",
                   "DHI_Min",
                   "DHI_Var",
                   "Elevation",
                   "TPI",
                   "Roughness",
                   "gHM",
                   "ForHght",
                   "SSL",
                   "SCV",
                   "FWOS",
                   "SWE")
plot(env.st)

## Write the rasterbrick for quick loading
writeRaster(env.st, filename = "R:/Users/skeyser/PhD/NASA Project/Marten_Data/EnvLayers_1km_MartenSuitability.tif", overwrite = T)
