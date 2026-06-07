## ----------------------------------------------------------
##
## Script name: Great Lakes Marten Habitat Suitability
##
## Script purpose: Model marten distribution/habitat suitability
## across the Great Lakes. Model will be used to delineate habitat
## patches for marten for IBM.
##
## Author: Spencer R Keyser
##
## Date Created: 2024-10-02
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

## Data manipulation/plotting
library(dplyr)
library(ggplot2)
library(data.table)

## Spatial
library(sf)
library(terra)
library(rasterVis)

## SDM
library(SDMtune)
library(flexsdm)
library(plotROC)
library(zeallot)
library(ecospat)


## -----------------------------------------------------------

## -----------------------------------------------------------
##
## Begin Section: Species Distribution Modelling
##
## -----------------------------------------------------------

## Cost Data for spatial reference
cost <- rast(here("./Data/Circuitscape/Input/predcosts.txt"))

## Shapefile for the US
us <- st_read(here("./Data/US_Map/states.shp"))
us <- us |>
  filter(STATE_NAME %in% c("Michigan", "Minnesota", "Wisconsin"))

## Marten Genotype data that serves as an occurrence dataset
marten_data <- data.table::fread(here("./Data/Individuals/dapc_dat.csv"))
marten_data <- as.data.frame(marten_data)
marten.sf <- marten_data |>
  rename("Lat" = "latitude-Location")  |>
  rename("Long" = "longitude-Location") |>
  select(Lat, Long, RecoveryArea) |>
  filter(!RecoveryArea %in% c("DouglasCo WI")) |>
  st_as_sf(coords = c("Long", "Lat"), crs = 4326, remove = F)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Environmental Predictor Data
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load in the data layers for marten at 1 km resolution
pred <- rast(here("./Data/EnvLayers_1km_MartenSuitability_APIS.tif"))

## Quick plotting
plot(pred)

## Load in ROI
us <- st_read(here("./Data/US_Map/states.shp"))
us <- us |>
  filter(STATE_NAME %in% c("Michigan", "Minnesota", "Wisconsin")) |>
  st_transform(crs = crs(pred))

## SpatVector
gl.v <- vect(us)

## Crop the predictors to the extent
pred.gl <- crop(pred, gl.v)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Assessing colinearity
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pearson_var <- flexsdm::correct_colinvar(pred, method = c("pearson", th = "0.7"))
pearson_var$cor_table
pearson_var$cor_variables

corrplot::corrplot(pearson_var$cor_table)

## Vars to keep
v.keep <- c("MAT", #add
            "MATColdQ",
            "Pct_Con",
            "Pct_Dec",
            "Pct_Mix",
            #"CoverClass", #add
            "PcpAnnual",
            "Pct_Can", #add
            "DHI_Min",
            "DHI_Cum", #add
            "DHI_Var", #add
            "Elevation",
            "TPI",
            "Roughness",
            "Human_Mod",
            "For_Hght",
            "SSL", #add
            "SCV", #add
            "FWOS",
            "SWE")

pred.mod <- pred[[v.keep]]
plot(pred.mod)

## Focalize to smooth the surfaces
pred.mod[["Pct_Can"]] <- pred.mod[["Pct_Can"]]/100

pearson_var <- flexsdm::correct_colinvar(pred.mod, method = c("pearson", th = "0.8"))
pearson_var$cor_table
pearson_var$cor_variables

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Predictor spatial projection
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Make the spatial predictors 1km
pred.mod <- project(pred.mod, "EPSG:5070")
ras.template <- rast(xmin = xmin(pred.mod),
                     xmax = xmax(pred.mod),
                     ymin = ymin(pred.mod),
                     ymax = ymax(pred.mod),
                     resolution = 1000)

pred.mod <- resample(pred.mod, ras.template)
plot(pred.mod)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Background sampling
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Extract the presence locations
marten.pts <- marten.sf |>
  st_transform(crs = "EPSG:5070") |>
  dplyr::mutate(Long.dd = sf::st_coordinates(geometry)[,1],
                Lat.dd = sf::st_coordinates(geometry)[,2]) |>
  select(Long = Long.dd, Lat = Lat.dd) |>
  st_drop_geometry() |>
  as.data.frame()

## Sample BG points
## Thickening methods from Vollering et al., 2018
bg.pts <- flexsdm::sample_background(marten.pts,
                                     x = "Long",
                                     y = "Lat",
                                     n = 10000,
                                     rlayer = pred.mod,
                                     method = "thickening")

## Only take the spatial coords
bg.pts <- bg.pts[,1:2]


## Prepare the SWD object
## Arbitrary species column
## Presences from thinned marten.pts
## Absences from thickened BG points
## Environmental data as 1km predictor stack
marten.dat <- prepareSWD(species = "Marten",
                         p = marten.pts,
                         a = bg.pts,
                         env = pred.mod)


## Add presences to the bg
marten.dat <- addSamplesToBg(marten.dat)

## Train/Test split
## Taking a hybrid 80/20 train/test
## 10-fold CV for the training data
## Withold the testing data for validation
c(train, test) %<-% trainValTest(marten.dat,
                                 test = 0.2,
                                 only_presence = T,
                                 seed = 54321)


## Distribution of variables in training data
str(train)
hist(train@data$MATColdQ)
hist(train@data$Pct_Con)
hist(train@data$Pct_Mix)
hist(train@data$Pct_Dec)
hist(train@data$SWE)

## Create folds for CV
## 10 folds
folds <- randomFolds(train,
                     k = 10,
                     only_presence = T,
                     seed = 25)

## Train model with the target species information
## Use the 10-folds for the training data
train.mod.default <- train(method = "Maxnet",
                           data = train,
                           folds = folds)

## Combine the CV
## Retrains the data using the full dataset and keeps model hyperparameters
train.mod.default <- combineCV(train.mod.default)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Hyperparameter Tuning
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Model Tuning parameters
## Vary feature combinations and reg
h <- list(reg = seq(0.25, 2, 0.25),
          fc = c("l",
                 "lq",
                 "lh",
                 "lp",
                 "lpq",
                 "lpqh"))


## Perform a random search of tuning parameters
## Tuning parameters:
## 1. Regularization multiplier
## 2. Feature constraints

## Random Search
tune_mods <- gridSearch(train.mod.default,
                        hypers = h,
                        metric = "auc",
                        test = test)

## Take the best model
bm.ind <- which.max(tune_mods@results$test_AUC)
final.model <- combineCV(tune_mods@models[[bm.ind]])

## Examine the PPM for the final model against the test dataset
## Default, untuned model validation
test.auc.default <- auc(train.mod.default, test = test)
test.tss.default <- tss(train.mod.default, test = test)

## Tuned model validation
test.auc <- auc(final.model, test = test)
test.tss <- tss(final.model, test = test)

## Best performing model validation
mart_sdm_val <- data.frame(Species = "Martes americana",
                           Reg = final.model@model@reg,
                           FC = final.model@model@fc,
                           AUC = test.auc,
                           TSS = test.tss)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Variable Importance and functional responses
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variable importance
vi <- SDMtune::varImp(bias.train)
vi <- SDMtune::varImp(final.model)

## Create a dataframe for VI with better labels
vi <- vi %>%
  mutate(PredPretty = case_when(Variable == "SSL" ~ "Snow Season Length",
                                Variable == "SCV" ~ "Snow Cover Variability",
                                Variable == "FWOS" ~ "Frozen Ground w/o Snow",
                                Variable == "For_Hght" ~ "Forest Height",
                                Variable == "MAT" ~ "Mean Annual T",
                                Variable == "MATColdQ" ~ "Min T Coldest Month",
                                Variable == "PcpAnnual" ~ "Annual Precipitation",
                                Variable == "Elevation" ~ "Mean Elevation",
                                Variable == "Roughness" ~ "Terrain Roughness",
                                Variable == "TPI" ~ "Topographic Position Index",
                                Variable == "Human_Mod" ~ "Human Modification",
                                Variable == "DHI_Min" ~ "Minimum Productivity",
                                Variable == "DHI_Cum" ~ "Cumulative Productivity",
                                Variable == "DHI_Var" ~ "Variability in Productivity",
                                Variable == "SWE" ~ "Snow Water Equivalent",
                                Variable == "Pct_Con" ~ "Percent Conifer",
                                Variable == "Pct_Dec" ~ "Percent Deciduous",
                                Variable == "Pct_Mix" ~ "Percent Mixed Forest",
                                Variable == "Pct_Can" ~ "Percent Canopy Cover",
                                Variable == "CoverClass" ~ "Landcover Class"
  ))


## Code for SI Figure
## Plot VI
vi_plot <- ggplot(data = vi, aes(x = reorder(PredPretty, Permutation_importance), y = Permutation_importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("Predictor") +
  ylab("Variable Importance") +
  theme_bw()

ggsave(plot = vi_plot,
       filename = here::here("Figures/VI_Plot_RegionalSDM_Marten_BestMod.jpg"),
       height = 8, width = 8,
       dpi = 600)

## -----------------------------------------------------------
##
## End Section: Modeling fitting and validation
##
## -----------------------------------------------------------

## -----------------------------------------------------------
##
## Begin Section: Model predictions
##
## -----------------------------------------------------------
## Bring in the cost layer for projections
## Use for constraining the prediction surface
cost <- rast(here("./Data/Inputs/Cost_Water_1km_2.asc"))
cost <- classify(cost, cbind(500, 1000, NA))
plot(cost)

## Crop the environmental prediction surface to Lydia's cost surface
env.pred <- crop(pred.mod, cost, mask = T)

## Focalize the predictors for continuous surface using means
env.pred <- focal(env.pred, w = 9, fun = "mean", na.policy = "only")
plot(env.pred)

## Predict across space with the fitted model
pred.map.exp <- predict(final.model,
                        data = env.pred,
                        type = "exponential")

## Predict across space with the transformed output
pred.map <- predict(final.model,
                    data = env.pred,
                    type = "cloglog")

## Mask to the cost layer
pred.map <- crop(pred.map, cost, mask = T)
plot(pred.map)

## Calculate Brier Score for Ben Murley
observed <- test@pa
coords <- st_as_sf(test@coords, coords = c("X", "Y"), crs = "EPSG:5070")
ggplot(coords) + geom_sf()
preds <- terra::extract(pred.map, vect(coords))
bs_df <- data.frame(Obs = observed, Pred = preds$lyr1)
bs_df <- bs_df[complete.cases(bs_df),]
brier <- mean((bs_df$Pred - bs_df$Obs)^2)

## Boyce-Index
## Get the test presence and all points
pred_coords <- test@coords
pres_coords <- test@coords[test@pa == 1, ]

## Extract the values from the predicted surface
all_suit <- terra::extract(pred.map, vect(st_as_sf(data.frame(pred_coords),
                                                           coords = c("X", "Y"),
                                                           crs = "EPSG:5070")),
                                   ID = FALSE)[,1]

pres_suit <- terra::extract(pred.map, vect(st_as_sf(data.frame(pres_coords),
                                                           coords = c("X", "Y"),
                                                           crs = "EPSG:5070")),
                                   ID = FALSE)[,1]

## Filter NAs from irregular prediction surface
all_suit <- all_suit[!is.na(all_suit)]
pres_suit <- pres_suit[!is.na(pres_suit)]

## Calculate the Boyce index
boyce_index <- ecospat.boyce(all_suit, pres_suit, nclass = 0,
                             window.w = "default", res = 100,
                             PEplot = TRUE,
                             method = "kendall")
boyce_index

## Plot creation for SI Figure 
## GGplot PE plot
pe <- ggplot() +
  geom_point(aes(x = boyce_index$HS,
                 y = boyce_index$F.ratio)) +
  theme_bw() +
  xlab("Predicted Suitability") +
  ylab("Predicted/Expecred Ratio")

ggsave(plot = pe, filename = here::here("Figures/BoyceIndexPlot.jpg"),
       height = 8, width = 8, dpi = 600)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Continuous mapping
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sdm.df <- as.data.frame(pred.map, xy = T)
colnames(sdm.df) <- c("x", "y", "occ")

sdm.map <- ggplot(data = sdm.df) +
  geom_raster(aes(x=x, y=y, fill=occ)) +
  geom_sf(data = us, fill = NA, size = 3, color = "black") +
  scale_fill_viridis_c(option = "G") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18,
                                    hjust = 0.5)) +
  guides(fill = guide_colorbar(title = "American Marten Probability of Occurrence",
                               title.position = "top",
                               barwidth = 20))

## Write the continuous predictions
writeRaster(pred.map, filename = here("Maxent_Output/Marten_BestMod_PredSurface.tif"))

marten.proj <- marten.sf |> st_transform(crs = "EPSG:5070")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Binary predictions
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Predictions using a 0.5 cutoff
binary.map <-
  sdm.df |>
  mutate(occ = as.factor(ifelse(occ >= 0.5, 1, 0))) |>
  ggplot() +
  geom_raster(aes(x=x, y=y, fill=occ)) +
  geom_sf(data = marten.proj, aes(color = RecoveryArea)) +
  scale_fill_viridis_d() +
  theme_void() +
  labs(fill = "Marten Patches")

# ggsave(plot = binary.map,
#        filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/SDM_Binary_1km_Marten_Genetics.jpg",
#        height = 8,
#        width = 12,
#        dpi = 600)

## Calculate thresholds for the binarizations
## cloglog
ths <- thresholds(final.model,
                  type = "cloglog")
## exponential
ths.exp <- thresholds(final.model,
                      type = "exponential")

## Equal sens. + spec. threshold
binary.map.ths <-
  sdm.df |>
  mutate(occ = as.factor(ifelse(occ >= ths[2,2], 1, 0))) |>
  ggplot() +
  geom_raster(aes(x=x, y=y, fill=occ)) +
  scale_fill_viridis_d() +
  theme_void() +
  labs(fill = "Marten Patches")


#writeRaster(bias.map, filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/MartenContinuous.asc")

## Produce and save the binary version of the SDM for later use
## Using threshold: Maximum training sensitivity plus specificity
plotPA(pred.map,
       th = ths[2,2],
       filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/All_Marten_Binary_Thresh_APIS.tif",
       overwrite = T
)

plotPA(pred.map,
       th = 0.5,
       filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/All_Marten_Binary_50_APIS.tif",
       overwrite = T
)

plotPA(pred.map,
       th = 0.75,
       filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/All_Marten_Binary_75_APIS.tif",
       overwrite = T
)

plotPA(pred.map,
       th = 0.9,
       filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/All_Marten_Binary_90_APIS.tif",
       overwrite = T
)


## -----------------------------------------------------------
##
## End Section: Model predictions
##
## -----------------------------------------------------------



## Bias Map patch IDs
patch <- rast("R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/All_Marten_Binary_Thresh_APIS.tif")
plot(patch)

## Make patches using 4 neighbors (no-diagonal)
patch.id <- terra::patches(patch, 4, zeroAsNA = T, allowGaps = FALSE)
plot(patch.id)

## Calculate Area for Each patch
z = cellSize(patch.id,unit="km") |>
  zonal(patch.id, sum, as.raster = T)

## Restrict by Carly's estimate of minimum patch size
## Start with 50 km2 as a reference for now
z.t <- ifel(z < 1, NA, patch.id)
plot(z.t)
table(z.t[!is.na(values(z.t))])

## Number of unique patches
length(unique(values(z.t)))

## Make this into a shapefile
z.v <- as.polygons(z.t, aggregate = T)
z.v <- st_as_sf(z.v)

## Write this to file
st_write(z.v, "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/MartenPatchPolygons_NoAreaThresh.shp")

ggplot(data = z.v) +
  geom_sf(data = us, fill = NA, color = "black", size = 2) +
  geom_sf(aes(fill = patches)) +
  scale_fill_viridis_c() +
  theme_bw()


mapview::mapView(z.t, zcol = "patches")

## Segregate the patches to individual layers
z.seg <- segregate(z.t)

## 22 distinct patches
nlyr(z.seg)
plot(z.seg)

## Take intial starting patch for the UP
up.patch <- z.seg["302"]

## For RangeShifter we need all non-patches that are matrix to be zero
background <- mask(patch, z.t, inverse = F, updatevalue = 0)
background <- rast(patch,
                   vals = 0)
background <- background + z.t
plot(background)

background <- classify(z.t, cbind(NA, 0))
plot(background)

background.up <- classify(up.patch, cbind(NA, 0))
plot(background.up)

## Mask by pixels with dispersal costs
background <- mask(background, cost)
plot(background)

## Mask up background
background.up <- mask(background.up, cost)

## Mask by the cost layer
patch.final <- crop(background, cost, mask = T)

patch <- mask(background, patch, updatevalue = 0)
plot(patch)

patch <- writeRaster(patch.rs,
                     filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/PatchID.asc",
                     NAflag = -9999)

## Pull in the genetic clusters
gen <- read.csv("R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/dapc_dat.csv")


patch.rs <- crop(patch, cost.rs)
plot(patch.rs)

patch.rs <- mask(patch, bias.map.rs)
plot(patch.rs)


## Look at Lydia's resistance surfacehttp://127.0.0.1:46769/graphics/plot_zoom_png?width=1920&height=1017
cost <- rast("R:/Users/skeyser/PhD/NASA Project/Spencer connectivity data/Circuitscape/Input/predcosts.txt")
plot(cost)

## Make the cost boundary a polygon
cost.rc <- classify(cost, cbind(0, Inf, 1))
plot(cost.rc)
cost.rc <- focal(cost.rc,
                 w = 15,
                 "max",
                 na.policy = "only",
                 fillvalue = NA)
plot(cost.rc)


cost.poly <- terra::as.polygons(cost.rc,
                                dissolve = T)
plot(cost.poly)

cost.poly.fill <- terra::fillHoles(cost.poly)
plot(cost.poly.fill)
values(cost.poly.fill) <- 1000

## Make new polygon back into raster
fill.ras <- terra::rasterize(cost.poly.fill, cost)

## Template raster matching Lydia's cost raster
c.temp <- rast(nrows = nrow(cost),
               ncols = ncol(cost),
               xmin = xmin(cost),
               xmax = xmax(cost),
               ymin = ymin(cost),
               ymax = ymax(cost),
               resolution = res(cost),
               vals = 1000)

## Update
cost.water <- classify(cost, cbind(NA, 1000))
plot(cost.water)

## Crop by the boundary of the image
cost.water <- mask(cost.water, cost.poly.fill)
plot(log(cost.water))

## Mask the patches by this cost layer
ext(patch) <- ext(cost.water)
cost.resample <- resample(cost.water, patch)
patch.mask <- mask(patch, cost.resample)

background.rs <- resample(background, patch.mask)
patch.id.mask <- mask(background.rs, patch.mask)
plot(patch.id.mask)

background.up <- resample(background.up, patch.mask)
up.init <- mask(background.up, patch.mask)
up.init <- classify()
plot(up.init)


## Patch ID
patch.id.mask <- patches(patch.mask, 4, allowGaps = F, zeroAsNA = F)
plot(patch.id.mask)

patch.up <- ifel(patch.id.mask != 302, 0, 1)
plot(patch.up)

terra::writeRaster(patch.id.mask,
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/Patch_ID_Test.asc",
                   NAflag = -9999,
                   overwrite = T)

terra::writeRaster(patch.mask,
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/GenoSuitable_Patch_Test.asc",
                   NAflag = -9999,
                   overwrite = T)

terra::writeRaster(patch.up,
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/UPInitMarten.asc",
                   NAflag = -9999,
                   overwrite = T)

terra::writeRaster(cost.water,
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/Cost_Water.asc",
                   NAflag = -9999,
                   overwrite = T)

terra::writeRaster(cost.resample,
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/Cost_Water_1km.asc",
                   NAflag = -9999,
                   overwrite = T)


## Change the extent of bias.map to costs
bias.map.crop <- crop(bias.map, cost)
plot(bias.map.crop)


## Match the resolution of the costs
cost.rs <- resample(cost, bias.map.crop)
plot(cost.rs)

## Bias Map
plot(bias.map.crop)
bias.map.rs <- mask(bias.map.crop, cost.rs)
plot(bias.map.rs)


## Percent Forest
pct_for <- pred.red[[3]] + pred.red[[4]] + pred.red[[5]]

## Crop to the extent
pct_for_crop <- crop(pct_for, cost.rs, mask = T)

## Scale ot 100 for RangeShifter
pct_for_crop <- pct_for_crop*100
plot(pct_for_crop)
table(is.na(values(pct_for_crop)))

## Forest Patches Greater than 80%
## Pct Forest Map
mat <- c(0, 80, 0,
         80, 100, 1)

rcl <- matrix(mat, ncol = 3, byrow = T)

for_patch <- classify(pct_for_crop, rcl)
plot(for_patch)

for_patch <- terra::patches(for_patch, 4, zeroAsNA = T)
plot(for_patch)

background <- rast(ymax = ymax(for_patch),
                   ymin = ymin(for_patch),
                   xmax = xmax(for_patch),
                   xmin = xmin(for_patch),
                   resolution = 4000,
                   crs = crs(for_patch),
                   vals = 0)
values(background) <- ifelse(is.na(values(for_patch)), NA, 0)
plot(background)


patch <- mask(for_patch, background)
plot(patch)

plot(patch)
plot(pct_for_crop)

table(is.na(values(patch)))
table(is.na(values(pct_for_crop)))

table(is.na(patch[pct_for_crop[]>0]))

patch.m <- mask(patch, pct_for_crop)

test <- c(bias.map.rs, cost.rs)
plot(test)

values(bias.map.rs) <- round(values(bias.map.rs), digits = 3)
values(cost.rs) <- signif(values(cost.rs), digits = 3)

terra::writeRaster(cost.rs,
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/Costs4km.asc",
                   NAflag = -9999,
                   overwrite = T)

raster::writeRaster(round(AllData2a*100), format="ascii", filename = "Inputs/climate_suitabilitya", NAflag = -9, overwrite = T, bylayer = T, datatype = "INT2U")

terra::writeRaster(round(bias.map.rs*100),
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/MartHabSuit.asc",
                   NAflag = -9999,
                   overwrite = T)


terra::writeRaster(patch,
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/Forest_Patch.asc",
                   NAflag = -9999,
                   overwrite = T)


terra::writeRaster(pct_for,
                   filename = "R:/Users/skeyser/PhD/NASA Project/R_projects/Marten_MetaPop/RangeShiftr_Inits/Pct_For.asc",
                   NAflag = -9999,
                   overwrite = T)
