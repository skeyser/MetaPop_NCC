## ----------------------------------------------------------
##
## Script name: Marten RangeShiftR - RCM
##
## Script purpose: RangeShiftR simulations integrating
## climate change scenarios from RCMs
##
## Author: Spencer R Keyser
##
## Date Created: 2024-03-21
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
## Data manipulation
library(dplyr)
library(data.table)

## Directory management
library(here)

## Spatial
library(terra)
library(sf)

## Plotting
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(grid)
library(gridExtra)

## RS
library(RangeShiftR)

## -----------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: RS Directory Creation
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dirpath = here("Data/RS_Test/")

## -----------------------------------------------------------
##
## Begin Section: Landscape Parameters
##
## -----------------------------------------------------------

##------------------------------------------------------------
## Initialize the landscape
##------------------------------------------------------------
## ***********************************************************
##
## Section Notes:
## 1. LandscapeFile (file path) - (ascii format)
## 2. Resolution (integer) - Spatial resolution (Cell size in meters)
## 3. HabPercent (logical) - (F default) for discrete habitat classes
## and (T) for habitat percentage or suitability (expects 0-100 scale)
## 4. Nhabitats (integer) - # of habitats in the imported landscape
## 5. K_or_DensDep (vector or integer) - # of ind/ha
## If combined with StageStructured value provided will be used as
## strength of demographic density dependence (b^-1). If not,
## then value is interpreted as K. If HabPercent == F then
## K_or_DensDep is required for each landcover class otherwise
## value reaches a maximum in cells with 100% habitat and all
## other cells receive fractional support following K0rDensDep.
## 6. PatchFile (file path) - Path to patch file
## 7. CostsFile (file path) - Path to cost maos used for dispersal model
## 8. DynamicLandYears (Integer) - Years corresponding to switches in the
## habitat layers.
## 9. SpDistFile (file path) - Species Initial Distribution Map
## 10. SpDistResolution - Cell size in meters for the initial distribution
## Size can be coarser than the landscape
## ***********************************************************

## Load cost file and make matrix
cost <- rast(here("Data/Inputs/Cost_Water_1km_2.asc"))
plot(cost)
cost <- terra::as.matrix(cost, wide = TRUE)
cost[is.nan(cost)] <- NA
dim(cost)

## Replicate cost layers for the number of climate projections
cost <- list(cost, cost, cost)

## Demographic scaling layers
## Contemporary demographics
## Updated with adjusted winter length
adSurv <- rast(here("Data/Inputs/AdultSurvLayer_AllPatch_Current_WLadj.asc"))
cref <- crs(adSurv)
extAdSurv <- ext(adSurv)
origin_x <- extAdSurv[1]  # xmin
origin_y <- extAdSurv[3]  # ymin
adSurv <- adSurv * 100
plot(adSurv)

## Load and format the survival layers
surv.format.fun <- function(r.path){
  adSurv <- rast(r.path)
  cref <- crs(adSurv)
  extAdSurv <- ext(adSurv)
  origin_x <- extAdSurv[1]  # xmin
  origin_y <- extAdSurv[3]  # ymin
  adSurv <- adSurv * 100
  return(adSurv)
}

## Run the function on all scenarios
## 2 time points x 6 scenarios
## Modified for the WL adjustments
{
adSurv.MidAccess <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_MidAccess_WLadj.asc"))
adSurv.LateAccess <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_LateAccess_WLadj.asc"))
adSurv.MidCNRM <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_MidCNRM_WLadj.asc"))
adSurv.LateCNRM <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_LateCNRM_WLadj.asc"))
adSurv.MidGFDL <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_MidGFDL_WLadj.asc"))
adSurv.LateGFDL <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_LateGFDL_WLadj.asc"))
adSurv.MidIPSL <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_MidIPSL_WLadj.asc"))
adSurv.LateIPSL <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_LateIPSL_WLadj.asc"))
adSurv.MidMIROC5 <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_MidMIROC5_WLadj.asc"))
adSurv.LateMIROC5 <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_LateMIROC5_WLadj.asc"))
adSurv.MidMRI <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_MidMRI_WLadj.asc"))
adSurv.LateMRI <- surv.format.fun(r.path = here("Data/Inputs/AdultSurvLayer_AllPatch_LateMRI_WLadj.asc"))
}

## Copy over the adult survivals to the SA class
{
saSurv <- adSurv
saSurv.MidAccess <- adSurv.MidAccess
saSurv.LateAccess <- adSurv.LateAccess
saSurv.MidCNRM <- adSurv.MidCNRM
saSurv.LateCNRM <- adSurv.LateCNRM
saSurv.MidGFDL <- adSurv.MidGFDL
saSurv.LateGFDL <- adSurv.LateGFDL
saSurv.MidIPSL <- adSurv.MidIPSL
saSurv.LateIPSL <- adSurv.LateIPSL
saSurv.MidMIROC5 <- adSurv.MidMIROC5
saSurv.LateMIROC5 <- adSurv.LateMIROC5
saSurv.MidMRI <- adSurv.MidMRI
saSurv.LateMRI <- adSurv.LateMRI
}

## Landscape file
lscp <- app(adSurv, fun = function(x) ifelse(x > 0, 100, 0))
plot(lscp)
lscp <- terra::as.matrix(lscp, wide = T)

## Triple the list to ignore the errors
lscp <- list(lscp, lscp, lscp)

## Patch
patchShp <- st_read(here("Data/MartenPatchPolygons_NoAreaThresh.shp"))
patchShp <- patchShp |>
  mutate(Area = units::set_units(st_area(geometry), "km^2"))

## Load the Cost Data
template <- rast(here("Data/Inputs/Cost_Water_1km_2.asc"))
template[template >= 1] <- 0

## Rasterize
patches <- rasterize(vect(patchShp), template, field = "patches", background = 0)
patches <- crop(patches, template, mask = T)
plot(patches)

## Check patch indices with the output survival data
s.data <- read.csv(here("Data/Patch_Survival/all_patches_annual_survival_current_future_projections_2013_2022_wl_adjusted.csv"))
s.data.pid <- unique(s.data$patches)

## TRUE?
length(intersect(as.vector(unique(values(patches))), s.data.pid)) == length(s.data.pid)

p <- patches(adSurv, zeroAsNA = T)
p[is.na(p)] <- 0
plot(p)
p <- mask(p, adSurv)
plot(p)

patch <- terra::as.matrix(patches, wide = T)
patch[is.nan(patch)] <- NA
dim(patch)

## Same as landscape
patch <- list(patch, patch, patch)
str(patch)

## Function to prep layers
rs_prepDemog_fun <- function(adSurv,
                             adMidSurv,
                             adLateSurv,
                             saMidSurv,
                             saLateSurv){
  ndemo <- 1
  nc <- dim(adSurv)[2]
  nr <- dim(adSurv)[1]

  ## Kit survival layer (same as adult survival for mother dependency)
  demogKit <- array(terra::as.matrix(adSurv, wide = T), dim = c(nr, nc, ndemo))
  demogSA <- array(terra::as.matrix(saSurv, wide = T), dim = c(nr, nc, ndemo))
  demogAd <- array(terra::as.matrix(adSurv, wide = T), dim = c(nr, nc, ndemo))

  ## Combine: Kit, SA, Adult layers
  demog <- abind::abind(demogKit, demogSA, demogAd, along = 3)
  demog[is.nan(demog)] <- NA

  ## Mid-century
  demogKit.MidSurv <- array(terra::as.matrix(adMidSurv, wide = T), dim = c(nr, nc, ndemo))
  demogSA.MidSurv <- array(terra::as.matrix(saMidSurv, wide = T), dim = c(nr, nc, ndemo))
  demogAd.MidSurv <- array(terra::as.matrix(adMidSurv, wide = T), dim = c(nr, nc, ndemo))
  demog.MidSurv <- abind::abind(demogKit.MidSurv, demogSA.MidSurv, demogAd.MidSurv, along = 3)
  demog.MidSurv[is.nan(demog.MidSurv)] <- NA

  ## Late-century
  demogKit.LateSurv <- array(terra::as.matrix(adLateSurv, wide = T), dim = c(nr, nc, ndemo))
  demogSA.LateSurv <- array(terra::as.matrix(saLateSurv, wide = T), dim = c(nr, nc, ndemo))
  demogAd.LateSurv <- array(terra::as.matrix(adLateSurv, wide = T), dim = c(nr, nc, ndemo))
  demog.LateSurv <- abind::abind(demogKit.LateSurv, demogSA.LateSurv, demogAd.LateSurv, along = 3)
  demog.LateSurv[is.nan(demog.LateSurv)] <- NA

  demog <- list(demog, demog.MidSurv, demog.LateSurv)
  return(demog)
}

## Apply the function to all 6 scenarios
{
  demog_nocc <- rs_prepDemog_fun(adSurv = adSurv,
                                 adMidSurv = adSurv,
                                 saMidSurv = adSurv,
                                 adLateSurv = adSurv,
                                 saLateSurv = adSurv)

  demog_access <- rs_prepDemog_fun(adSurv = adSurv,
                                   adMidSurv = adSurv.MidAccess,
                                   saMidSurv = saSurv.MidAccess,
                                   adLateSurv = adSurv.LateAccess,
                                   saLateSurv = saSurv.LateAccess)

  demog_CNRM <- rs_prepDemog_fun(adSurv = adSurv,
                                 adMidSurv = adSurv.MidCNRM,
                                 saMidSurv = saSurv.MidCNRM,
                                 adLateSurv = adSurv.LateCNRM,
                                 saLateSurv = saSurv.LateCNRM)

  demog_GFDL <- rs_prepDemog_fun(adSurv = adSurv,
                                 adMidSurv = adSurv.MidGFDL,
                                 saMidSurv = saSurv.MidGFDL,
                                 adLateSurv = adSurv.LateGFDL,
                                 saLateSurv = saSurv.LateGFDL)

  demog_IPSL <- rs_prepDemog_fun(adSurv = adSurv,
                                 adMidSurv = adSurv.MidIPSL,
                                 saMidSurv = saSurv.MidIPSL,
                                 adLateSurv = adSurv.LateIPSL,
                                 saLateSurv = saSurv.LateIPSL)

  demog_MIROC5 <- rs_prepDemog_fun(adSurv = adSurv,
                                   adMidSurv = adSurv.MidMIROC5,
                                   saMidSurv = saSurv.MidMIROC5,
                                   adLateSurv = adSurv.LateMIROC5,
                                   saLateSurv = saSurv.LateMIROC5)

  demog_MRI <- rs_prepDemog_fun(adSurv = adSurv,
                                adMidSurv = adSurv.MidMRI,
                                saMidSurv = saSurv.MidMRI,
                                adLateSurv = adSurv.LateMRI,
                                saLateSurv = saSurv.LateMRI)
}

## Build the 6 Landscape params
ndemlyr <- 3L
k_orig <- 0.004 # Assumes 1 marten per ha
k <- c(1.2, 1.9, 2.4, 0.4, 0.43, 0.012, 0.042, 0.036) # Buskirk, Williams, Skalski, Clare
k_mean <- mean(k) / 100
{
  land_nocc <- ImportedLandscape(LandscapeMatrix = lscp,
                                   PatchMatrix = patch,
                                   CostsMatrix = cost,
                                   Resolution = 1000,
                                   HabPercent = T,
                                   K_or_DensDep = k_mean,
                                   demogScaleLayersMatrix = demog_nocc, #Just adult survival
                                   nrDemogScaleLayers = ndemlyr,
                                   DynamicLandYears = c(0L, 40L, 80L),
                                   OriginCoords = c(origin_x, origin_y)
  )

  land_access <- ImportedLandscape(LandscapeMatrix = lscp,
                                   PatchMatrix = patch,
                                   CostsMatrix = cost,
                                   Resolution = 1000,
                                   HabPercent = T,
                                   K_or_DensDep = k_mean,
                                   demogScaleLayersMatrix = demog_access, #Just adult survival
                                   nrDemogScaleLayers = ndemlyr,
                                   DynamicLandYears = c(0L, 40L, 80L),
                                   OriginCoords = c(origin_x, origin_y)
  )

  land_cnrm <- ImportedLandscape(LandscapeMatrix = lscp,
                                 PatchMatrix = patch,
                                 CostsMatrix = cost,
                                 Resolution = 1000,
                                 HabPercent = T,
                                 K_or_DensDep = k_mean,
                                 demogScaleLayersMatrix = demog_CNRM, #Just adult survival
                                 nrDemogScaleLayers = ndemlyr,
                                 DynamicLandYears = c(0L, 40L, 80L),
                                 OriginCoords = c(origin_x, origin_y)
  )
  land_gfdl <- ImportedLandscape(LandscapeMatrix = lscp,
                                 PatchMatrix = patch,
                                 CostsMatrix = cost,
                                 Resolution = 1000,
                                 HabPercent = T,
                                 K_or_DensDep = k_mean,
                                 demogScaleLayersMatrix = demog_GFDL, #Just adult survival
                                 nrDemogScaleLayers = ndemlyr,
                                 DynamicLandYears = c(0L, 40L, 80L),
                                 OriginCoords = c(origin_x, origin_y)
  )
  land_ipsl <- ImportedLandscape(LandscapeMatrix = lscp,
                                 PatchMatrix = patch,
                                 CostsMatrix = cost,
                                 Resolution = 1000,
                                 HabPercent = T,
                                 K_or_DensDep = k_mean,
                                 demogScaleLayersMatrix = demog_IPSL, #Just adult survival
                                 nrDemogScaleLayers = ndemlyr,
                                 DynamicLandYears = c(0L, 40L, 80L),
                                 OriginCoords = c(origin_x, origin_y)
  )
  land_miroc5 <- ImportedLandscape(LandscapeMatrix = lscp,
                                   PatchMatrix = patch,
                                   CostsMatrix = cost,
                                   Resolution = 1000,
                                   HabPercent = T,
                                   K_or_DensDep = k_mean,
                                   demogScaleLayersMatrix = demog_MIROC5, #Just adult survival
                                   nrDemogScaleLayers = ndemlyr,
                                   DynamicLandYears = c(0L, 40L, 80L),
                                   OriginCoords = c(origin_x, origin_y)
  )
  land_mri <- ImportedLandscape(LandscapeMatrix = lscp,
                                PatchMatrix = patch,
                                CostsMatrix = cost,
                                Resolution = 1000,
                                HabPercent = T,
                                K_or_DensDep = k_mean,
                                demogScaleLayersMatrix = demog_MRI, #Just adult survival
                                nrDemogScaleLayers = ndemlyr,
                                DynamicLandYears = c(0L, 40L, 80L),
                                OriginCoords = c(origin_x, origin_y)
  )
}

## -----------------------------------------------------------
##
## Scenario a: a sexual model with mate finding
##
## -----------------------------------------------------------

## Create the population transition matrixcrop
## Column 1 is juvenile, adult1, and adult2 stages
## This is similar to the three stage model
## we will use for marten
## Columns = Current Stage
## Rows = Effect of current stage on the next stage
## Stage classes: Kit (0-1), Juv (1-2), SA (2-3), Adult (3-9)
## Fecundity SA: 2.96/2
## Pregnancy rates: 56% SY and 79% 3+
## Fecundity accounting for pregnancy
## F SY = 0.56 * (2.96/2) (assuming Female only model 50/50 sex ratio)
## F +SY = 0.79 * (2.96/2)

## Juv survival is the hardest vital rate to get
## Wikston S1 = 0.71 (0.5 year post weaning) * Adult survival
## Manlick S1 = 0.39
## Average = 0.55
## SA Survival SNF = 0.81
## Adult Survival = Demographic Scaling Layer (1.00 max)

## Stages
## 1. Juv
## 2. Yearling
## 3. Subadult
## 4. Adult

nstage <- 4

## Trans Mat Wikston
(trans_mat2 <- matrix(c(0,0.71,0,0,
                       0,0,0.84,0,
                       0.83,0,0,1,
                       1.17,0,0,1
), nrow = nstage, byrow = F))


surv_mat_fonly <- c(1,NA,1,1)
#Changed from c(NA,NA,1,1)

## We can package this transition matrixs into a Stage Structured model obj
stg <- StageStructure(Stages = 4, # Same for marten potentially
                      TransMatrix = trans_mat2, # transition prob mat
                      MaxAge = 13, # max age in yrs
                      RepSeasons = 1,
                      RepInterval = 0,
                      SurvSched = 0, # when should survival and development occur
                      SurvDensDep = F,
                      FecDensDep = T, # make fecundity density dependent (1/b parameter)
                      SurvLayer = surv_mat_fonly)

demo <- Demography(StageStruct = stg,
                   ReproductionType = 0) # simple sexual model


getLocalisedEquilPop(demog = demo, DensDep_values = c(0.004, k_mean, 0.01))

## Add in dispersal components
## ***********************************************************
##
## Section Notes:
## Emigration()
## EmigProb takes the following structure in this example
## Column 1: Stage Classes (0,1,2) in three stage class model
## Marten Stage Classes - Kit -> Juvenile -> Adult
## Column 2: D0 = stage-specific max emigration probability
## Kit = 0, Juvenile (only dispersing stage class = 0.5), Adult = 0
## Column 3: alpha = stage-specific slope (how severely do individuals respond to density
## 1 equals a linear response, 10 is a logistic response)
## Column 4: beta = stage-specific inflection point (what is the scale of the response to relative
## population density (N/K or bN))
## These values are used for the following emigration probability (d)
## estimated as follows for (1) density dependent and (2) stage-strcutured
## (1) d(i,t) = D0 / (1 + exp(-asub>E (N(i,t) / K(i,t) - Betasub>E)))
## (2) d(i,t) = D0 / (1 + exp(-asub>E (b(i,t)*N(i,t) - Betasub>E)))
##
## SMS()
## 1. PR - Perceptual Range (# of cells)
## 2. PRMethod - Method to evaluate the effective cost for moving a step
##  1 = Arithmetic mean
##  2 = Harmonic mean
##  3 = Weighted mean
## 3. MemSize - Size of memory (# of steps over which to calculate the indivduals
## movement)
## 2. DP - Directional Persistence (tendency to follow a correlated random walk; must be >= 1
## think of this as an autocorrelated parameter)
##
## ***********************************************************

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Transfer
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Body size scaling for PR Mech and Zollner 2002
bm_f <- c(660, 694, 713, 600, 698, 637, 950, 1080, 860, 775, 901, 901, 889, 900, 817, 765, 550, 513, 550, 663, 672, 530, 430, 370, 890, 595)
bs <- mean(bm_f) #1kg for female marten
s <- 0.53839 - 0.00052*bs
ad <- c(63, 69)
fh <- 18
d <- (ad / s) * (fh/15.5)

## SMS dispersal
## Johnson et al., 2009
## x1 = 0, y1 = 0
## x2 = 45, y2 = 2.79
## Transform Hazard
hz2prob <- function(haz){
  prob <- 1 - exp(-haz)
  return(prob)
}

hazRate <- 2.79/45
stepMort <- hz2prob(haz = hazRate)

trans_sms <- SMS(PR = 1, # Perceptual Range (cells) 1-km
                 PRMethod = 2,
                 MemSize = 3,
                 DP = 5, # Directional persistence
                 Costs = "file", # Set for each habitat class in this example
                 StepMort = stepMort) # per step mortality probability (can be constant or habitat-specific)

trans_sms_nomort <- SMS(PR = 1, # Perceptual Range (cells) 1-km
                 PRMethod = 2,
                 MemSize = 3,
                 DP = 5, # Directional persistence
                 Costs = "file", # Set for each habitat class in this example
                 StepMort = 0) # per step mortality probability (can be constant or habitat-specific)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Emigration
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eprob <- cbind(0:3,
               c(1,0,0,0),
               c(10,0,0,0),
               c(1,0,0,0))

nodisp <- cbind(0:3,
                c(0,0,0,0),
                c(0,0,0,0),
                c(0,0,0,0))


emig_full <- Emigration(DensDep = T,
                        StageDep = T,
                        SexDep = F,
                        EmigProb = eprob
)

emig_none <- Emigration(DensDep = T,
                        StageDep = T,
                        SexDep = F,
                        EmigProb = nodisp
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Settlement
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
settle_sms <- Settlement(MaxSteps = 214,
                         StageDep = F,
                         SexDep = F,
                         Settle = 0,
                         DensDep = F
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Dispersal
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp_sms <- Dispersal() + emig_full + settle_sms + trans_sms
disp_none <- Dispersal() + emig_none + settle_sms + trans_sms
disp_nomort <- Dispersal() + emig_full + settle_sms + trans_sms_nomort

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Check parameters
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Look at our defined demographics and emigration probs
par(mfrow = c(2,2))
plotProbs(demo@StageStruct)
plotProbs(disp_sms@Emigration)
plotProbs(disp_none@Emigration)
plotProbs(disp_nomort@Emigration)

## -----------------------------------------------------------
##
## Begin Section: Initialisation
##
## -----------------------------------------------------------

## Now that we have defined the demographics and dispersal models
init_full_k <- Initialise(InitType = 0,
                         FreeType = 1,
                         InitDens = 0,
                         PropStages = c(0.0,0.30,0.2,0.5))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Simulation
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Spatial demo + SMS
RepNo <- 100


##Simulation = 666 for runs with WL adjustment
sim_ssms <- Simulation(Simulation = 666,
                       Replicates = RepNo,
                       Years = 100,
                       CreatePopFile = T,
                       ReturnPopMatrix = F,
                       ReturnPopDataFrame = T,
                       ReturnStages = T,
                       SMSHeatMap = T)

## Package simulation parameters up for the 6 scenarios
## No Climate change
s_nocc <- RSsim(batchnum = 000,
                  land = land_nocc,
                  demog = demo,
                  dispersal = disp_sms,
                  simul = sim_ssms,
                  init = init_full_k,
                  seed = 123)

## Access
s_access <- RSsim(batchnum = 111,
                land = land_access,
                demog = demo,
                dispersal = disp_sms,
                simul = sim_ssms,
                init = init_full_k,
                seed = 123)
##CNRM
s_cnrm <- RSsim(batchnum = 222,
                land = land_cnrm,
                demog = demo,
                dispersal = disp_sms,
                simul = sim_ssms,
                init = init_full_k,
                seed = 123)
## GFDL
s_gfdl <- RSsim(batchnum = 333,
                land = land_gfdl,
                demog = demo,
                dispersal = disp_sms,
                simul = sim_ssms,
                init = init_full_k,
                seed = 123)
## IPSL
s_ipsl <- RSsim(batchnum = 444,
                land = land_ipsl,
                demog = demo,
                dispersal = disp_sms,
                simul = sim_ssms,
                init = init_full_k,
                seed = 123)
## MIROC5
s_miroc5 <- RSsim(batchnum = 555,
                  land = land_miroc5,
                  demog = demo,
                  dispersal = disp_sms,
                  simul = sim_ssms,
                  init = init_full_k,
                  seed = 123)
## MRI
s_mri <- RSsim(batchnum = 666,
               land = land_mri,
               demog = demo,
               dispersal = disp_sms,
               simul = sim_ssms,
               init = init_full_k,
               seed = 123)

## Validate
all(c(
  validateRSparams(s_nocc),
  validateRSparams(s_access),
  validateRSparams(s_cnrm),
  validateRSparams(s_gfdl),
  validateRSparams(s_ipsl),
  validateRSparams(s_miroc5),
  validateRSparams(s_mri)
))

## Run the simulations
r_nocc <- RunRS(s_nocc, dirpath)
r_access <- RunRS(s_access, dirpath)
r_cnrm <- RunRS(s_cnrm, dirpath)
r_gfdl <- RunRS(s_gfdl, dirpath)
r_ipsl <- RunRS(s_ipsl, dirpath)
r_miroc5 <- RunRS(s_miroc5, dirpath)
r_mri <- RunRS(s_mri, dirpath)

## Keep only the simulation runs for the plotting script
keep.list <- c("s_nocc",
               "s_access",
               "s_cnrm",
               "s_gfdl",
               "s_ipsl",
               "s_miroc5",
               "s_mri")

rm(list = setdiff(ls(), keep.list))

## Save the simulation information for reference
#save.image(file = here("Data/RS_ParamMaster/RS_RCMs_Params_WLadj.RData"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Max dispersal no mortality
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Spatial demo + SMS
RepNo <- 100

##Simulation = 555 for prior runs without adjustment for WL
##Simulation = 666 for runs with WL adjustment

sim_ssms_nomort <- Simulation(Simulation = 6661,
                              Replicates = RepNo,
                              Years = 100,
                              CreatePopFile = T,
                              ReturnPopMatrix = F,
                              ReturnPopDataFrame = T,
                              ReturnStages = T,
                              SMSHeatMap = T)

## Package simulation parameters up for the 6 scenarios
## No CC
s_nocc_nomort <- RSsim(batchnum = 0001,
                         land = land_nocc,
                         demog = demo,
                         dispersal = disp_nomort,
                         simul = sim_ssms_nomort,
                         init = init_full_k,
                         seed = 123)
## Access
s_access_nomort <- RSsim(batchnum = 1111,
                  land = land_access,
                  demog = demo,
                  dispersal = disp_nomort,
                  simul = sim_ssms_nomort,
                  init = init_full_k,
                  seed = 123)
##CNRM
s_cnrm_nomort <- RSsim(batchnum = 2221,
                land = land_cnrm,
                demog = demo,
                dispersal = disp_nomort,
                simul = sim_ssms_nomort,
                init = init_full_k,
                seed = 123)
## GFDL
s_gfdl_nomort <- RSsim(batchnum = 3331,
                land = land_gfdl,
                demog = demo,
                dispersal = disp_nomort,
                simul = sim_ssms_nomort,
                init = init_full_k,
                seed = 123)
## IPSL
s_ipsl_nomort <- RSsim(batchnum = 4441,
                land = land_ipsl,
                demog = demo,
                dispersal = disp_nomort,
                simul = sim_ssms_nomort,
                init = init_full_k,
                seed = 123)
## MIROC5
s_miroc5_nomort <- RSsim(batchnum = 5551,
                  land = land_miroc5,
                  demog = demo,
                  dispersal = disp_nomort,
                  simul = sim_ssms_nomort,
                  init = init_full_k,
                  seed = 123)
## MRI
s_mri_nomort <- RSsim(batchnum = 6661,
               land = land_mri,
               demog = demo,
               dispersal = disp_nomort,
               simul = sim_ssms_nomort,
               init = init_full_k,
               seed = 123)

## Validate
all(c(
  validateRSparams(s_nocc_nomort),
  validateRSparams(s_access_nomort),
  validateRSparams(s_cnrm_nomort),
  validateRSparams(s_gfdl_nomort),
  validateRSparams(s_ipsl_nomort),
  validateRSparams(s_miroc5_nomort),
  validateRSparams(s_mri_nomort)
))

## Run the simulations
r_nocc_nomort <- RunRS(s_nocc_nomort, dirpath)
r_access_nomort <- RunRS(s_access_nomort, dirpath)
r_cnrm_nomort <- RunRS(s_cnrm_nomort, dirpath)
r_gfdl_nomort <- RunRS(s_gfdl_nomort, dirpath)
r_ipsl_nomort <- RunRS(s_ipsl_nomort, dirpath)
r_miroc5_nomort <- RunRS(s_miroc5_nomort, dirpath)
r_mri_nomort <- RunRS(s_mri_nomort, dirpath)

## Keep only the simulation runs for the plotting script
keep.list.nomort <- c(
  "s_nocc_nomort",
  "s_access_nomort",
  "s_cnrm_nomort",
  "s_gfdl_nomort",
  "s_ipsl_nomort",
  "s_miroc5_nomort",
  "s_mri_nomort")

rm(list = setdiff(ls(), keep.list.nomort))

## Save the simulation information for reference
#save.image(file = here("Data/RS_ParamMaster/RS_RCMs_Params_WLadj_NoMort.RData"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: No emigration
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Spatial demo + SMS
RepNo <- 100

##Simulation = 555 for prior runs without adjustment for WL
##Simulation = 666 for runs with WL adjustment

sim_ssms_noe <- Simulation(Simulation = 6662,
                           Replicates = RepNo,
                           Years = 100,
                           CreatePopFile = T,
                           ReturnPopMatrix = F,
                           ReturnPopDataFrame = T,
                           ReturnStages = T,
                           SMSHeatMap = T)

## Package simulation parameters up for the 6 scenarios
s_nocc_noe <- RSsim(batchnum = 0002,
                    land = land_nocc,
                    demog = demo,
                    dispersal = disp_none,
                    simul = sim_ssms_noe,
                    init = init_full_k,
                    seed = 123)
## Access
s_access_noe <- RSsim(batchnum = 1112,
                  land = land_access,
                  demog = demo,
                  dispersal = disp_none,
                  simul = sim_ssms_noe,
                  init = init_full_k,
                  seed = 123)
##CNRM
s_cnrm_noe <- RSsim(batchnum = 2222,
                land = land_cnrm,
                demog = demo,
                dispersal = disp_none,
                simul = sim_ssms_noe,
                init = init_full_k,
                seed = 123)
## GFDL
s_gfdl_noe <- RSsim(batchnum = 3332,
                land = land_gfdl,
                demog = demo,
                dispersal = disp_none,
                simul = sim_ssms_noe,
                init = init_full_k,
                seed = 123)
## IPSL
s_ipsl_noe <- RSsim(batchnum = 4442,
                land = land_ipsl,
                demog = demo,
                dispersal = disp_none,
                simul = sim_ssms_noe,
                init = init_full_k,
                seed = 123)
## MIROC5
s_miroc5_noe <- RSsim(batchnum = 5552,
                  land = land_miroc5,
                  demog = demo,
                  dispersal = disp_none,
                  simul = sim_ssms_noe,
                  init = init_full_k,
                  seed = 123)
## MRI
s_mri_noe <- RSsim(batchnum = 6662,
               land = land_mri,
               demog = demo,
               dispersal = disp_none,
               simul = sim_ssms_noe,
               init = init_full_k,
               seed = 123)

## Validate
all(c(
  validateRSparams(s_nocc_noe),
  validateRSparams(s_access_noe),
  validateRSparams(s_cnrm_noe),
  validateRSparams(s_gfdl_noe),
  validateRSparams(s_ipsl_noe),
  validateRSparams(s_miroc5_noe),
  validateRSparams(s_mri_noe)
))

## Run the simulations
r_nocc_noe <- RunRS(s_nocc_noe, dirpath)
r_access_noe <- RunRS(s_access_noe, dirpath)
r_cnrm_noe <- RunRS(s_cnrm_noe, dirpath)
r_gfdl_noe <- RunRS(s_gfdl_noe, dirpath)
r_ipsl_noe <- RunRS(s_ipsl_noe, dirpath)
r_miroc5_noe <- RunRS(s_miroc5_noe, dirpath)
r_mri_noe <- RunRS(s_mri_noe, dirpath)

## Keep only the simulation runs for the plotting script
keep.list <- c(
  "s_nocc_noe",
  "s_access_noe",
  "s_cnrm_noe",
  "s_gfdl_noe",
  "s_ipsl_noe",
  "s_miroc5_noe",
  "s_mri_noe")

rm(list = setdiff(ls(), keep.list))

## Save the simulation information for reference
#save.image(file = here("Data/RS_ParamMaster/RS_RCMs_Params_WLadj_NoEmig.RData"))
