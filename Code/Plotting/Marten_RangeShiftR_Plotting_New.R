## ----------------------------------------------------------
##
## Script name: RangeShifter RCM Plotting
##
## Script purpose: Plotting output from RCM runs for RangeShifter Models
##
## Author: Spencer R Keyser
##
## Date Created: 2025-07-02
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
library(purrr)
library(ggplot2)
library(RangeShiftR)
library(patchwork)
library(paletteer)
library(ggrepel)
library(sf)
library(terra)
library(tidyterra)
library(giscoR)
library(landscapemetrics)
library(rnaturalearth)
library(rnaturalearthdata)
## -----------------------------------------------------------

load(file = here("Data/RS_ParamMaster/RS_RCMs_Params_WLadj.RData"))
load(file = here("Data/RS_ParamMaster/RS_RCMs_Params_WLadj_NoMort.RData"))
load(file = here("Data/RS_ParamMaster/RS_RCMs_Params_WLadj_NoEmig.RData"))

## ***********************************************************
##
## Section Notes: Users need to specify the ABSOLUTE directory
## path here...RangeShiftR does not work with RELATIVE paths
## To run the script change the my_root string to be user
## specific.
##
## ***********************************************************
my_root <- "User/Specifc/Root/To/Downloaded/Rproj"

proj_dir <- "MetaPop_NCC/Data/RS_Test/"

dirpath = paste0(my_root, proj_dir)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Load the population simulation files
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## No CC
pop_nocc <- readPop(s_nocc, dirpath)
pop_nocc$Scenario <- "NoCC"
pop_nocc$Exp <- "Baseline"

## Access
pop_access <- readPop(s_access, dirpath)
pop_access$Scenario <- "Access"
pop_access$Exp <- "Baseline"

## CNRM
pop_cnrm <- readPop(s_cnrm, dirpath)
pop_cnrm$Scenario <- "CNRM"
pop_cnrm$Exp <- "Baseline"

## GFDL
pop_gfdl <- readPop(s_gfdl, dirpath)
pop_gfdl$Scenario <- "GFDL"
pop_gfdl$Exp <- "Baseline"

## IPSL
pop_ipsl <- readPop(s_ipsl, dirpath)
pop_ipsl$Scenario <- "IPSL"
pop_ipsl$Exp <- "Baseline"

## MIROC5
pop_miroc5 <- readPop(s_miroc5, dirpath)
pop_miroc5$Scenario <- "MIROC5"
pop_miroc5$Exp <- "Baseline"

## MRI
pop_mri <- readPop(s_mri, dirpath)
pop_mri$Scenario <- "MRI"
pop_mri$Exp <- "Baseline"

## Bind into a massive DF
pop_mart <- rbind(
  pop_nocc,
  pop_access,
  pop_cnrm,
  pop_gfdl,
  pop_ipsl,
  pop_miroc5,
  pop_mri)

## No Perstep Costs
## No CC
pop_nocc_nm <- readPop(s_nocc_nomort, dirpath)
pop_nocc_nm$Scenario <- "NoCC"
pop_nocc_nm$Exp <- "NoMort"

## Access
pop_access_nm <- readPop(s_access_nomort, dirpath)
pop_access_nm$Scenario <- "Access"
pop_access_nm$Exp <- "NoMort"

## CNRM
pop_cnrm_nm <- readPop(s_cnrm_nomort, dirpath)
pop_cnrm_nm$Scenario <- "CNRM"
pop_cnrm_nm$Exp <- "NoMort"

## GFDL
pop_gfdl_nm <- readPop(s_gfdl_nomort, dirpath)
pop_gfdl_nm$Scenario <- "GFDL"
pop_gfdl_nm$Exp <- "NoMort"

## IPSL
pop_ipsl_nm <- readPop(s_ipsl_nomort, dirpath)
pop_ipsl_nm$Scenario <- "IPSL"
pop_ipsl_nm$Exp <- "NoMort"

## MIROC5
pop_miroc5_nm <- readPop(s_miroc5_nomort, dirpath)
pop_miroc5_nm$Scenario <- "MIROC5"
pop_miroc5_nm$Exp <- "NoMort"

## MRI
pop_mri_nm <- readPop(s_mri_nomort, dirpath)
pop_mri_nm$Scenario <- "MRI"
pop_mri_nm$Exp <- "NoMort"

## Bind into a massive DF
pop_mart_nm <- rbind(
  pop_nocc_nm,
  pop_access_nm,
  pop_cnrm_nm,
  pop_gfdl_nm,
  pop_ipsl_nm,
  pop_miroc5_nm,
  pop_mri_nm)

## No Emigration
## No CC
pop_nocc_ne <- readPop(s_nocc_noe, dirpath)
pop_nocc_ne$Scenario <- "NoCC"
pop_nocc_ne$Exp <- "NoEmig"

## Access
pop_access_ne <- readPop(s_access_noe, dirpath)
pop_access_ne$Scenario <- "Access"
pop_access_ne$Exp <- "NoEmig"

## CNRM
pop_cnrm_ne <- readPop(s_cnrm_noe, dirpath)
pop_cnrm_ne$Scenario <- "CNRM"
pop_cnrm_ne$Exp <- "NoEmig"

## GFDL
pop_gfdl_ne <- readPop(s_gfdl_noe, dirpath)
pop_gfdl_ne$Scenario <- "GFDL"
pop_gfdl_ne$Exp <- "NoEmig"

## IPSL
pop_ipsl_ne <- readPop(s_ipsl_noe, dirpath)
pop_ipsl_ne$Scenario <- "IPSL"
pop_ipsl_ne$Exp <- "NoEmig"

## MIROC5
pop_miroc5_ne <- readPop(s_miroc5_noe, dirpath)
pop_miroc5_ne$Scenario <- "MIROC5"
pop_miroc5_ne$Exp <- "NoEmig"

## MRI
pop_mri_ne <- readPop(s_mri_noe, dirpath)
pop_mri_ne$Scenario <- "MRI"
pop_mri_ne$Exp <- "NoEmig"

## Bind into a massive DF
pop_mart_ne <- rbind(
  pop_nocc_ne,
  pop_access_ne,
  pop_cnrm_ne,
  pop_gfdl_ne,
  pop_ipsl_ne,
  pop_miroc5_ne,
  pop_mri_ne)

## Bind all of the variables together and clean
pop_mart <- rbind(pop_mart,
                  pop_mart_nm,
                  pop_mart_ne)

## Remove
rm(list = setdiff(ls(), "pop_mart"))
gc()
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Load in the spatial data for patches
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Bring in some state shapefiles for coarse assignment
gl.state <- st_read(here("./Data/US_Map/states.shp"))
gl.state <- gl.state |> filter(STATE_NAME %in% c("Michigan",
                                                 "Minnesota",
                                                 "Wisconsin"))

## Load the patch shapefile
patchShp <- st_read(here("./Data/MartenPatchPolygons_NoAreaThresh.shp"))

## Get area
patchShp <- patchShp |>
  mutate(Area = units::set_units(st_area(geometry), "km^2"))

patchShp$Lat <- st_coordinates(st_centroid(st_transform(patchShp, crs = 4326))$geometry)[,2]

## Function to calculate nearest neighbor distance
calc_nn_distance <- function(sf_obj) {
  nn_distances <- numeric(nrow(sf_obj))

  for(i in 1:nrow(sf_obj)) {
    # Calculate distances from patch i to all other patches
    distances <- st_distance(sf_obj[i,], sf_obj[-i,])
    # Get minimum distance
    nn_distances[i] <- min(distances)
  }

  return(as.numeric(nn_distances))
}

## Apply the function
patchShp <- st_transform(patchShp, crs = "EPSG:3174")
patchShp$nn_sf <- calc_nn_distance(patchShp)

## Plot
ggplot(data = patchShp) +
  geom_sf(aes(fill = nn_sf)) +
  scale_fill_viridis_c(name = "NN Distance (m)") +
  theme_void()

# Distance-weighted connectivity (Hanski-style)
calculate_connectivity_index <- function(patches, alpha = 0.001) {

  # Get centroids for distance calculations
  centroids <- st_centroid(patches)
  coords <- st_coordinates(centroids)

  # Calculate all pairwise distances
  dist_matrix <- as.matrix(dist(coords))
  area_vector <- st_area(patches)

  connectivity <- map_dbl(1:nrow(patches), function(i) {
    # For focal patch i, sum contributions from all other patches
    other_indices <- (1:nrow(patches))[-i]
    distances <- dist_matrix[i, other_indices]
    areas <- as.numeric(area_vector[other_indices])

    sum(areas * exp(-alpha * distances))
  })

  return(connectivity)
}

patchShp$con_ind <- calculate_connectivity_index(patchShp)

hist(con_ind)

## Transform the boundary layer
gl.state <- gl.state |> st_transform(crs = st_crs(patchShp))

# Spatial join
patches_in_states <- st_join(patchShp, gl.state, largest = TRUE)
patches_in_states <- patches_in_states |> arrange(patches)

# To summarize:
# Count patches per state
gl.state.pid <- patches_in_states |>
  group_by(STATE_NAME) |>
  summarize(
    patch_count = n(),
    patch_ids = list(unique(patches))
  )


## Merge popMart with spatial stuff
patch.geo <- patches_in_states |>
  dplyr::select(PatchID = patches,
                Area,
                Iso = nn_sf,
                Lat,
                Region = STATE_NAME) |>
  st_drop_geometry()

## Make the DF for plotting
pop_mart <- pop_mart |>
  left_join(patch.geo, by = "PatchID") |>
  mutate(Scenario = if_else(Scenario == "Access", "ACCESS", Scenario))
glimpse(pop_mart)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Abundance plots
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define colors for the 6 scenarios
no_scen <- 7
my_colors <- rev(c("#552000FF",
                   "#C17D17FF",
                   "#F8B150FF",
                   "#93C6E1FF",
                   "#5F93ACFF",
                   "#00344AFF",
                   "#C5C6C7"))

## Summarize the population trends
abun_sum <- pop_mart |>
  mutate(Scenario = factor(Scenario, levels = c(
    "NoCC",
    "GFDL",
    "MRI",
    "MIROC5",
    "CNRM",
    "ACCESS",
    "IPSL"))) |>
  mutate(Exp = factor(Exp, levels = c("NoEmig",
                                      "Baseline",
                                      "NoMort"),
                      labels = c("No Connectivity",
                                 "Limited Connectivity",
                                 "Max Connectivity"))) |>
  group_by(Year, Rep, Scenario, Exp) |>
  summarize(A = sum(NInd)) |>
  ungroup() |>
  group_by(Year, Scenario, Exp) |>
  summarise(Amean = mean(A),
            AUCI = quantile(A, probs = 0.975),
            ALCI = quantile(A, probs = 0.025),
            .groups = 'keep') |>
  ungroup()

# Filter for discrete time points (already in your code)
abun_pts <- abun_sum |>
  filter(Year %in% seq(0, 100, 20)) |>
  mutate(highlight = ifelse(Scenario %in% c("NoCC", "GFDL", "IPSL"), "highlight", "other"))

# Define dodge so points & error bars line up
pd <- position_dodge(width = 12)  # adjust dodge for spacing

## Fig 1A
Abun_Dots <- ggplot() +
  ## Intermediate scenarios
  geom_point(data = filter(abun_pts, highlight == "other"),
             aes(x = Year, y = Amean, color = Scenario),
             position = pd, alpha = 0.5, size = 2) +
  geom_errorbar(data = filter(abun_pts, highlight == "other"),
                aes(x = Year, y = Amean, color = Scenario,
                    ymin = ALCI, ymax = AUCI),
                position = pd, alpha = 0.5, width = 0) +
  geom_line(data = filter(abun_pts, highlight == "other"),
            aes(x = Year, y = Amean, color = Scenario, group = Scenario),
            position = pd, linewidth = 0.3, linetype = "dashed", alpha = 0.5) +
  ## Exemplars
  geom_line(data = filter(abun_pts, highlight == "highlight"),
            aes(x = Year, y = Amean, color = Scenario, group = Scenario),
            position = pd, linewidth = 0.5, linetype = "solid", alpha = 1) +
  geom_point(data = filter(abun_pts, highlight == "highlight"),
             aes(x = Year, y = Amean, color = Scenario),
             position = pd, size = 3) +
  geom_errorbar(data = filter(abun_pts, highlight == "highlight"),
                aes(x = Year, y = Amean, color = Scenario,
                    ymin = ALCI, ymax = AUCI),
                position = pd, width = 0) +
  scale_color_manual(values = my_colors,
                     breaks = c("NoCC", "GFDL", "MRI", "MIROC5", "CNRM", "ACCESS", "IPSL")) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme_bw() +
  labs(x = "Years from present",
       y = "Mean abundance (# ind.)",
       color = "GCM") +
  facet_wrap(~Exp) +
  theme(axis.title = element_text(size = 14, family = "sans"),
        axis.text = element_text(size = 12, family = "sans"),
        legend.text = element_text(size = 12, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, family = "sans", face = "bold"))

## Summary stats for change
abun_pts |>
  group_by(Scenario, Exp) |>
  summarise(PctChange = (last(Amean) - first(Amean))/ first(Amean),
            PctChangeAU = (min(AUCI) - max(AUCI))/ max(AUCI),
            PctChangeAL = (min(ALCI) - max(ALCI))/ max(ALCI)) |>
  print(n = 50)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Occupancy plots
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Occupancy plot
occ_sum <- pop_mart |>
  mutate(Scenario = factor(Scenario, levels = c(
    "NoCC",
    "GFDL",
                                                "MRI",
                                                "MIROC5",
                                                "CNRM",
                                                "ACCESS",
                                                "IPSL"))) |>
  mutate(Exp = factor(Exp, levels = c("NoEmig",
                                      "Baseline",
                                      "NoMort"),
                      labels = c("No Connectivity",
                                 "Limited Connectivity",
                                 "Max Connectivity"))) |>
  group_by(Year, Rep, Exp, Scenario) |>
  mutate(PatchOcc = ifelse(NInd > 0, 1, 0)) |>
  summarize(Nocc = sum(PatchOcc)) |>
  ungroup() |>
  group_by(Year, Scenario, Exp) |>
  summarise(Occmean = mean(Nocc),
            OccUCI = quantile(Nocc, probs = 0.975),
            OccLCI = quantile(Nocc, probs = 0.025),
            .groups = 'keep')

## Filter for discrete time points (already in your code)
occ_pts <- occ_sum |>
  filter(Year %in% seq(0, 100, 20)) |>
  mutate(highlight = ifelse(Scenario %in% c("NoCC", "GFDL", "IPSL"), "highlight", "other"))

occ_pts |>
  group_by(Scenario) |>
  summarise(PctChangeOcc = (min(Occmean) - max(Occmean))/ max(Occmean),
            PctChangeOccU = (min(OccUCI) - max(OccUCI))/ max(OccUCI),
            PctChangeOccL = (min(OccLCI) - max(OccLCI))/ max(OccLCI))

# Define dodge so points & error bars line up
pd <- position_dodge(width = 12)  # adjust dodge for spacing

# Final plot
Occ_Dots <- ggplot() +
  ## Intermediate scenarios
  geom_point(data = filter(occ_pts, highlight == "other"),
             aes(x = Year, y = Occmean, color = Scenario),
             position = pd, alpha = 0.5, size = 2) +
  geom_errorbar(data = filter(occ_pts, highlight == "other"),
                aes(x = Year, y = Occmean, color = Scenario,
                    ymin = OccLCI, ymax = OccUCI),
                position = pd, alpha = 0.5, width = 0) +
  geom_line(data = filter(occ_pts, highlight == "other"),
            aes(x = Year, y = Occmean, color = Scenario, group = Scenario),
            position = pd, linewidth = 0.3, linetype = "dashed", alpha = 0.5) +
  ## Exemplars
  geom_line(data = filter(occ_pts, highlight == "highlight"),
            aes(x = Year, y = Occmean, color = Scenario, group = Scenario),
            position = pd, linewidth = 0.5, linetype = "solid", alpha = 1) +
  geom_point(data = filter(occ_pts, highlight == "highlight"),
             aes(x = Year, y = Occmean, color = Scenario),
             position = pd, size = 3) +
  geom_errorbar(data = filter(occ_pts, highlight == "highlight"),
                aes(x = Year, y = Occmean, color = Scenario,
                    ymin = OccLCI, ymax = OccUCI),
                position = pd, width = 0) +
  scale_color_manual(values = my_colors,
                     breaks = c("NoCC", "GFDL", "MRI", "MIROC5", "CNRM", "ACCESS", "IPSL")) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme_bw() +
  labs(x = "Years from present",
       y = "Mean # of occupied patches",
       color = "GCM") +
  theme(axis.title = element_text(size = 14, family = "sans"),
        axis.text = element_text(size = 12, family = "sans"),
        legend.text = element_text(size = 12, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  facet_wrap(~Exp)


mart_demo_patch <- Abun_Dots / Occ_Dots + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(plot = mart_demo_patch,
       filename = here("./Figures/Fig2.jpg"),
       height = 8, width = 10, dpi = 800)

## -----------------------------------------------------------
##
## Begin Section: Summary statistics for pop declines
##
## -----------------------------------------------------------

glimpse(abun_sum)

abun_sum |>
  group_by(Scenario, Exp) |>
  mutate(Total_Diff = )

# Compare each climate scenario to NoCC at Year 100
baseline <- abun_sum %>%
  filter(Scenario == "NoCC", Year %in% c(60, 100)) %>%
  group_by(Exp, Year) %>%
  summarise(Baseline_Abundance = Amean,
            Baseline_Abundance_U = AUCI,
            Baseline_Abundance_L = ALCI,
            .groups = "drop")

tot_change <- abun_sum %>%
  filter(Year %in% c(60, 100) & Scenario != "NoCC") %>%
  left_join(baseline, by = c("Exp", "Year")) %>%
  mutate(
    Percent_Change_vs_NoCC = ((Amean - Baseline_Abundance) / Baseline_Abundance) * 100,
    Percent_Change_vs_NoCC_U = ((AUCI - Baseline_Abundance_U) / Baseline_Abundance_U) * 100,
    Percent_Change_vs_NoCC_L = ((ALCI - Baseline_Abundance_L) / Baseline_Abundance_L) * 100
  )

print(tot_change, n = 50)

baseline.o <- occ_sum %>%
  filter(Scenario == "NoCC", Year %in% c(20, 60, 100)) %>%
  group_by(Exp, Year) %>%
  summarise(Baseline_Occ = Occmean, .groups = "drop")

tot_change_occ <- occ_sum %>%
  filter(Year %in% c(60, 100) & Scenario != "NoCC") %>%
  left_join(baseline.o, by = c("Exp", "Year")) %>%
  mutate(
    Percent_Change_vs_NoCC = ((Occmean - Baseline_Occ) / Baseline_Occ) * 100
  )

print(tot_change_occ, n = 50)

## Dispesal scenarios
baseline_within_gcm <- abun_sum %>%
  filter(Exp == "No Connectivity", Year %in% c(100)) %>%
  group_by(Scenario, Year) %>%
  summarise(Baseline_NoConn = Amean,
            Baseline_NoConn_U = AUCI,
            Baseline_NoConn_L = ALCI,
            .groups = "drop")

# Compare Limited and Max Connectivity to No Connectivity within each GCM
within_gcm_change <- abun_sum %>%
  filter(Year %in% c(100) & Exp != "No Connectivity") %>%
  left_join(baseline_within_gcm, by = c("Scenario", "Year")) %>%
  mutate(
    Percent_Change_vs_NoConn = ((Amean - Baseline_NoConn) / Baseline_NoConn) * 100,
    Percent_Change_vs_NoConn_U = ((AUCI - Baseline_NoConn_U) / Baseline_NoConn_U) * 100,
    Percent_Change_vs_NoConn_L = ((ALCI - Baseline_NoConn_L) / Baseline_NoConn_L) * 100
  )

# Same for occupancy
baseline_within_gcm_occ <- occ_sum %>%
  filter(Exp == "No Connectivity", Year %in% c(100)) %>%
  group_by(Scenario, Year) %>%
  summarise(Baseline_NoConn_Occ = Occmean, .groups = "drop")

within_gcm_change_occ <- occ_sum %>%
  filter(Year %in% c(100) & Exp != "No Connectivity") %>%
  left_join(baseline_within_gcm_occ, by = c("Scenario", "Year")) %>%
  mutate(
    Percent_Change_vs_NoConn = ((Occmean - Baseline_NoConn_Occ) / Baseline_NoConn_Occ) * 100
  )

# Key distillation metrics
summary_stats <- within_gcm_change_occ %>%
  filter(Year == 100) %>%
  group_by(Scenario) %>%
  summarise(
    Max_Connectivity_Benefit = max(Percent_Change_vs_NoConn),
    Mean_Connectivity_Benefit = mean(Percent_Change_vs_NoConn),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_Connectivity_Benefit))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Patch Extinction Probability trends
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RepNo <- 100
Calc_PatchExtProb <- function(pop_df, RepNb = RepNo) {
  patch_extprob <- pop_df |>
    group_by(Rep, Year, PatchID, Scenario, Exp) |>
    summarise(sumPop = sum(NInd), .groups='keep') %>%
    group_by(Year, PatchID, Scenario, Exp) %>%
    # Average extinction probability (1 minus the proportion of replicates with surviving populations)
    summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb) %>%
    # Make sure that data frame is filled until last year of simulation
    #right_join(tibble(Year = seq_len(s@simul@Years)), by='Year') %>%
    mutate(extProb = ifelse(is.na(extProb), 1, extProb))

  return(patch_extprob)
}

## Calculate
patch_extprob <- Calc_PatchExtProb(pop_mart)
glimpse(patch_extprob)

patch_extprob <- patch_extprob |> left_join(patch.geo, by = "PatchID")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Extinction Probability Maps
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# States
states <- ne_states(country = "united states of america", returnclass = "sf") %>%
  filter(name %in% c("Michigan", "Wisconsin", "Minnesota")) %>%
  st_transform(crs = st_crs(patchShp))%>%
  st_crop(st_bbox(patchShp))

prov <- ne_states(country = "canada", returnclass = "sf") %>%
  filter(name %in% c("Ontario")) %>%
  st_transform(crs = st_crs(patchShp))%>%
  st_crop(st_bbox(patchShp))

# Great Lakes
lakes <- ne_download(scale = 50, type = "lakes", category = "physical", returnclass = "sf") %>%
  st_transform(crs = st_crs(patchShp)) %>%
  st_crop(st_bbox(patchShp))

## Calculate extinction probabilities first, then join spatial data
ext_summary <- patch_extprob |>
  filter(PatchID != 0 & !is.na(Region)) |>
  filter(Scenario %in% c("NoCC", "GFDL", "IPSL") & Year %in% c(100)) |>
  group_by(Scenario, PatchID, Year, Exp) |>
  summarise(
    median_prob = median(extProb),
    .groups = "drop"
  ) |>
  mutate(
    Exp = factor(Exp,
                 levels = c("NoEmig", "Baseline", "NoMort"),
                 labels = c("No\nConnectivity", "Limited\nConnectivity", "Max\nConnectivity"))
  ) |>
  mutate(
    Scenario = factor(Scenario,
                 levels = c("NoCC", "GFDL", "IPSL"),
                 labels = c("No Change", "GFDL", "IPSL"))
  )


ext.map <- ext_summary |>
  left_join(patchShp, by = c("PatchID" = "patches")) |>
  filter(Scenario %in% c("GFDL", "IPSL") & Exp %in% c("No\nConnectivity", "Max\nConnectivity")) |>
  st_as_sf() |>
  ggplot() +
  # Basemap layers
  geom_sf(data = states, fill = "grey98", color = "grey80", size = 0.3) +
  geom_sf(data = prov, fill = "gray98", color = "grey80", size = 0.3) +
  geom_sf(data = lakes, fill = "aliceblue", color = "lightsteelblue", size = 0.2) +
  # Your data - no borders around patches
  geom_sf(aes(fill = median_prob), color = NA, size = 0.3) +
  scale_fill_gradient2(name = "Extinction Probability",
                       #low = "#023198", mid = "#C0E9C2", high = "#7E1700",
                       low = "#003f5c",      # Dark blue (safe)
                       mid = "#ffa600",      # Orange (moderate risk)
                       high = "#d62728",
                       midpoint = 0.5,
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, margin = margin(5, 5, 5, 5)),
    legend.position = "bottom",
      plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
    strip.text.y = element_text(angle = 270)
  ) +
  guides(fill = guide_colorbar(barwidth = 10)) +
  facet_grid(vars(Exp), vars(Scenario))

ggsave(plot = ext.map, filename = here("./Figures/ExtProb_Year100_GFDL_IPSL_2x2.jpg"),
       height = 10, width = 8, dpi = 800)


## Extinction Probability maps for remaining scenarios for SI
ext_summary_si <- patch_extprob |>
  filter(PatchID != 0 & !is.na(Region)) |>
  filter(Year %in% c(100)) |>
  group_by(Scenario, PatchID, Year, Exp) |>
  summarise(
    median_prob = median(extProb),
    .groups = "drop"
  ) |>
  mutate(
    Exp = factor(Exp,
                 levels = c("NoEmig", "Baseline", "NoMort"),
                 labels = c("No\nConnectivity", "Limited\nConnectivity", "Max\nConnectivity"))
  ) |>
  mutate(
    Scenario = factor(Scenario,
                      levels = c("NoCC", "GFDL", "MRI", "MIROC5", "CNRM", "ACCESS", "IPSL"),
                      labels = c("No Change", "GFDL", "MRI", "MIROC5", "CNRM", "ACCESS", "IPSL"))
  )

ext.map.si <- ext_summary_si |>
  left_join(patchShp, by = c("PatchID" = "patches")) |>
  st_as_sf() |>
  ggplot() +
  # Basemap layers
  geom_sf(data = states, fill = "grey98", color = "grey80", size = 0.3) +
  geom_sf(data = prov, fill = "gray98", color = "grey80", size = 0.3) +
  geom_sf(data = lakes, fill = "aliceblue", color = "lightsteelblue", size = 0.2) +
  # Your data - no borders around patches
  geom_sf(aes(fill = median_prob), color = NA, size = 0.3) +
  scale_fill_gradient2(name = "Extinction Probability",
                       #low = "#023198", mid = "#C0E9C2", high = "#7E1700",
                       low = "#003f5c",      # Dark blue (safe)
                       mid = "#ffa600",      # Orange (moderate risk)
                       high = "#d62728",
                       midpoint = 0.5,
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, margin = margin(5, 5, 5, 5)),
    legend.position = "bottom",
    plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
    strip.text.y = element_text(angle = 270)
  ) +
  guides(fill = guide_colorbar(barwidth = 10)) +
  facet_grid(vars(Scenario), vars(Exp))

ggsave(plot = ext.map.si, filename = here("./Figures/ExtProb_Year100_Remaining_Scens_SI.jpg"),
       height = 8, width = 8, dpi = 800)


## -----------------------------------------------------------
##
## End Section: Extinction Probability Maps
##
## -----------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Abundance Maps
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Calculate extinction probabilities first, then join spatial data
abun_map <- pop_mart |>
  mutate(Scenario = factor(Scenario, levels = c(
    "NoCC",
    "GFDL",
    "MRI",
    "MIROC5",
    "CNRM",
    "Access",
    "IPSL"))) |>
  mutate(Exp = factor(Exp, levels = c("NoEmig",
                                      "Baseline",
                                      "NoMort"),
                      labels = c("No\nConnectivity",
                                 "Limited\nConnectivity",
                                 "Max\nConnectivity"))) |>
  filter(Scenario %in% c("NoCC", "GFDL", "IPSL") & Year == 100) |>
  group_by(PatchID, Scenario, Exp) |>
  summarize(Amean = mean(NInd)) |>
  ungroup()

a.map <- abun_map |>
  left_join(patchShp, by = c("PatchID" = "patches")) |>
  mutate(
    Scenario = factor(Scenario,
                      levels = c("NoCC", "GFDL", "IPSL"),
                      labels = c("No Change", "GFDL", "IPSL"))
  ) |>
  st_as_sf() |>
  ggplot() +
  # Basemap layers
  geom_sf(data = states, fill = "grey98", color = "grey80", size = 0.3) +
  geom_sf(data = prov, fill = "gray98", color = "grey80", size = 0.3) +
  geom_sf(data = lakes, fill = "aliceblue", color = "lightsteelblue", size = 0.2) +
  # Your data - no borders around patches
  geom_sf(aes(fill = Amean), color = NA, size = 0.3) +
  scale_fill_gradientn(name = "Abundance (mean # ind.)",
                       colors = c("#d73027", "#fee08b", "#4575b4"),
                       values = scales::rescale(c(50, 500, max(abun_map$Amean))),
                       breaks = c(50, 500, 2000),
                       labels = c("50", "500", "2000")
                       ) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, margin = margin(5, 5, 5, 5)),
    legend.position = "bottom",
    plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
    strip.text.y = element_text(angle = 270)
  ) +
  guides(fill = guide_colorbar(barwidth = 10)) +
  facet_grid(vars(Exp), vars(Scenario))

ggsave(plot = a.map, filename = here("./Figures/AbunMap_Year100_GFDL_IPSL_NoCC.jpg"),
       height = 8, width = 8, dpi = 800)


## Abundance Map for the Remaining Scenarios for SI
abun_map_si <- pop_mart |>
  mutate(Scenario = factor(Scenario, levels = c(
    "NoCC",
    "GFDL",
    "MRI",
    "MIROC5",
    "CNRM",
    "ACCESS",
    "IPSL"))) |>
  mutate(Exp = factor(Exp, levels = c("NoEmig",
                                      "Baseline",
                                      "NoMort"),
                      labels = c("No\nConnectivity",
                                 "Limited\nConnectivity",
                                 "Max\nConnectivity"))) |>
  filter(Scenario %in% c("MRI", "MIROC5", "CNRM", "ACCESS") & Year == 100) |>
  group_by(PatchID, Scenario, Exp) |>
  summarize(Amean = mean(NInd)) |>
  ungroup()

a.map.si <- abun_map_si |>
  left_join(patchShp, by = c("PatchID" = "patches")) |>
  mutate(
    Scenario = factor(Scenario,
                      levels = c("MRI", "MIROC5", "CNRM", "ACCESS"),
                      labels = c("MRI", "MIROC5", "CNRM", "ACCESS"))
  ) |>
  st_as_sf() |>
  ggplot() +
  # Basemap layers
  geom_sf(data = states, fill = "grey98", color = "grey80", size = 0.3) +
  geom_sf(data = prov, fill = "gray98", color = "grey80", size = 0.3) +
  geom_sf(data = lakes, fill = "aliceblue", color = "lightsteelblue", size = 0.2) +
  # Your data - no borders around patches
  geom_sf(aes(fill = Amean), color = NA, size = 0.3) +
  scale_fill_gradientn(name = "Abundance (mean # ind.)",
                       colors = c("#d73027", "#fee08b", "#4575b4"),
                       values = scales::rescale(c(50, 500, max(abun_map_si$Amean))),
                       breaks = c(50, 500, 2000),
                       labels = c("50", "500", "2000")
  ) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, margin = margin(5, 5, 5, 5)),
    legend.position = "bottom",
    plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
    strip.text.y = element_text(angle = 270)
  ) +
  guides(fill = guide_colorbar(barwidth = 10)) +
  facet_grid(vars(Scenario), vars(Exp))

ggsave(plot = a.map.si, filename = here("./Figures/AbunMap_Year100_Other_Scens_SI.jpg"),
       height = 8, width = 8, dpi = 800)


## Abundance Maps

sdl_p_sims_areawt <- patch_extprob |>
  filter(PatchID != 0 & !is.na(Region)) |>
  mutate(Scenario = factor(Scenario, levels = c(
    "NoCC",
    "GFDL",
                                                "MRI",
                                                "MIROC5",
                                                "CNRM",
                                                "ACCESS",
                                                "IPSL"))) |>
  group_by(Region, Year, PatchID, Area, Scenario, Exp) |>
  mutate(Area = as.numeric(Area)) |>
  summarise(
    median_prob = median(extProb),
    .groups = "drop"
  ) |>
  ggplot(aes(x=Year)) +
  geom_line(aes(y=median_prob, color=log(Area), group=PatchID),
            alpha=0.75) +
  geom_smooth(data = . %>%
                group_by(Region, Year, Scenario) %>%
                summarise(weighted_mean = weighted.mean(median_prob, w = Area)),
              aes(y=weighted_mean),
              color="black", size=1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_color_viridis_c(name = "ln(Patch Area)", option = "rocket") +
  #scale_color_distiller(name = "Log Patch Area") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  theme_bw() +
  ylab("Extinction Probability") +
  xlab("Years from present") +
  facet_grid(vars(Exp), vars(Scenario)) +
  theme(axis.title = element_text(size = 14, family = "sans"),
        axis.text = element_text(size = 12, family = "sans"),
        legend.text = element_text(size = 12, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        strip.text = element_text(size = 14, family = "sans", face = "bold"),
        strip.background = element_blank())

sdl_p_sims_areawt

## Heat Map
sdl_p_sims_heatmap <- patch_extprob %>%
  filter(PatchID != 0 & !is.na(Region)) %>%
  mutate(Scenario = factor(Scenario, levels = c("GFDL",
                                                "MRI",
                                                "MIROC5",
                                                "CNRM",
                                                "ACCESS",
                                                "IPSL"))) %>%
  mutate(Exp = factor(Exp, levels = c("NoEmig",
                                      "Baseline",
                                      "NoMort"),
                      labels = c("No Emigration",
                                 "Dispersal: Mortality",
                                 "Dispersal: No Mortality"))) |>
  group_by(Year, PatchID, Scenario, Exp) %>%
  summarise(
    median_prob = median(extProb),
    .groups = "drop"
  ) %>%
  # Bin into 10-year intervals and 0.1 extinction probability intervals:
  mutate(
    Year_bin = floor(Year / 10) * 10,                 # e.g., 0-9 → 0, 10-19 → 10, etc.
    Ext_bin  = cut(median_prob,
                   breaks = seq(0, 1, by = 0.1),
                   include.lowest = TRUE,
                   right = FALSE)
  ) %>%
  group_by(Scenario, Exp, Year_bin, Ext_bin) %>%
  summarise(n_patches = n_distinct(PatchID), .groups = "drop") %>%
  # Convert Ext_bin to numeric bin centers for plotting:
  mutate(
    Ext_bin_num = as.numeric(sub("\\[|\\)", "", sub(",.*", "", Ext_bin))) + 0.05
  ) |>
  tidyr::complete(
    Scenario,
    Exp,
    Year_bin = seq(0, 100, by = 10),                # your x bins
    Ext_bin_num = seq(0.05, 0.95, by = 0.1),         # y bin midpoints
    fill = list(n_patches = 0)                      # fill missing with 0
  )

# Plot heat map
patch_heatmap <- ggplot(sdl_p_sims_heatmap, aes(x = Year_bin, y = Ext_bin_num, fill = n_patches)) +
  geom_tile(na.rm = T) +
  scale_fill_viridis_c(name = "# Patches", option = "mako") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(
    x = "Years from present",
    y = "Extinction Probability"
  ) +
  facet_grid(vars(Exp), vars(Scenario)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, family = "sans"),
    axis.text = element_text(size = 12, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 14, family = "sans"),
    strip.text = element_text(size = 14, family = "sans", face = "bold"),
    strip.background = element_blank()
  )

ggsave(plot = patch_heatmap, filename = "R:/Users/skeyser/PhD/NASA Project/Writing/Marten_Chapter/NCC_Sub/MetaPatchHeatMap.jpg",
       height = 8, width = 10, dpi = 800)

sdl_p_sims_heatmap_prop <- patch_extprob %>%
  filter(PatchID != 0 & !is.na(Region)) %>%
  mutate(
    Scenario = factor(Scenario, levels = c("GFDL","MRI","MIROC5","CNRM","ACCESS","IPSL")),
    Exp = factor(Exp, levels = c("NoEmig","Baseline","NoMort"),
                 labels = c("No Emigration","Dispersal: Mortality","Dispersal: No Mortality"))
  ) %>%
  group_by(Year, PatchID, Scenario, Exp) %>%
  summarise(median_prob = median(extProb), .groups = "drop") %>%
  mutate(
    Year_bin = floor(Year / 10) * 10,
    Ext_bin  = cut(median_prob, breaks = seq(0, 1, by = 0.1),
                   include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(Scenario, Exp, Year_bin, Ext_bin) %>%
  summarise(n_patches = n_distinct(PatchID), .groups = "drop") %>%
  mutate(Ext_bin_num = as.numeric(sub("\\[|\\)", "", sub(",.*", "", Ext_bin))) + 0.05) %>%
  # Compute proportions
  group_by(Scenario, Exp, Year_bin) %>%
  mutate(prop_patches = n_patches / sum(n_patches)) %>%
  ungroup() |>
  tidyr::complete(
    Scenario,
    Exp,
    Year_bin = seq(0, 100, by = 10),                # your x bins
    Ext_bin_num = seq(0.05, 0.95, by = 0.1),         # y bin midpoints
    fill = list(prop_patches = 0)                      # fill missing with 0
  )


# Plot proportion heatmap
patch_prop_heatmap <- ggplot(sdl_p_sims_heatmap_prop, aes(x = Year_bin, y = Ext_bin_num, fill = prop_patches)) +
  geom_tile(na.rm = TRUE) +
  scale_fill_viridis_c(name = "Proportion of patches", option = "mako", limits = c(0,1)) +
  # scale_fill_gradient2(name = "Proportion of patches",
  #                      limits = c(0,1),
  #                      low = "blue",
  #                      mid = "white",
  #                      high = "red",
  #                      midpoint = 0.5) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  labs(x = "Years from present", y = "Extinction Probability") +
  facet_grid(vars(Exp), vars(Scenario)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, family = "sans"),
    axis.text = element_text(size = 12, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 14, family = "sans"),
    strip.text = element_text(size = 14, family = "sans", face = "bold"),
    strip.background = element_blank()
  )

heatmap_joint <- patch_heatmap / patch_prop_heatmap

ggsave(plot = heatmap_joint, filename = "R:/Users/skeyser/PhD/NASA Project/Writing/Marten_Chapter/NCC_Sub/MetaPatchHeatMapJoint.jpg",
       height = 12, width = 12, dpi = 800)


threshold <- 0.5  # extinction probability threshold

heat_df <- extProb %>%
  group_by(PatchID, Year) %>%
  summarise(
    n_models_over50 = sum(ExtProb >= threshold),     # count models >= threshold
    total_models = n_distinct(Model),
    .groups = "drop"
  ) %>%
  mutate(frac_models_over50 = n_models_over50 / total_models)  # optional fraction

## Proportion of sites exceeding 50% threshold across scenarios
exceed <- patch_extprob |>
  filter(Year %in% c(20, 60, 80)) |>
  filter(PatchID != 0 & !is.na(Region)) |>
  mutate(Scenario = factor(Scenario, levels = c("GFDL",
                                                "MRI",
                                                "MIROC5",
                                                "CNRM",
                                                "ACCESS",
                                                "IPSL"))) |>
  group_by(Region, Year, PatchID, Scenario) |>
  summarise(
    median_prob = median(extProb),
    .groups = "drop"
  ) |>
  group_by(Year, Scenario) |>
  mutate(Exceed = ifelse(median_prob > 0.5, 1, 0)) |>
  summarise(
    n_patches      = n_distinct(PatchID),
    n_exceed       = sum(Exceed, na.rm = TRUE),
    prop_exceed    = n_exceed / n_patches,
    .groups = "drop"
  )



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Patch extinction by landscape context
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Area vs. change in extinction prob
ext_area <- patch_extprob |>
  filter(PatchID != 0 & !is.na(Region)) |>
  mutate(Scenario = factor(Scenario, levels = c("GFDL",
                                                "MRI",
                                                "MIROC5",
                                                "CNRM",
                                                "ACCESS",
                                                "IPSL"))) |>
  mutate(Exp = factor(Exp, levels = c("NoEmig",
                                      "Baseline",
                                      "NoMort"),
                      labels = c("No Connectivity",
                                 "Limited Connectivity",
                                 "Max Connectivity"))) |>
  group_by(PatchID, Area, Region, Scenario, Exp) |>
  summarise(
    ext_initial = extProb[Year == min(Year)],
    ext_final = extProb[Year == max(Year)],
    delta_ext = ext_final - ext_initial,
    .groups = "drop"
  )

# Create the plot
range(ext_area$ext_initial, na.rm = T)

ext_mod <- ext_area |>
  mutate(Trials = 100,
         Extinctions = ext_final * 100)

area_mod <- glm(cbind(Extinctions, Trials - Extinctions) ~ log(Area), family = "quasibinomial", data = ext_mod)
summary(area_mod)

# Create a list to store results
models <- list()
scenarios <- expand.grid(
  climate = c("GFDL", "MRI", "MIROC5", "CNRM", "ACCESS", "IPSL"),
  dispersal = c("No Connectivity", "Limited Connectivity", "Max Connectivity")
)

# Fit separate models
for(i in 1:nrow(scenarios)) {
  subset_data <- ext_mod[ext_mod$Scenario == scenarios$climate[i] &
                            ext_mod$Exp == scenarios$dispersal[i], ]

  models[[paste(scenarios$climate[i], scenarios$dispersal[i], sep="_")]] <-
    glm(cbind(Extinctions, Trials - Extinctions) ~ log(Area),
        family = "quasibinomial", data = subset_data)
}

# Function to find area for 50% extinction probability
find_critical_area <- function(model) {
  # Logistic model: logit(p) = intercept + slope * log(Area)
  # For p = 0.5, logit(0.5) = 0
  # So: 0 = intercept + slope * log(Area)
  # Therefore: log(Area) = -intercept / slope

  intercept <- coef(model)[1]
  slope <- coef(model)[2]

  if(is.na(slope) || slope == 0) return(NA)

  log_area_50 <- -intercept / slope
  area_50 <- exp(log_area_50)

  return(area_50)
}

# Apply to all models
critical_areas <- map_dbl(models, find_critical_area)

# Create summary table
critical_area_table <- data.frame(
  Scenario = scenarios$climate,
  Dispersal = scenarios$dispersal,
  Critical_Area_ha = critical_areas
) %>%
  separate(names(critical_areas), into = c("Climate", "Dispersal_temp"), sep = "_", remove = FALSE) %>%
  select(-Dispersal_temp) %>%
  mutate(Critical_Area_ha = round(Critical_Area_ha, 1))

# Extract coefficients
coeffs <- map_dfr(models, ~tidy(.x), .id = "scenario") %>%
  filter(term == "log(Area)")

# Compare area effects across scenarios
ggplot(coeffs, aes(x = scenario, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = estimate - 1.96*std.error,
                    ymax = estimate + 1.96*std.error)) +
  theme(axis.text.x = element_text(angle = 45))

ext_area_plot <- ext_area |>
  filter(Scenario %in% c("GFDL", "IPSL")) |>
  ggplot(aes(x = log(as.numeric(Area)), y = ext_final)) +
  geom_point(aes(color = Region), size = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_bw() +
  scale_color_viridis_d() +
  labs(x = "Log Patch Area (km²)",
       y = "Extinction Probability") +
  scale_y_continuous(limits = c(0, 1)) +
  geom_smooth(method = "glm",
              method.args = list(family = "quasibinomial"),
              color = "black", se = T) +
  facet_grid(vars(Exp), vars(Scenario)) +
  theme(axis.title = element_text(size = 14, family = "sans"),
        axis.text = element_text(size = 12, family = "sans"),
        legend.text = element_text(size = 12, family = "sans"),
        legend.title = element_text(size = 14, family = "sans"),
        strip.text = element_text(size = 14, family = "sans", face = "bold"),
        strip.background = element_blank())

ggsave(plot = ext_area_plot, filename = "R:/Users/skeyser/PhD/NASA Project/Writing/Marten_Chapter/NCC_Sub/ExtProbVSarea_AllPatch_WLadj_DispScen_GLM_GFDL_IPSL_V2.jpg",
       height = 8, width = 8, dpi = 800)

## Isolation
ext_iso <- patch_extprob |>
  filter(PatchID != 0 & !is.na(Region)) |>
  group_by(PatchID, Region, Iso, Scenario, Exp) |>
  summarise(
    ext_initial = extProb[Year == min(Year)],
    ext_final = extProb[Year == max(Year)],
    delta_ext = ext_final - ext_initial,
    .groups = "drop"
  ) |>
  mutate(
    dist_bin = case_when(
      (Iso/1000) < 2 ~ "Close (<1.5)",
      (Iso/1000) < 5 ~ "Medium (1.5-2.5)",
      (Iso/1000) < 10 ~ "Far (2.5-3.5)",
      TRUE ~ "Very Far (>3.5)"
    )
  )

iso_mod <- ext_iso |>
  filter(Scenario == "IPSL" & Exp == "Baseline") |>
  mutate(Trials = 100,
         Extinctions = ext_final * 100)

imod <- lm(ext_final ~ dist_bin, data = iso_mod)
summary(imod)

ext_iso_plot <- ggplot(data = ext_iso, aes(x = dist_bin, y = ext_final)) +
  geom_point(aes(color = Region), size = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_smooth(method = "gam",
              method.args = list(family = "quasibinomial"),
              color = "black", se = T) +
  theme_bw() +
  scale_color_viridis_d() +
  labs(x = "Log Dist. Nearest Neighbor (km)",
       y = "Extinction Probability") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(vars(Scenario), vars(Exp))

ggsave(plot = ext_iso_plot, filename = "R:/Users/skeyser/PhD/NASA Project/Writing/Marten_Chapter/NCC_Sub/ExtProbVSiso_AllPatch_WLadj.jpg",
       height = 8, width = 12, dpi = 800)

## Latitude centroid
ext_lat <- patch_extprob |>
  filter(PatchID != 0 & !is.na(Region)) |>
  group_by(PatchID, Region, Lat, Scenario, Exp) |>
  summarise(
    ext_initial = extProb[Year == min(Year)],
    ext_final = extProb[Year == max(Year)],
    delta_ext = ext_final - ext_initial,
    .groups = "drop"
  )

ext_lat_plot <- ggplot(data = ext_lat, aes(x = Lat, y = ext_final)) +
  geom_point(aes(color = Region), size = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_smooth(method = "glm",
              method.args = list(family = "quasibinomial"),
              color = "black", se = T) +
  theme_bw() +
  scale_color_viridis_d() +
  labs(x = "Latitude",
       y = "Extinction Probability") +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(vars(Scenario), vars(Exp))

ggsave(plot = ext_is_plot, filename = "R:/Users/skeyser/PhD/NASA Project/Writing/Marten_Chapter/NCC_Sub/ExtProbVSiso_AllPatch_WLadj.jpg",
       height = 8, width = 12, dpi = 800)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection:
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Roi
roi <- st_as_sfc(st_bbox(patchShp))

## Raster plotting for the future scenarios
mcsnow <- rast(here("Data/MidCentury_PctChange.tif"))
lcsnow <- rast(here("Data/LateCentury_PctChange.tif"))
plot(lcsnow)
plot(mcsnow)

## Clip these by the shapefile
gl <- gl.state |> st_transform(crs = crs(lcsnow))
roi <- roi |> st_transform(crs = crs(lcsnow))

mcsnow.c <- crop(mcsnow, roi, mask = T)
mcsnow.c <- crop(mcsnow.c, gl, mask = T)

lcsnow.c <- crop(lcsnow, roi, mask = T)
lcsnow.c <- crop(lcsnow.c, gl, mask = T)
plot(lcsnow.c)

gl.c <- st_crop(gl, roi)


## Regional Maps
mcsnow.df <- terra::as.data.frame(mcsnow.c, xy = T)
mcgrid <- mcsnow.df |>
  pivot_longer(cols = contains("mid21"), names_to = "Scenario", values_to = "pct_change") |>
  mutate(Scenario = factor(case_when(Scenario == "access_mid21" ~ "ACCESS",
                              Scenario == "cnrm_mid21" ~ "CNRM",
                              Scenario == "gfdl_mid21" ~ "GFDL",
                              Scenario == "ipsl_mid21" ~ "IPSL",
                              Scenario == "miroc5_mid21" ~ "MIROC5",
                              Scenario == "mri_mid21" ~ "MRI"),
         levels = c("GFDL",
                    "MRI",
                    "MIROC5",
                    "CNRM",
                    "ACCESS",
                    "IPSL"))) |>
  ggplot() +
  geom_raster(aes(x=x, y=y, fill = pct_change)) +
  theme_void() +
  scale_fill_gradient2(low = scales::muted("red"),
                       mid = "white",
                       high = scales::muted("blue"),
                       midpoint = 0,
                       limits = c(-68, 14),
                       breaks = c(-60, -40, -20, 0, 20)) +
  labs(fill = "% Change SWE \n(Relative to Baseline 1990-2009)") +
  facet_wrap(~Scenario) +
  ggtitle("Mid-century Snow Projections (2030-2049)")

lcsnow.df <- terra::as.data.frame(lcsnow.c, xy = T)
lcgrid <- lcsnow.df |>
  pivot_longer(cols = contains("late21"), names_to = "Scenario", values_to = "pct_change") |>
  mutate(Scenario = factor(case_when(Scenario == "access_late21" ~ "ACCESS",
                                     Scenario == "cnrm_late21" ~ "CNRM",
                                     Scenario == "gfdl_late21" ~ "GFDL",
                                     Scenario == "ipsl_late21" ~ "IPSL",
                                     Scenario == "miroc5_late21" ~ "MIROC5",
                                     Scenario == "mri_late21" ~ "MRI"),
                           levels = c("GFDL",
                                      "MRI",
                                      "MIROC5",
                                      "CNRM",
                                      "ACCESS",
                                      "IPSL"))) |>
  ggplot() +
  geom_raster(aes(x=x, y=y, fill = pct_change)) +
  theme_void() +
  scale_fill_gradient2(low = scales::muted("red"),
                       mid = "white",
                       high = scales::muted("blue"),
                       midpoint = 0,
                       limits = c(-68, 14),
                       breaks = c(-60, -40, -20, 0, 20)) +
  labs(fill = "% Change SWE \n(Relative to Baseline 1990-2009)") +
  facet_wrap(~Scenario) +
  ggtitle("Late-century Snow Projections (2080-2099)")

snowgrid <- mcgrid / lcgrid + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom",
                                                                      legend.title = element_text(vjust = 0.5,
                                                                                                  hjust = 0.5))

ggsave(snowgrid, filename = here("Figures/SI_Fig_SnowProjectionGrid.jpg"),
       dpi = 600, height = 8, width = 8)

## Make the raster a df
lcsnow.df <- as.data.frame(lcsnow.c, xy = TRUE)

ggplot(data = lcsnow.df) +
  geom_raster(aes(x=x, y=y, fill = ipsl_late21)) +
  scale_fill_gradientn(colors = c("#c6dbef",
                                  "#9ecae1","#6baed6","#3182bd","#08519c")) +
  theme_bw() +
  labs(fill = "IPSL Late-Century \nSWE (% Change)") +
  theme_void()

ggplot(data = mcsnow.df) +
  geom_raster(aes(x=x, y=y, fill = gfdl_late21)) +
  scale_fill_gradientn(colors = c("#f7fbff","#deebf7","#c6dbef",
                                  "#9ecae1","#6baed6","#3182bd","#08519c")) +
  geom_sf(data = gl.c, size = 1, color = "darkgray", fill = NA) +
  theme_bw() +
  labs(fill = "GFDL Late-Century \nSWE (% Change)") +
  theme_void()


library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# --- Read rasters ---

baseline <- rast(here("Data/Baseline_PctChange.tif"))
mcsnow <- rast(here("Data/MidCentury_PctChange.tif"))
lcsnow <- rast(here("Data/LateCentury_PctChange.tif"))

# ROI from patch shapefile
roi <- st_as_sfc(st_bbox(patchShp)) |> st_transform(crs = crs(mcsnow))
gl <- gl.state |> st_transform(crs = crs(mcsnow))
gl <- st_crop(gl, roi)

# --- Crop/mask rasters to ROI ---
bsnow.c <- crop(baseline, roi, mask = TRUE) |> crop(gl, mask = TRUE)
mcsnow.c <- crop(mcsnow, roi, mask = TRUE) |> crop(gl, mask = TRUE)
lcsnow.c <- crop(lcsnow, roi, mask = TRUE) |> crop(gl, mask = TRUE)

# --- Convert to dataframes ---
bsnow.df <- as.data.frame(bsnow.c, xy = TRUE)
mcsnow.df <- as.data.frame(mcsnow.c, xy = TRUE)
lcsnow.df <- as.data.frame(lcsnow.c, xy = TRUE)

# --- Combine mid and late century for GFDL + IPSL ---
snow_df <- bind_rows(
  bsnow.df |>
    transmute(x, y,
              Value = ipsl_late20, Scenario = "IPSL", Horizon = "Baseline"),
  mcsnow.df %>%
    transmute(x, y,
              Value = ipsl_mid21, Scenario = "IPSL", Horizon = "Mid-century"),
  lcsnow.df %>%
    transmute(x, y,
              Value = ipsl_late21, Scenario = "IPSL", Horizon = "Late-century")
)

# --- Ensure facet order ---
snow_df <- snow_df %>%
  mutate(
    Horizon = forcats::fct_relevel(Horizon, "Baseline", "Mid-century", "Late-century"),
  )

# Separate baseline and change data
baseline_df <- snow_df %>% filter(Horizon == "Baseline")
change_df   <- snow_df %>% filter(Horizon != "Baseline")

ggplot() +
  # --- Baseline SWE
  geom_raster(data = baseline_df,
              aes(x = x, y = y, fill = Value)) +
  scale_fill_viridis_c(name = expression("SWE (kg/m"^2*")"), option = "G") +
  # If you want fixed limits:
  # scale_fill_viridis_c(name = "SWE (mm)", option = "C", limits = c(0, max_SWE))

  new_scale_fill() +  # This lets us add a second fill scale

  # --- Percent change layers
  geom_raster(data = change_df,
              aes(x = x, y = y, fill = Value)) +
  scale_fill_gradientn(
    colours = c("#ffffff", "#fddbc7", "#ef8a62", "#b2182b"),
    values = scales::rescale(c(0, -10, -30, -50)),
    limits = c(-60, 10),
    name   = "SWE % Change"
  ) +

  facet_grid(~Horizon) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.margin = unit(c(0, 0, 0, 0), "lines")
  )
ggsave(filename = "R:/Users/skeyser/PhD/NASA Project/Writing/Marten_Chapter/NCC_Sub/FutureSnowMap.jpg",
       height = 6, width = 8, dpi = 800)

library(dplyr)
library(forcats)

library(dplyr)
library(forcats)

# Make sure all rows have the same column name: Value
snow_df <- bind_rows(
  bsnow.df %>%
    transmute(x, y,
              Value = ipsl_late20, Scenario = "IPSL", Horizon = "Baseline"),
  mcsnow.df %>%
    transmute(x, y,
              Value = ipsl_mid21, Scenario = "IPSL", Horizon = "Mid-century"),
  lcsnow.df %>%
    transmute(x, y,
              Value = ipsl_late21, Scenario = "IPSL", Horizon = "Late-century")
) %>%
  mutate(Horizon = forcats::fct_relevel(Horizon, "Baseline", "Mid-century", "Late-century"))

# Now create the Category column
snow_df <- snow_df %>%
  mutate(
    Category = case_when(
      Horizon == "Baseline" & Value < 150 ~ "<150",
      Horizon == "Baseline" & Value < 200 ~ "150–200",
      Horizon == "Baseline" & Value < 250 ~ "200–250",
      Horizon == "Baseline"               ~ ">250",

      Horizon != "Baseline" & Value < -50 ~ "< -50%",
      Horizon != "Baseline" & Value < -40 ~ "-50% to -40%",
      Horizon != "Baseline" & Value < -30 ~ "-40% to -30%",
      Horizon != "Baseline" & Value < -20 ~ "-30% to -20%",
      Horizon != "Baseline" & Value < -10 ~ "-20% to -10%",
      Horizon != "Baseline"               ~ "> -10%"
    ),
    Category = factor(Category, levels = c("<150", "150–200", "200–250", ">250",
                                           "< -50%", "-50% to -40%", "-40% to -30%", "-30% to -20%", "-20% to -10%", "> -10%"))
  )

baseline_df <- snow_df %>% filter(Horizon == "Baseline")
change_df   <- snow_df %>% filter(Horizon != "Baseline")

ggplot() +
  # --- Baseline SWE colors: Blues ---
  geom_raster(data = baseline_df, aes(x = x, y = y, fill = Category)) +
  scale_fill_manual(
    name = "SWE (kg/m²)",
    values = c(
      "<150"    = "#c6dbef",
      "150–200" = "#6baed6",
      "200–250" = "#428AC2",
      ">250"    = "#08306b"
    ),
    drop = FALSE
  ) +

  new_scale_fill() +

  # --- Percent change colors: Reds ---
  geom_raster(data = change_df, aes(x = x, y = y, fill = Category)) +
  scale_fill_manual(
    name = "SWE % Change",
    values = c(
      "< -50%"       = "#bd0026",
      "-50% to -40%" = "#f03b20",
      "-40% to -30%" = "#fd8d3c",
      "-30% to -20%" = "#feb24c",
      "-20% to -10%" = "#fed976",
      "-10%" = "#ffffb2"

    ),
    drop = FALSE
  ) +

  facet_grid(~Horizon) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "bottom",      # keep both legends together
    legend.box = "vertical",         # stack them vertically
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width  = unit(2.0, "cm"),
    strip.text = element_text(size = 12, face = "bold"),
    #panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
  )

ggsave(filename = "R:/Users/skeyser/PhD/NASA Project/Writing/Marten_Chapter/NCC_Sub/FutureSnowMap.jpg",
       height = 6, width = 10, dpi = 800)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Global Inset Map
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# Orthographic projection centered on Great Lakes
ortho_crs <- "+proj=ortho +lat_0=44 +lon_0=-85 +R=6371000 +units=m +no_defs"

# Countries
world <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(crs = ortho_crs)

# Make graticule in EPSG:4326
grat <- st_graticule(lat = seq(-80, 80, 40),
                     lon = seq(-180, 180, 20),
                     crs = 4326) %>%
  st_transform(ortho_crs) %>%
  st_segmentize(dfMaxLength = 10000)  # densify so vertices ~10 km apart

# Build sphere mask *in projected coords*
sphere <- st_sfc(st_point(c(0,0)), crs = ortho_crs) %>%
  st_buffer(6371000)

# Clip densified graticule
grat <- st_intersection(grat, sphere)

# ROI bounding box
roi <- st_as_sfc(
  st_bbox(c(xmin = -92, ymin = 42, xmax = -82, ymax = 49), crs = 4326)
) |> st_transform(crs = ortho_crs)

# Great Lakes polygons
lakes <- rnaturalearth::ne_download(scale = 50, type = "lakes", category = "physical",
                                    returnclass = "sf")

great_lakes <- lakes %>%
  filter(name %in% c("Lake Superior", "Lake Michigan", "Lake Huron", "Lake Erie", "Lake Ontario")) %>%
  st_transform(crs = ortho_crs)

# Mask land to sphere boundary
world <- st_make_valid(world) |> st_crop(st_bbox(sphere))
world_masked <- st_intersection(world, sphere)

ggplot() +
  # Light ocean background (sphere)
  geom_sf(data = sphere, fill = "white", color = "grey70") +
  # Graticules
  geom_sf(data = grat, color = "grey70", linewidth = 0.25) +
  # Land masked to sphere
  geom_sf(data = world_masked, fill = "grey20", color = "black", size = 0.2) +
  # Great Lakes overlay
  geom_sf(data = great_lakes, fill = "white", color = "black", size = 0.25) +
  # ROI outline
  geom_sf(data = roi, fill = NA, color = "red", linewidth = 0.5) +
  # Keep whole sphere visible
  coord_sf(crs = ortho_crs,
           xlim = c(-6371000, 6371000),  # ~full radius
           ylim = c(-6371000, 6371000)) +
  theme_void()

ggsave(filename = "R:/Users/skeyser/PhD/NASA Project/Writing/Marten_Chapter/NCC_Sub/GlobalInset.jpg",
       height = 6, width = 6, dpi = 600)

## Basemap with snow season length
wl <- rast("C:/Users/skeyser/Downloads/SSL_WHIs_mean20032020/SSL_MOD_Daily_cUSA_mean20032020_WGS84_14d.tif")
gl.t <- st_transform(gl.state, crs = crs(wl)) |> vect()

wl.c <- crop(wl, gl.t)
plot(wl.c)
