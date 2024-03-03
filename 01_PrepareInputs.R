# This script loads and prepares the available scenario data from the various 
# raw datasets. The raw datasets need to be downloaded separately.

# For each dataset the respective layer for 2050 is saved as well as a 
# summarized time series summarized as statistical moments at global scale.

# -------------- #
# Load the function sourcing script
source('00_Function.R')
# -------------- #

#### Biodiversity - Nature Conservation ----
# Terrestrial
# Load Jung et al. for 30% expansion
ras_naturecons <- rast(paste0(path_rawdata,"BiodiversityOnly/50km/minshort_speciestargets_biome.id_withPA_esh50km_repruns10_ranked.tif"))
ras_naturecons <- project(ras_naturecons, background)

new <- terra::deepcopy(ras_naturecons)
new[new>30] <- 0 # Get top 30%
new[new>0] <- 1
writeRaster(new, paste0(path_processed, "biodiv_natcons_terr_low.tif"),overwrite=T)

new <- terra::deepcopy(ras_naturecons)
new[new>50] <- 0 # Get top 30%
new[new>0] <- 1
writeRaster(new, paste0(path_processed, "biodiv_natcons_terr_high.tif"), overwrite=T)

# Save a summary, in this case simply the baseline and assumption of full protection
# by 2050
out <- data.frame(nexus = "Biodiversity", entrypoint = "Nature conservation",
                  realm = "terrestrial",
                  scenario = c("low", "low", "high", "high"),
                  year = c(2020, 2050),
                  mean = c(0.19, .30, 0.19, .5)
                  )
write.csv2(out, paste0(path_processed, "biodiv_natcons_terr.csv"))

# Marine
# Marine nature conservation 2050 from the estimates of Sala et al.
ras_marinecons <- rast(paste0(path_rawdata,"SalaEtAl_marineCons/biodiversity_ranking_2050.tif"))
ras_marinecons <- project(ras_marinecons, background)
ar <- terra::cellSize(ras_marinecons, unit = "km")
ar <- terra::mask(ar, background,inverse=TRUE) # Get only marine area

new <- terra::deepcopy(ras_marinecons)
# terra::global(new * ar, "sum", na.rm = TRUE) / terra::global(ar, "sum", na.rm = TRUE)
new[new<0.95] <- 0 # We get the top 5% of values for the highest ranked priority cells
new[new>0] <- 1
writeRaster(new, paste0(path_processed, "biodiv_natcons_marine_low.tif"),overwrite= TRUE)

new <- terra::deepcopy(ras_marinecons)
new[new<0.8] <- 0
new[new>0] <- 1
writeRaster(new, paste0(path_processed, "biodiv_natcons_marine_high.tif"),overwrite= TRUE)

#### Biodiversity - Biodiversity targets ----
# Here we use estimates from BES-SIM
# Community intactness
ifname <- paste0(path_rawdata, "EBV_BES-SIM/pereira_comcom_id28_20231212_v2.nc")
ebv_datacubepaths(ifname)
ebv_trend(filepath = ifname,datacubepath = "scenario_1/metric_4/ebv_cube",entity = 1)
new <- ebvcube::ebv_read(filepath = ifname,datacubepath = "scenario_1/metric_4/ebv_cube",
                  entity = 1,timestep = 2:3,
                  type = "r")
terra::time(new) <- as.POSIXct(c("2015-01-01", "2050-01-01"))
new <- project(new, background)

out <- data.frame(nexus = "Biodiversity", entrypoint = "Biodiversity targets",
                  realm = "terrestrial",
                  scenario = c("high", "high"),
                  year = c(2020, 2050)
) |> dplyr::mutate(
  mean = terra::global(new, "mean", na.rm = TRUE)[,1],
  median =  terra::global(new, median, na.rm = TRUE)[,1],
  min = terra::global(new, "min", na.rm = TRUE)[,1],
  max = terra::global(new, "max", na.rm = TRUE)[,1]
)
write.csv2(out, paste0(path_processed, "biodiv_biodtarget_terr_high.csv"))

new <- ibis.iSDM::predictor_transform(new, option = "norm") # Normalize to counter reprojection issues
new <- 1 - new # Invert to ensure 
writeRaster(new[[2]], paste0(path_processed, "biodiv_biodtarget_terr_high.tif"),overwrite= TRUE)

# --- #
# Low ambition
ebv_datacubepaths(ifname)
ebv_trend(filepath = ifname,datacubepath = "scenario_3/metric_4/ebv_cube",entity = 1)
new <- ebvcube::ebv_read(filepath = ifname,datacubepath = "scenario_3/metric_4/ebv_cube",
                         entity = 1,timestep = 2:3,
                         type = "r")
terra::time(new) <- as.POSIXct(c("2015-01-01", "2050-01-01"))
new <- project(new, background)

out <- data.frame(nexus = "Biodiversity", entrypoint = "Biodiversity targets",
                  realm = "terrestrial",
                  scenario = c("low", "low"),
                  year = c(2020, 2050)
) |> dplyr::mutate(
  mean = terra::global(new, "mean", na.rm = TRUE)[,1],
  median =  terra::global(new, median, na.rm = TRUE)[,1],
  min = terra::global(new, "min", na.rm = TRUE)[,1],
  max = terra::global(new, "max", na.rm = TRUE)[,1]
)
write.csv2(out, paste0(path_processed, "biodiv_biodtarget_terr_low.csv"))

new <- ibis.iSDM::predictor_transform(new, option = "norm") # Normalize to counter reprojection issues
new <- 1 - new # Invert to ensure 
writeRaster(new[[2]], paste0(path_processed, "biodiv_biodtarget_terr_low.tif"),overwrite= TRUE)

#### Water - Quality ----
# Maybe from here? https://zenodo.org/records/7811612

#### Water - Supply ----
# Total potential water supply
# Here using a cWATM run for a single GCM
# High ambition
wat <- read_ncdf("raw_data/isimip/cwatm_gfdl-esm4_w5e5_ssp126_2015soc-from-histsoc_default_qtot_global_monthly_2015_2100.nc", proxy= FALSE)

# Get reference year 2015
ref <- ibis.iSDM:::stars_to_raster(wat, which= 1:12)
ref <- Reduce("+", ref)
ref <- terra::project(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(wat, which= 421:432)
new <- Reduce("+", new)
new <- terra::project(new, background)

# Relative change
out <- (new - ref) / ref
out <- predictor_transform(out, option = "windsor",windsor_props = c(0,.99))
writeRaster(out, paste0(path_processed, "water_supply_terr_high.tif"), overwrite= TRUE)

# Summarize 
df <- as.data.frame(wat) |> tidyr::drop_na()
df$year <- lubridate::year(df$time)
df$qtot <- units::drop_units(df$qtot)
df <- df |> group_by(year) |> summarise(sum = sum(qtot),
                                        min = min(qtot), max = max(qtot),
                                        mean = mean(qtot), median = median(qtot))
gc()
df <- df |> dplyr::mutate(nexus = "Water", entrypoint = "Supply",
                          realm = "terrestrial",
                          scenario = "high")
write.csv2(df, paste0(path_processed, "water_supply_terr_high.csv"))

# --- # 
# low ambition
wat <- read_ncdf("raw_data/isimip/cwatm_gfdl-esm4_w5e5_ssp585_2015soc-from-histsoc_default_qtot_global_monthly_2015_2100.nc", proxy= FALSE)

# Get reference year 2015
ref <- ibis.iSDM:::stars_to_raster(wat, which= 1:12)
ref <- Reduce("+", ref)
ref <- terra::project(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(wat, which= 421:432)
new <- Reduce("+", new)
new <- terra::project(new, background)

# Relative change
out <- (new - ref) / ref
out <- predictor_transform(out, option = "windsor",windsor_props = c(0,.99))
writeRaster(out, paste0(path_processed, "water_supply_terr_low.tif"), overwrite= TRUE)

# Summarize 
df <- as.data.frame(wat) |> tidyr::drop_na()
df$year <- lubridate::year(df$time)
df$qtot <- units::drop_units(df$qtot)
df <- df |> group_by(year) |> summarise(sum = sum(qtot),
                                        min = min(qtot), max = max(qtot),
                                        mean = mean(qtot), median = median(qtot))
gc()
df <- df |> dplyr::mutate(nexus = "Water", entrypoint = "Supply",
                          realm = "terrestrial",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "water_supply_terr_low.csv"))

#### Water - Demand ----
# Total actual water withdrawal
# Here using a cWATM run for a single GCM
# High ambition
wat <- read_ncdf("raw_data/isimip/cwatm_gfdl-esm4_w5e5_ssp126_2015soc-from-histsoc_default_atotww_global_monthly_2015_2100.nc", proxy= FALSE)

# Get reference year 2015
ref <- ibis.iSDM:::stars_to_raster(wat, which= 1:12)
ref <- Reduce("+", ref)
ref <- terra::project(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(wat, which= 421:432)
new <- Reduce("+", new)
new <- terra::project(new, background)

# Relative change
out <- (new - ref) / ref
out <- predictor_transform(out, option = "windsor",windsor_props = c(0,.99))
writeRaster(out, paste0(path_processed, "water_demand_terr_high.tif"), overwrite= TRUE)

# Summarize 
df <- as.data.frame(wat) |> tidyr::drop_na()
df$year <- lubridate::year(df$time)
df$atotww <- units::drop_units(df$atotww)
df <- df |> group_by(year) |> summarise(sum = sum(atotww),
                                        min = min(atotww), max = max(atotww),
                                        mean = mean(atotww), median = median(atotww))
gc()
df <- df |> dplyr::mutate(nexus = "Water", entrypoint = "Demand",
                          realm = "terrestrial",
                          scenario = "high")
write.csv2(df, paste0(path_processed, "water_demand_terr_high.csv"))

# --- # 
# low ambition
wat <- read_ncdf("raw_data/isimip/cwatm_gfdl-esm4_w5e5_ssp585_2015soc-from-histsoc_default_atotww_global_monthly_2015_2100.nc", proxy= FALSE)

# Get reference year 2015
ref <- ibis.iSDM:::stars_to_raster(wat, which= 1:12)
ref <- Reduce("+", ref)
ref <- terra::project(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(wat, which= 421:432)
new <- Reduce("+", new)
new <- terra::project(new, background)

# Relative change
out <- (new - ref) / ref
out <- predictor_transform(out, option = "windsor",windsor_props = c(0,.99))
writeRaster(out, paste0(path_processed, "water_demand_terr_low.tif"), overwrite= TRUE)

# Summarize 
df <- as.data.frame(wat) |> tidyr::drop_na()
df$year <- lubridate::year(df$time)
df$atotww <- units::drop_units(df$atotww)
df <- df |> group_by(year) |> summarise(sum = sum(atotww),
                                        min = min(atotww), max = max(atotww),
                                        mean = mean(atotww), median = median(atotww))
gc()
df <- df |> dplyr::mutate(nexus = "Water", entrypoint = "Demand",
                          realm = "terrestrial",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "water_demand_terr_low.csv"))

#### Health - Infectious diseases ----

#### Food - Supply ----
n <- read_ncdf("raw_data/isimip/lpjml_gfdl-esm4_w5e5_ssp126_2015soc_default_biom-swh-firr_global_annual-gs_2015_2100.nc",proxy = F,
               make_time = FALSE) 
n <- stars::st_set_dimensions(n, 'time',values = seq.Date(as.Date('2015-01-01'), as.Date('2100-01-01'),by = "year"))

# Get reference year
ref <- ibis.iSDM:::stars_to_raster(n, which = 1)[[1]]
ref <- terra::project(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(n, which= 36)[[1]]
new <- terra::project(new, background)

# Relative change
out <- (new - ref) / ref
out <- predictor_transform(out, option = "windsor",windsor_props = c(0,.99))
writeRaster(out, paste0(path_processed, "food_supply_terr_high.tif"), overwrite= TRUE)

# Summarize 
df <- as.data.frame(n) |> tidyr::drop_na()
df$year <- lubridate::year(df$time)
df$biom.swh.firr <- units::drop_units(df$biom.swh.firr)
df <- df |> group_by(year) |> summarise(sum = sum(biom.swh.firr),
                                        min = min(biom.swh.firr), max = max(biom.swh.firr),
                                        mean = mean(biom.swh.firr), median = median(biom.swh.firr))
gc()
df <- df |> dplyr::mutate(nexus = "Food", entrypoint = "Supply",
                          realm = "terrestrial",
                          scenario = "high")
write.csv2(df, paste0(path_processed, "food_supply_terr_high.csv"))

# --- #
# Low ambition
n <- read_ncdf("raw_data/isimip/lpjml_gfdl-esm4_w5e5_ssp585_2015soc_default_biom-swh-firr_global_annual-gs_2015_2100.nc",proxy = F,
               make_time = FALSE) 
n <- stars::st_set_dimensions(n, 'time',values = seq.Date(as.Date('2015-01-01'), as.Date('2100-01-01'),by = "year"))

# Get reference year
ref <- ibis.iSDM:::stars_to_raster(n, which = 1)[[1]]
ref <- terra::project(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(n, which= 36)[[1]]
new <- terra::project(new, background)

# Relative change
out <- (new - ref) / ref
out <- predictor_transform(out, option = "windsor",windsor_props = c(0,.99))
writeRaster(out, paste0(path_processed, "food_supply_terr_low.tif"), overwrite= TRUE)

# Summarize 
df <- as.data.frame(n) |> tidyr::drop_na()
df$year <- lubridate::year(df$time)
df$biom.swh.firr <- units::drop_units(df$biom.swh.firr)
df <- df |> group_by(year) |> summarise(sum = sum(biom.swh.firr),
                                        min = min(biom.swh.firr), max = max(biom.swh.firr),
                                        mean = mean(biom.swh.firr), median = median(biom.swh.firr))
gc()
df <- df |> dplyr::mutate(nexus = "Food", entrypoint = "Supply",
                          realm = "terrestrial",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "food_supply_terr_low.csv"))

# ---- #
# Marine Food supply of total catch from Boats model
n <- read_ncdf("raw_data/isimip/boats_gfdl-esm4_nobasd_ssp126_2015soc-from-histsoc_default_tc_global_monthly_2015_2100.nc",proxy = F,
               make_time = F) 

# Get reference year
ref <- ibis.iSDM:::stars_to_raster(n, which = 1:12)
ref <- Reduce("+", ref)
ref <- terra::project(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(n, which= 421:432)
new <- Reduce("+", new)
new <- terra::project(new, background)

# Relative change
out <- (new - ref) / ref
writeRaster(out, paste0(path_processed, "food_supply_mar_high.tif"), overwrite= TRUE)

# Summarize 
n <- read_ncdf("raw_data/isimip/boats_gfdl-esm4_nobasd_ssp126_2015soc-from-histsoc_default_tc_global_monthly_2015_2100.nc",proxy = F,
               make_time = T) 
df <- as.data.frame(n) |> tidyr::drop_na()
df$time <- as.character(df$time)
df$year <- lubridate::year(df$time)
df$tc <- units::drop_units(df$tc)
df <- df |> group_by(year) |> summarise(sum = sum(tc),
                                        min = min(tc), max = max(tc),
                                        mean = mean(tc), median = median(tc))
gc()
df <- df |> dplyr::mutate(nexus = "Food", entrypoint = "Supply",
                          realm = "marine",
                          scenario = "high")
write.csv2(df, paste0(path_processed, "food_supply_mar_high.csv"))
rm(df);gc()

# --- #
# Load ambition
n <- read_ncdf("raw_data/isimip/boats_gfdl-esm4_nobasd_ssp585_2015soc-from-histsoc_default_tc_global_monthly_2015_2100.nc",proxy = F,
               make_time = F) 

# Get reference year
ref <- ibis.iSDM:::stars_to_raster(n, which = 1:12)
ref <- Reduce("+", ref)
ref <- terra::project(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(n, which= 421:432)
new <- Reduce("+", new)
new <- terra::project(new, background)

# Relative change
out <- (new - ref) / ref
writeRaster(out, paste0(path_processed, "food_supply_mar_low.tif"), overwrite= TRUE)

# Summarize 
n <- read_ncdf("raw_data/isimip/boats_gfdl-esm4_nobasd_ssp585_2015soc-from-histsoc_default_tc_global_monthly_2015_2100.nc",proxy = F,
               make_time = T) 
df <- as.data.frame(n) |> tidyr::drop_na()
df$time <- as.character(df$time)
df$year <- lubridate::year(df$time)
df$tc <- units::drop_units(df$tc)
df <- df |> group_by(year) |> summarise(sum = sum(tc),
                                        min = min(tc), max = max(tc),
                                        mean = mean(tc), median = median(tc))
gc()
df <- df |> dplyr::mutate(nexus = "Food", entrypoint = "Supply",
                          realm = "marine",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "food_supply_mar_low.csv"))
rm(df);gc()

#### Climate - Impacts ----
# Use multi-sectoral risk data from IIASA hotspot explorer https://hotspots-explorer.org/
# These data are described in full in https://doi.org/10.1088/1748-9326/aabf45
# Here we use the multi-sectoral climate change impact risk
ras_climaterisk  <- stars::read_ncdf("raw_data/Hotspot/multisector-all-data-SSP1_1p5.nc") |> 
  dplyr::select(Multisector_SSP1_1p5_Score)

# To SpatRaster
ras_climaterisk <- rast(ras_climaterisk) # predictor_transform(option = "norm")
ras_climaterisk <- terra::project(ras_climaterisk, background)

## Convert to vulnerability
#ras_climaterisk <- (ras_climaterisk * -1) |> predictor_transform(method = "norm")
writeRaster(ras_climaterisk, paste0(path_processed, "climate_impact_terr_high.tif"), overwrite= TRUE)

df <- data.frame(year = c(2020, 2050),
                 sum = c(0, terra::global(ras_climaterisk, "sum", na.rm = TRUE)[,1]),
                 mean = c(0, terra::global(ras_climaterisk, "mean", na.rm = TRUE)[,1]),
                 min = c(0, terra::global(ras_climaterisk, "min", na.rm = TRUE)[,1]),
                 max = c(0, terra::global(ras_climaterisk, "max", na.rm = TRUE)[,1])
                 )
df <- df |> dplyr::mutate(nexus = "Climate", entrypoint = "Impact",
                          realm = "terrestrial",
                          scenario = "high")
write.csv2(df, paste0(path_processed, "climate_impact_terr_high.csv"))

# Low ambition
ras_climaterisk  <- stars::read_ncdf("raw_data/Hotspot/multisector-all-data-SSP3_3p0.nc") |> 
  dplyr::select(Multisector_SSP3_3p0_Score)

# To SpatRaster
ras_climaterisk <- rast(ras_climaterisk) # predictor_transform(option = "norm")
ras_climaterisk <- terra::project(ras_climaterisk, background)
# Convert to vulnerability
#ras_climaterisk <- (ras_climaterisk * -1) |> predictor_transform(method = "norm")
writeRaster(ras_climaterisk, paste0(path_processed, "climate_impact_terr_low.tif"), overwrite= TRUE)

df <- data.frame(year = c(2020, 2050),
                 sum = c(0, terra::global(ras_climaterisk, "sum", na.rm = TRUE)[,1]),
                 mean = c(0, terra::global(ras_climaterisk, "mean", na.rm = TRUE)[,1]),
                 min = c(0, terra::global(ras_climaterisk, "min", na.rm = TRUE)[,1]),
                 max = c(0, terra::global(ras_climaterisk, "max", na.rm = TRUE)[,1])
)
df <- df |> dplyr::mutate(nexus = "Climate", entrypoint = "Impact",
                          realm = "terrestrial",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "climate_impact_terr_low.csv"))


#### Climate - Mitigation ----
# Bioenergy production
# https://portal.geobon.org/ebv-detail?id=60
ifname <- paste0(path_rawdata, "EBV_BES-SIM/pereira_ecoser_id60_20231204_v2.nc")
ebv_datacubepaths(ifname)
ebv_trend(filepath = ifname,datacubepath = "scenario_1/metric_6/ebv_cube",entity = 1)

new <- ebvcube::ebv_read(filepath = ifname,datacubepath = "scenario_1/metric_6/ebv_cube",
                         entity = 1,timestep = 2:3,
                         type = "r")
terra::time(new) <- as.POSIXct(c("2015-01-01", "2050-01-01"))
new <- project(new, background)

out <- data.frame(nexus = "Climate", entrypoint = "Climate mitigation",
                  realm = "terrestrial",
                  scenario = c("high", "high"),
                  year = c(2020, 2050)
) |> dplyr::mutate(
  mean = terra::global(new, "mean", na.rm = TRUE)[,1],
  median =  terra::global(new, median, na.rm = TRUE)[,1],
  min = terra::global(new, "min", na.rm = TRUE)[,1],
  max = terra::global(new, "max", na.rm = TRUE)[,1]
)
write.csv2(out, paste0(path_processed, "climate_climatemitig_terr_high.csv"))

# Calculate relative change
o <- (new[[2]] - new[[1]]) / new[[1]]
o <- terra::clamp(o,lower = -1, upper = 1) # Clamped
writeRaster(o, paste0(path_processed, "climate_climatemitig_terr_high.tif"),overwrite= TRUE)

# Low ambition
ebv_datacubepaths(ifname)
ebv_trend(filepath = ifname,datacubepath = "scenario_3/metric_6/ebv_cube",entity = 1)

new <- ebvcube::ebv_read(filepath = ifname,datacubepath = "scenario_3/metric_6/ebv_cube",
                         entity = 1,timestep = 2:3,
                         type = "r")
terra::time(new) <- as.POSIXct(c("2015-01-01", "2050-01-01"))
new <- project(new, background)

out <- data.frame(nexus = "Climate", entrypoint = "Climate mitigation",
                  realm = "terrestrial",
                  scenario = c("low", "low"),
                  year = c(2020, 2050)
) |> dplyr::mutate(
  mean = terra::global(new, "mean", na.rm = TRUE)[,1],
  median =  terra::global(new, median, na.rm = TRUE)[,1],
  min = terra::global(new, "min", na.rm = TRUE)[,1],
  max = terra::global(new, "max", na.rm = TRUE)[,1]
)
write.csv2(out, paste0(path_processed, "climate_climatemitig_terr_low.csv"))

# Calculate relative change
o <- (new[[2]] - new[[1]]) / new[[1]]
o <- terra::clamp(o,lower = -1, upper = 1) # Clamped
writeRaster(o, paste0(path_processed, "climate_climatemitig_terr_low.tif"),overwrite= TRUE)


#### Climate - Adaptation ----
# Marine Costal risk
# https://portal.geobon.org/ebv-detail?id=63
ifname <- paste0(path_rawdata, "EBV_BES-SIM/pereira_ecoser_id63_20231207_v2.nc")
ebv_datacubepaths(ifname)
ebv_trend(filepath = ifname,datacubepath = "scenario_1/metric_6/ebv_cube",entity = 1)
new <- ebvcube::ebv_read(filepath = ifname,datacubepath = "scenario_1/metric_6/ebv_cube",
                         entity = 1,timestep = 2:3,
                         type = "r")
terra::time(new) <- as.POSIXct(c("2015-01-01", "2050-01-01"))
new <- project(new, background)

# Convert to vulnerability
#new <- (new * -1) |> predictor_transform(method = "norm")

out <- data.frame(nexus = "Climate", entrypoint = "Climate adaptation",
                  realm = "marine",
                  scenario = c("high", "high"),
                  year = c(2020, 2050)
) |> dplyr::mutate(
  mean = terra::global(new, "mean", na.rm = TRUE)[,1],
  median =  terra::global(new, median, na.rm = TRUE)[,1],
  min = terra::global(new, "min", na.rm = TRUE)[,1],
  max = terra::global(new, "max", na.rm = TRUE)[,1]
)
write.csv2(out, paste0(path_processed, "climate_climateadapt_mar_high.csv"))

# Calculate relative change
o <- (new[[2]] - new[[1]]) / new[[1]]
writeRaster(o, paste0(path_processed, "climate_climateadapt_mar_high.tif"),overwrite= TRUE)

# Low ambition
ebv_datacubepaths(ifname)
ebv_trend(filepath = ifname,datacubepath = "scenario_3/metric_6/ebv_cube",entity = 1)
new <- ebvcube::ebv_read(filepath = ifname,datacubepath = "scenario_3/metric_6/ebv_cube",
                         entity = 1,timestep = 2:3,
                         type = "r")
terra::time(new) <- as.POSIXct(c("2015-01-01", "2050-01-01"))
new <- project(new, background)

# Convert to vulnerability
#new <- (new * -1) |> predictor_transform(method = "norm")

out <- data.frame(nexus = "Climate", entrypoint = "Climate adaptation",
                  realm = "marine",
                  scenario = c("low", "low"),
                  year = c(2020, 2050)
) |> dplyr::mutate(
  mean = terra::global(new, "mean", na.rm = TRUE)[,1],
  median =  terra::global(new, median, na.rm = TRUE)[,1],
  min = terra::global(new, "min", na.rm = TRUE)[,1],
  max = terra::global(new, "max", na.rm = TRUE)[,1]
)
write.csv2(out, paste0(path_processed, "climate_climateadapt_mar_low.csv"))

# Calculate relative change
o <- (new[[2]] - new[[1]]) / new[[1]]
writeRaster(o, paste0(path_processed, "climate_climateadapt_mar_low.tif"),overwrite= TRUE)

