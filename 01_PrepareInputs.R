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
                  realm = "terrestrial", indicator = "Terrestrial protected areas",
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

# NOTE: Marine indicators not summarized since trends will be similar as for terrestrial.

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
                  realm = "terrestrial",indicator = "Biotic intactness",
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
                  realm = "terrestrial",indicator = "Biotic intactness",
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
                          realm = "terrestrial", indicator = "Total potential water supply",
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
                          realm = "terrestrial",indicator = "Total potential water supply",
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
                          realm = "terrestrial",indicator = "Total actual water withdrawal",
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
                          realm = "terrestrial",indicator = "Total actual water withdrawal",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "water_demand_terr_low.csv"))

#### Health - Heat mortality ----
h <- read_ncdf("raw_data/isimip/trm-tsukuba_gfdl-esm2m_ewembi_rcp26_2005soc_an-tot-heat-all_global_annual_2006_2099.nc4")
l <- read_ncdf("raw_data/isimip/trm-tsukuba_gfdl-esm2m_ewembi_rcp60_2005soc_an-tot-heat-all_global_annual_2006_2099.nc4")

# Get reference year 2015
ref <- ibis.iSDM:::stars_to_raster(h, which= 10)[[1]]
ref <- terra::project(ref, background)
ref <- terra::mask(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(h, which= 45)[[1]]
new <- terra::project(new, background)
new <- terra::mask(new, background)

# Relative change
out <- (new - ref)# / ref
writeRaster(out, paste0(path_processed, "health_heatmort_terr_high.tif"), overwrite= TRUE)

# Summarize 
df <- as.data.frame(h) |> tidyr::drop_na()
df$year <- lubridate::year(df$time)
df <- df |> dplyr::filter(an.tot.heat.all>0) |> 
  group_by(year) |> summarise(sum = sum(an.tot.heat.all),
                                        min = min(an.tot.heat.all), max = max(an.tot.heat.all),
                                        mean = mean(an.tot.heat.all), median = median(an.tot.heat.all)) |> 
  # Use the sum in this case as we are dealing with a count
  # Overwrite mean as this is used by default in future processing.
  mutate(mean = sum)
gc()
df <- df |> dplyr::mutate(nexus = "Health", entrypoint = "Green and blue space",
                          realm = "terrestrial",indicator = "Heat mortality",
                          scenario = "high")
write.csv2(df, paste0(path_processed, "health_heatmort_terr_high.csv"))

# Low ambition
# Get reference year 2015
ref <- ibis.iSDM:::stars_to_raster(l, which= 10)[[1]]
ref <- terra::project(ref, background)
ref <- terra::mask(ref, background)

# Get 2050
new <- ibis.iSDM:::stars_to_raster(l, which= 45)[[1]]
new <- terra::project(new, background)
new <- terra::mask(new, background)

# Relative change
out <- (new - ref)# / ref
writeRaster(out, paste0(path_processed, "health_heatmort_terr_low.tif"), overwrite= TRUE)

# Summarize 
df <- as.data.frame(l) |> tidyr::drop_na()
df$year <- lubridate::year(df$time)
df <- df |> dplyr::filter(an.tot.heat.all>0) |> 
  group_by(year) |> summarise(sum = sum(an.tot.heat.all),
                              min = min(an.tot.heat.all), max = max(an.tot.heat.all),
                              mean = mean(an.tot.heat.all), median = median(an.tot.heat.all)) |> 
  # Use the sum in this case as we are dealing with a count
  # Overwrite mean as this is used by default in future processing.
  mutate(mean = sum)
gc()
df <- df |> dplyr::mutate(nexus = "Health", entrypoint = "Green and blue space",
                          realm = "terrestrial",indicator = "Heat mortality",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "health_heatmort_terr_low.csv"))

#### Health - Infectious diseases ----
# Load for Malaria
# LTS Length of the transmission season for historical and projected period
mal_h <- read_ncdf("raw_data/OSF_Health/Malaria outputs/LCMI/lcmi_ensmean_historical_lts_1970_1999.nc4") 
mal_f <- read_ncdf("raw_data/OSF_Health/Malaria outputs/LCMI/lcmi_ensmean_rcp26_ssp1soc_lts_2040_2069.nc4")
deng_h <- read_ncdf("raw_data/OSF_Health/Dengue outputs/DGM/dgm_ensmean_historical_lts_1970_1999.nc4")
deng_f <- read_ncdf("raw_data/OSF_Health/Dengue outputs/DGM/dgm_ensmean_rcp26_ssp1soc_lts_2040_2069.nc4")

# Get reference year
ref1 <- ibis.iSDM:::stars_to_raster(mal_h, which = 1)[[1]]
ref1 <- terra::project(ref1, background)
ref2 <- ibis.iSDM:::stars_to_raster(deng_h, which = 1)[[1]]
ref2 <- terra::project(ref2, background)

# Get 2050
new1 <- ibis.iSDM:::stars_to_raster(mal_f, which= 1)[[1]]
new1 <- terra::project(new1, background)
new2 <- ibis.iSDM:::stars_to_raster(deng_f, which= 1)[[1]]
new2 <- terra::project(new2, background)

# Relative change
out1 <- (new1 - ref1)
out2 <- (new2 - ref2)
out <- mean(out1,out2)
writeRaster(out, paste0(path_processed, "health_zoo_terr_high.tif"), overwrite = TRUE)

df <- rbind(
  as.data.frame(mean(ref1,ref2)) |> tidyr::drop_na() |> mutate(year = 2020),
  as.data.frame(mean(new1,new2)) |> tidyr::drop_na() |> mutate(year = 2050)
)
df <- df |> group_by(year) |> summarise(sum = sum(lts),
                                        min = min(lts), max = max(lts),
                                        mean = mean(lts), median = median(lts))
df <- df |> dplyr::mutate(nexus = "Health", entrypoint = "Vector-borne diseases",
                          realm = "terrestrial", indicator = "Length of \ninfectious transmission season",
                          scenario = "high")
write.csv2(df, paste0(path_processed, "health_zoo_terr_high.csv"))

# Low ambition
mal_f <- read_ncdf("raw_data/OSF_Health/Malaria outputs/LCMI/lcmi_ensmean_rcp85_ssp5soc_lts_2040_2069.nc4")
deng_f <- read_ncdf("raw_data/OSF_Health/Dengue outputs/DGM/dgm_ensmean_rcp85_ssp5soc_lts_2040_2069.nc4")

# Get reference year
ref1 <- ibis.iSDM:::stars_to_raster(mal_h, which = 1)[[1]]
ref1 <- terra::project(ref1, background)
ref2 <- ibis.iSDM:::stars_to_raster(deng_h, which = 1)[[1]]
ref2 <- terra::project(ref2, background)

# Get 2050
new1 <- ibis.iSDM:::stars_to_raster(mal_f, which= 1)[[1]]
new1 <- terra::project(new1, background)
new2 <- ibis.iSDM:::stars_to_raster(deng_f, which= 1)[[1]]
new2 <- terra::project(new2, background)

# Relative change
out1 <- (new1 - ref1)
out2 <- (new2 - ref2)
out <- mean(out1,out2)
writeRaster(out, paste0(path_processed, "health_zoo_terr_low.tif"), overwrite = TRUE)

df <- rbind(
  as.data.frame(mean(ref1,ref2)) |> tidyr::drop_na() |> mutate(year = 2020),
  as.data.frame(mean(new1,new2)) |> tidyr::drop_na() |> mutate(year = 2050)
)
df <- df |> group_by(year) |> summarise(sum = sum(lts),
                                        min = min(lts), max = max(lts),
                                        mean = mean(lts), median = median(lts))
df <- df |> dplyr::mutate(nexus = "Health", entrypoint = "Vector-borne diseases",
                          realm = "terrestrial", indicator = "Length of \ninfectious transmission season",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "health_zoo_terr_low.csv"))

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
                          realm = "terrestrial", indicator = "Irrigated summer wheat yield",
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
                          realm = "terrestrial", indicator = "Irrigated summer wheat yield",
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
                          realm = "marine",indicator = "Total fish catch",
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
                          realm = "marine",indicator = "Total fish catch",
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
                          realm = "terrestrial",indicator = "Multi-sectoral climate risk",
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
                          realm = "terrestrial",indicator = "Multi-sectoral climate risk",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "climate_impact_terr_low.csv"))

# --------- #
# Marine
# Using median projections of the pH from here https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0259391/nc/median/
# High ambition
ph_ref <- read_ncdf("raw_data/Marine_Accidification/pHT_median_historical.nc")
ph_high <- read_ncdf("raw_data/Marine_Accidification/pHT_median_ssp126.nc")
#  West: -179.5 South: -89.5 East: 179.5 North: 89.5
sf::st_crs(ph_ref) <- sf::st_crs(4326)
sf::st_crs(ph_high) <- sf::st_crs(4326)
# HACKY : change dimensions manually
dim <- stars::st_dimensions(ph_ref)
dim$lon$from <- -179.5
dim$lon$to <- 179.5
dim$lat$from <- -89.5
dim$lat$to <- 89.5
st_dimensions(ph_ref) <- dim
dim <- stars::st_dimensions(ph_high)
dim$lon$from <- -179.5
dim$lon$to <- 179.5
dim$lat$from <- -89.5
dim$lat$to <- 89.5
st_dimensions(ph_high) <- dim

# Reference
ref <- ibis.iSDM:::stars_to_raster(ph_ref, which = 18)[[1]]
ref <- shift(ref,dx=200) # Shift to 0 meridian
ref <- terra::project(ref, background)

# Get 2050
# seq(2020,2100,10)
new <- ibis.iSDM:::stars_to_raster(ph_high, which= 4)[[1]]
new <- shift(new, dx=200) # Shift to 0 meridian
new <- terra::project(new, background)

# Relative change
# Inverse to highlight impacts of accidification
out <- ((new - ref) / ref) * -1
writeRaster(out, paste0(path_processed, "climate_impact_mar_high.tif"), overwrite= TRUE)

# Get projection
df <- as.data.frame(ph_high) |> tidyr::drop_na()
df$year <- df$time+5
df$pHT <- units::drop_units(df$pHT)
df <- df |> group_by(year) |> summarise(sum = sum(pHT),
                                        min = min(pHT), max = max(pHT),
                                        mean = mean(pHT), median = median(pHT))
gc()
df <- df |> dplyr::mutate(nexus = "Climate", entrypoint = "Impact",
                          realm = "marine", indicator = "Ocean pH values",
                          scenario = "high")
write.csv2(df, paste0(path_processed, "climate_impact_mar_high.csv"))

# --- #
# Low ambition
ph_low <- read_ncdf("raw_data/Marine_Accidification/pHT_median_ssp585.nc")
sf::st_crs(ph_low) <- sf::st_crs(4326)
# HACKY : change dimensions manually
dim <- stars::st_dimensions(ph_low)
dim$lon$from <- -179.5
dim$lon$to <- 179.5
dim$lat$from <- -89.5
dim$lat$to <- 89.5
st_dimensions(ph_low) <- dim

# Get 2050
# seq(2020,2100,10)
new <- ibis.iSDM:::stars_to_raster(ph_low, which= 4)[[1]]
new <- shift(new, dx=200) # Shift to 0 meridian
new <- terra::project(new, background)
# Relative change
# Inverse to highlight impacts of accidification
out <- ((new - ref) / ref) * -1
writeRaster(out, paste0(path_processed, "climate_impact_mar_low.tif"), overwrite= TRUE)

df <- as.data.frame(ph_low) |> tidyr::drop_na()
df$year <- df$time+5
df$pHT <- units::drop_units(df$pHT)
df <- df |> group_by(year) |> summarise(sum = sum(pHT),
                                        min = min(pHT), max = max(pHT),
                                        mean = mean(pHT), median = median(pHT))
gc()
df <- df |> dplyr::mutate(nexus = "Climate", entrypoint = "Impact",
                          realm = "marine",indicator = "Ocean pH values",
                          scenario = "low")
write.csv2(df, paste0(path_processed, "climate_impact_mar_low.csv"))

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
                  realm = "terrestrial",indicator = "Bioenergy production",
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
                  realm = "terrestrial",indicator = "Bioenergy production",
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
                  realm = "marine", indicator = "Coastal flooding risk",
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
                  realm = "marine", indicator = "Coastal flooding risk",
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

