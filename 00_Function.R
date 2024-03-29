# Generic packages and paths
# --------------- #
#### Packages ####
library(terra)
library(sf)
library(gdalUtilities)
library(ebvcube)
library(stars)
library(tidyterra)
library(exactextractr)
library(assertthat)
library(ibis.iSDM)
library(lubridate)
# --- #

# Path raw data
path_rawdata <- "raw_data/"
dir.create(path_rawdata,showWarnings = FALSE)

# Path data
path_data <- "data/"
dir.create(path_data,showWarnings = FALSE)

# Path processed data
path_processed <- "processed_data/"
dir.create(path_processed,showWarnings = FALSE)

# path figures
path_figures <- "figures/"
dir.create(path_figures, showWarnings = FALSE)

# Load background raster
background <- terra::rast("data/globalgrid_mollweide_50km.tif")
background[background>=0] <- 1
# IPBES uses the global robinson projection by default
proj_robin <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
background <- terra::project(background, proj_robin)

# IPBES regions shapefile (origin bTC analysis) reprojected to robinson 
regions_ipbes <- sf::st_read("raw_data/BendingTheCurve_supportingMaterial_final_28August2020/Analysis_R3_final/simplifiedIPBESsubregions_shapefile/IPBES_sr_landOnly_NoExcluded_simplified_withJoinOnAttributes_reProjectedWinkelTripel.shp",
                             quiet = TRUE) |> 
  dplyr::select(IPBES_regi, IPBES_sub) |> 
  sf::st_transform(crs = sf::st_crs(background))

# Generic colours
cols <- c("None" = "grey30", "Biodiversity" = "#C6D68A", "Food" = "#B65719",
          "Water" = "#4A928F", "Health" = "#791E32", "Climate" = "#BAB0C9")

# ---------- #
#### Checks ####
assertthat::assert_that(
  dir.exists(path_rawdata),
  dir.exists(path_processed),
  dir.exists(path_figures),
  is.Raster(background),
  inherits(regions_ipbes, 'sf')
)