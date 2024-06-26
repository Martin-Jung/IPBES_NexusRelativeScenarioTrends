library(ggplot2)
library(here)
library(patchwork)
library(scico)

# This script creates a spatial representation of the various extracted layers
# in earlier scripts. The IPBES style guide is followed where epossible
# Source: https://ict.ipbes.net/ipbes-ict-guide/data-and-knowledge-management/technical-guidelines/cartographic-guidelines

# Load in parameters
source('00_Function.R')

# Custom theme with relevant background
theme_mapgrey <- function (base_size = 18, base_family = "") {
  cowplot::theme_map(font_size = base_size, font_family = base_family) %+replace% 
    theme(
      plot.background = element_rect(fill = "#ededed",color = NA),
      panel.grid.major  = element_line(color = "white",linewidth = .5),
      panel.background = element_rect(fill = "#ededed",color = NA),
      panel.border = element_rect(color = NA, fill = NA),
      axis.line = element_line(color = "#ededed"),
      axis.ticks = element_line(color = "#ededed"),
      axis.text = element_text(color = "black")
    )
}

# Figure idea:
# Show at the top a figure of possible nexus interactions in the future
# also next to it a colour map of the number of studies per IPBES region

# -- Display actions as total cumulative rank

# Load all the loaded 
tif_files <- list.files(path_processed,full.names = TRUE)
tif_files <- tif_files[has_extension(tif_files, "tif")]

# Get low ammbition
tif_files_low <- tif_files[grep("_low", tif_files)]
tif_files_high <- tif_files[grep("_high", tif_files)]

assert_that(length(tif_files_high)==length(tif_files_low))

# Separate marine and terrestrial files
tif_files_low_mar <- tif_files_low[grep('mar', tif_files_low)]
tif_files_low_terr <- tif_files_low[grep('ter', tif_files_low)]
tif_files_high_mar <- tif_files_high[grep('mar', tif_files_high)]
tif_files_high_terr <- tif_files_high[grep('ter', tif_files_high)]

# Load in files
ras_low_mar <- terra::rast(tif_files_low_mar)
ras_low_terr <- terra::rast(tif_files_low_terr)
ras_high_mar <- terra::rast(tif_files_high_mar)
ras_high_terr <- terra::rast(tif_files_high_terr)

# Normalize all layers
ras_low_mar <- ras_low_mar |> predictor_transform(option = "norm")
ras_high_mar <- ras_high_mar |> predictor_transform(option = "norm")
ras_low_terr <- ras_low_terr |> predictor_transform(option = "norm")
ras_high_terr <- ras_high_terr |> predictor_transform(option = "norm")

# ------------------ #
#### Combine all and build single relative indicator figure ####

# Normalized score by addinig individual risk and exposure maps
out_low <- terra::app(c(ras_low_mar,ras_low_terr), 'mean',na.rm=TRUE)
out_high <- terra::app(c(ras_high_mar,ras_high_terr), 'mean',na.rm=TRUE)
# out_low <- sum(ras_low,na.rm = TRUE) |> predictor_transform(option = "norm")
# out_high <- sum(ras_high,na.rm = TRUE) |> predictor_transform(option = "norm")

names(out_low) <- "Low Ambition"
names(out_high) <- "High ambition"

# Transform IPBES regions and aggregate
wm <- regions_ipbes |> sf::st_transform(crs = sf::st_crs(background)) |> 
  group_by(IPBES_regi) |> dplyr::summarise()

g_combined <- ggplot() +
  tidyterra::geom_spatraster(data = c(out_low,out_high)) +
    facet_wrap(~lyr,ncol = 1) +
  geom_sf(data = wm, fill = NA, colour = "black", lwd = 1) +
  theme_mapgrey(base_size = 20) +
  scale_fill_gradientn(colours = scico(10, palette = 'lipari',direction = 1),na.value = NA) +
  # scale_fill_gradientn(colours = scico(10, palette = 'glasgow',direction = -1),na.value = NA) +
    guides(fill = guide_colorbar(title = "Normalized\nscore")) +
    theme(legend.position = "bottom",legend.justification = "center",legend.key.width = unit(.8,"in")) + # Previous 1.2 in
  labs(title = "Potential future co-occurring\nnexus interactions by 2050") +
  theme(plot.title = element_text(hjust = 0.5))

g_combined
# Make a figure
ggsave(plot = g_combined, paste0(path_figures, "Fig_MapRanking_overall.png"),
       width = 14,height = 10)
ggsave(plot = g_combined, paste0(path_figures, "Fig_MapRanking_overall.svg"),
       width = 14,height = 10)

# --- #
# Combine with patchwork
library(patchwork)

# gg <- g_combined #/ gt
# ggsave(filename = paste0(path_figures, "Combined.png"),plot = gg,width = 14,height = 20)
#cowplot::plot_grid(g_combined, gt,rel_widths = c(1,2),ncol = 1,labels = "AUTO")

# NOTE: This is done now externally to allow better control of elements!

# g_gaps <- ggplot() +
#   cowplot::theme_map(font_size = 18) +
#   geom_sf(data = ipbes, aes(fill = studycount)) +
#   scale_fill_viridis_c(guide = guide_colorbar(title = "Number of studies")) +
#   labs(title = "Research on\n future nexus interactions") +
#   theme(plot.title = element_text(hjust = 0.5))
# # Make a figure
# ggsave(plot = g_gaps, paste0(path_figures, "Fig_researchgaps.png"),width = 10,height = 8)

# ------------------ #
#### Alternative map visualization ####
# Don't make a cumulative sum but a global hexagonal grid and 
# colour each nexus element respectively
# Source: https://labs.ala.org.au/posts/2024-01-25_hex_point_maps/

# Make a global grid 
hex_grid <- st_make_grid(regions_ipbes |> sf::st_transform(4326),
                         cellsize = 10,
                         # cellsize = 1000000, # Resoluton of hex grid
                         what = "polygons",
                         square = FALSE,
                         flat_topped = TRUE) |> 
  st_as_sf() |> 
  st_set_geometry("hex_geometry") |> 
  tibble::rowid_to_column(var = "hex_id")

# Short plot
ggplot() +
  geom_sf(data = regions_ipbes |> sf::st_transform(4326),
  colour = "darkgrey",
  fill = NA,
  linewidth = 0.3) +
  geom_sf(data = hex_grid, 
          fill = NA, 
          col = "deepskyblue4", 
          linewidth = 0.2) +
  theme_void()

assert_that(length(tif_files_low)>0)

# Extract per hexagon the respective values
ras_biodiversity_low <- terra::rast(tif_files_low[grep("biodiv", tif_files_low)]) |> 
  predictor_transform(option = "norm") |> mean(na.rm=TRUE)
ras_food_low <- terra::rast(tif_files_low[grep("food",tif_files_low)]) |> 
  predictor_transform(option = "norm") |> mean(na.rm=TRUE)
ras_water_low <- terra::rast(tif_files_low[grep("water",tif_files_low)]) |> 
  predictor_transform(option = "norm") |> mean(na.rm=TRUE)
# ras_health_low <- terra::rast(tif_files_low[grep("health",tif_files_low)])
ras_climate_low <- terra::rast(tif_files_low[grep("climate",tif_files_low)]) |> 
  predictor_transform(option = "norm") |> mean(na.rm=TRUE)

# High
ras_biodiversity_high <- terra::rast(tif_files_high[grep("biodiv", tif_files_high)]) |> 
  predictor_transform(option = "norm") |> mean(na.rm=TRUE)
ras_food_high <- terra::rast(tif_files_high[grep("food",tif_files_high)]) |> 
  predictor_transform(option = "norm") |> mean(na.rm=TRUE)
ras_water_high <- terra::rast(tif_files_high[grep("water",tif_files_high)]) |> 
  predictor_transform(option = "norm") |> mean(na.rm=TRUE)
# ras_health_low <- terra::rast(tif_files_low[grep("health",tif_files_low)])
ras_climate_high <- terra::rast(tif_files_high[grep("climate",tif_files_high)]) |> 
  predictor_transform(option = "norm") |> mean(na.rm=TRUE)

out_low <- c(ras_biodiversity_low, ras_food_low, ras_water_low, ras_climate_low)
out_high <- c(ras_biodiversity_high, ras_food_high, ras_water_high, ras_climate_high)
names(out_high) <- names(out_low) <- c("Biodiversity", "Food", "Water", "Climate")

# Make smaller hexagon within
vertex_coords <- hex_grid |> 
  mutate(vertices = pmap(
    .l = list(x = hex_geometry),
    .f = function(x) {
      x |>
        st_buffer(dist = -2.5) |>         # STEP 2: set size of smaller hex
        st_coordinates() |>               # STEP 3: get vertex coordinates of smaller hex        
        as_tibble() |>                    # convert matrix to tibble  
        st_as_sf(coords = c("X", "Y")) |> # convert tibble to simple features
        select(-L1, -L2) |>               # remove unnecessary columns
        mutate(vertex_position = 1:7)     # STEP 4: number vertices 
    })) |> 
  unnest(cols = vertices)

# Get vertex coordinates
vertex_centroid_coords <- vertex_coords |> 
  mutate(geometry = ifelse(vertex_position == 7,      
                           st_centroid(hex_geometry), 
                           geometry)) |> 
  st_drop_geometry()

# Extract values per hexgrid
z1 <- terra::extract(x = out_low, hex_grid,
                    # |> sf::st_transform(crs = terra::crs(background)),
                    fun = "mean", na.rm = TRUE)

z2 <- terra::extract(x = out_high, hex_grid,
                    # |> sf::st_transform(crs = terra::crs(background)),
                    fun = "mean", na.rm = TRUE)

z1 <- exactextractr::exact_extract(out_low, hex_grid |> sf::st_transform(crs = terra::crs(background)),
                                  "mean")
z2 <- exactextractr::exact_extract(out_high, hex_grid |> sf::st_transform(crs = terra::crs(background)),
                                   "mean")
assert_that(nrow(z1) == nrow(hex_grid))
zz <- bind_rows(
  z1 |> 
  # dplyr::rename(hex_id = ID) |> 
  dplyr::mutate(hex_id = 1:nrow(z)) |>
  reshape2::melt(id.vars = "hex_id") |> 
  # Drop missing data
  tidyr::drop_na() |> 
  dplyr::mutate(variable = str_remove(variable, "mean.")) |>
  # Assign corners
  dplyr::mutate(vertex_position = case_when(
    variable == "Biodiversity" ~ 1,
    variable == "Food" ~ 3,
    variable == "Water" ~ 5,
    variable == "Climate" ~ 6
    )) |> 
  dplyr::mutate(ambition = "low"),
  z2 |> 
    # dplyr::rename(hex_id = ID) |> 
    dplyr::mutate(hex_id = 1:nrow(z)) |>
    reshape2::melt(id.vars = "hex_id") |> 
    # Drop missing data
    tidyr::drop_na() |> 
    dplyr::mutate(variable = str_remove(variable, "mean.")) |>
    # Assign corners
    dplyr::mutate(vertex_position = case_when(
      variable == "Biodiversity" ~ 1,
      variable == "Food" ~ 3,
      variable == "Water" ~ 5,
      variable == "Climate" ~ 6
    )) |> 
    dplyr::mutate(ambition = "high")
)

# Join in extracted indicator values
ind_points <- zz |>
  left_join(vertex_centroid_coords,
            by = join_by(vertex_position, hex_id)) |> 
  sf::st_as_sf() |> 
  sf::st_set_crs(value = sf::st_crs(4326))
  # sf::st_set_crs(value = sf::st_crs(terra::crs(background))) |> 
  # sf::st_transform(terra::crs(background))

# Map
wm <- regions_ipbes |> sf::st_transform(crs = sf::st_crs(background)) |> 
  group_by(IPBES_regi) |> dplyr::summarise()

gm1 <- ggplot() +
  geom_sf(data = wm |> st_transform(st_crs(4326)) , fill = NA, colour = "black") +
  geom_sf(data = hex_grid, fill = NA) +
  # Add points
  geom_sf(data = ind_points |> dplyr::filter(ambition == "low"),
          aes(colour = variable, size = value)) +
    scale_color_manual(values = cols) +
  # Guides
  guides(colour = guide_legend(title = ""),
         size = guide_legend(title = "Normalized\nscore")) +
  scale_size_binned_area(max_size = 3) +
  theme_mapgrey(base_size = 20) +
  labs(title = "Potential for future \nnexus interactions by 2050\n(low ambition)") +
  theme(plot.title = element_text(hjust = 0.5))
gm1

gm2 <- ggplot() +
  geom_sf(data = wm |> st_transform(st_crs(4326)) , fill = NA, colour = "black") +
  geom_sf(data = hex_grid, fill = NA) +
  # Add points
  geom_sf(data = ind_points |> dplyr::filter(ambition == "high"),
          aes(colour = variable, size = value)) +
  scale_color_manual(values = cols) +
  # Guides
  guides(colour = guide_legend(title = ""),
         size = guide_legend(title = "Normalized\nscore")) +
  scale_size_binned_area(max_size = 3) +
  theme_mapgrey(base_size = 20) +
  labs(title = "Potential for future \nnexus interactions by 2050\n(high ambition)") +
  theme(plot.title = element_text(hjust = 0.5))
gm2

pg <- cowplot::plot_grid(gm1, gm2+ guides(colour = "none", size = "none"),
                         align = "hv",nrow = 1)

cowplot::ggsave2(plot = pg, paste0(path_figures, "Fig_MapRanking_nexus.png"),
       width = 20,height = 18)
