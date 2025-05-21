library(ggplot2)
library(here)
library(patchwork)
library(scico)

# This script roughly estimates the number of layers spatially 
# and overall per nexus elements

# Load in parameters
source('00_Function.R')

df <- readxl::read_xlsx("C:/Users/tuete/United Nations/Teamkanal - Chapter 3/Explorative Figure/FigureData.xlsx")

# Maximum number of indicators across nexus element section and realm
max_scenarios <- dplyr::n_distinct(df$Entrypoint)*3 # 3 sections per element except for health

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

#### Full map of scenario data gaps ####

# Load all the loaded 
tif_files <- list.files(path_processed,full.names = TRUE)
tif_files <- tif_files[has_extension(tif_files, "tif")]
# Get low ambition as example
tif_files_low <- tif_files[grep("_low", tif_files)]

# Load terr and mar
ras_ter <- rast(tif_files_low[grep('terr', tif_files_low)])
ras_mar <- rast(tif_files_low[grep('mar', tif_files_low)])
names(ras_ter) <- tools::file_path_sans_ext(basename(tif_files_low[grep('terr', tif_files_low)]))
names(ras_mar) <- tools::file_path_sans_ext(basename(tif_files_low[grep('mar', tif_files_low)]))

# Summarize for terrestrial
n1 <- terra::not.na(ras_ter) |> sum()
n1 <- terra::mask(n1, background)

# For marine
n2 <- terra::not.na(ras_mar) |> sum()
n2[n2==0] <- NA
n2 <- terra::mask(n2, background,inverse = TRUE)

nn <- sum(n1, n2, na.rm = TRUE)
nn <- terra::mask(nn, sf::st_bbox(background) |> sf::st_as_sfc() |> sf::st_as_sf())
nn <- (max_scenarios - nn) / max_scenarios
names(nn) <- "Scenario data gaps"
  
g_gaps <- ggplot() +
  tidyterra::geom_spatraster(data = nn) +
  facet_wrap(~lyr) +
  # geom_sf(data = wm, fill = NA, colour = "black", lwd = 1) +
  theme_mapgrey(base_size = 20) +
  # scale_fill_gradientn(colours = scico(10, palette = 'lipari',direction = -1),na.value = "#ededed") +
  scale_fill_gradientn(colours = scico(10, palette = 'acton',direction = -1),na.value = "#ededed") +
  guides(fill = guide_colorbar(title = "Few                 Many")) +
  theme(legend.position = "right",legend.title.position = 'right',
        legend.ticks = element_blank(), legend.text = element_blank(),
        legend.title = element_text(size = 16,vjust = 1, angle=90), legend.justification = "center") +
  labs(title = "") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))
g_gaps
# Save the figure outputs as png
ggsave(plot = g_gaps, paste0(path_figures, "Fig_MapGaps_overall.png"),
       width = 8,height = 6, dpi = 400)
ggsave(plot = g_gaps, paste0(path_figures, "Fig_MapGaps_overall.eps"),
       width = 8,height = 6, dpi = 400)
ggsave(plot = g_gaps, paste0(path_figures, "Fig_MapGaps_overall.svg"),
       width = 8,height = 6, dpi = 400)

# -------------- #
#### Overview plot per nexus element ####
# Small stacked barplot that highlights the number of assessed indicators
# per scenario
library(tidyr)
library(dplyr)

# Load file and expand missing columns
df <- readxl::read_xlsx("C:/Users/tuete/United Nations/Teamkanal - Chapter 3/Explorative Figure/FigureData.xlsx") |> 
  # dplyr::filter(!is.na(Indicator), !is.na(Scenario)) |> 
  dplyr::select(Nexus:Realm, Scenario, Source)

# Add Covered variable if present
df <- df |> dplyr::mutate(covered = if_else(is.na(Scenario) , 0, 1) ) |>
  dplyr::select(Nexus:Realm, covered) |> distinct()

df$Nexus <- factor(df$Nexus, levels = names(cols))

df |> group_by(Nexus) |> summarise(co = sum(covered))

# Make a stacked barplot
gb <- ggplot(df |> filter(covered == 1), aes(x = Nexus, group = Entrypoint, fill = Nexus)) +
  theme_bw(base_size = 18) +
  geom_bar() +
  scale_fill_manual(values = cols)

gb

#### Experimental Correlation plot of indicators ####
library(ggcorrplot)

co <- terra::layerCor(ras_ter, "cor")

ggcorrplot(co$correlation, method = "square",
           hc.order = TRUE, type = "upper",outline.color = 'black')
