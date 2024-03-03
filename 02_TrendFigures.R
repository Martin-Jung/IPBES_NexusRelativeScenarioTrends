library(tidyverse)
library(ggplot2)
library(scales)
library(scico)
library(ggsci)

# Load parameters
source("00_Function.R")
# Relative change function
relChange <- function(v, fac = 100){ 
  if(v[1]==0) return(v * fac)
  ((v - v[1]) / v[1]) * fac 
}

# Create a new theme
theme_lightgrey <- function (base_size = 18, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      plot.background = element_rect(fill = "#ededed", color = NA),
      panel.grid.major  = element_line(color = "white",linewidth = .5),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "#ededed", fill = NA),
      axis.line = element_line(color = "#ededed"),
      axis.ticks = element_line(color = "#ededed"),
      axis.text = element_text(color = "black")
    )
}

# ---- #
# Load the csv files per indicator
csv_files <- list.files(path_processed,full.names = TRUE)
csv_files <- csv_files[has_extension(csv_files, "csv")]

# Load in files
df <- csv_files %>% 
  map_dfr(read_csv2)

# Note starting data was recoded to be at 2020 throughout!
o <- df |> filter(year %in% c(2020,2050)) |> 
  dplyr::group_by(nexus, entrypoint, scenario, realm) |> 
  # Add relative change to the mean
  dplyr::mutate(mean = relChange(mean)) |> dplyr::ungroup() |> 
  # Group and average across realms
  dplyr::group_by(nexus, entrypoint, scenario, year) |> 
  summarise(mean = mean(mean)) |> dplyr::ungroup()

# Relabel indicator names climate
o <- o |> mutate(entrypoint = case_when(entrypoint == "Climate adaptation" ~ "Climate adaptation (risk)",
                                   entrypoint == "Impact" ~ "Climate impacts (risk)",
                                   TRUE ~ entrypoint))

# --------------------------------- #
### Build time series plot ####

gt <- ggplot(data = o, aes(x = year, y = mean, group = scenario, colour = scenario)) + 
  # Theming
  theme_lightgrey(base_size = 18) +
  geom_hline(yintercept = 0, linetype = "dotted",linewidth = .75) +
  geom_line(linewidth = 1.5) +
  # Nature colour scale
  scale_color_npg() +
    guides(colour = guide_legend(title = "Policy ambition")) +
    theme(legend.position = "bottom") +
  facet_wrap(nexus~entrypoint,scales = "free_y") +
  #Remove y-axis labels
  theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
  # Axis labels
  ylab(label = expression("Low values" %->% "High values")) + xlab(label = "")
gt

ggsave(filename = paste0(path_figures, "IndicatorTrends.png"), plot = gt,
       width = 10,height = 9,dpi = 400)
