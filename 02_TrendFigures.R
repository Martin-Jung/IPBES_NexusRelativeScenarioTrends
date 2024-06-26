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
o <- df |> 
  filter(year >= 2020 & year <= 2050) |>
  # filter(year %in% c(2020,2050)) |>
  dplyr::group_by(nexus, entrypoint, indicator, scenario, realm) |> 
  # Add relative change to the mean
  dplyr::mutate(mean = relChange(mean)) |> dplyr::ungroup() |> 
  # Group and average across realms
  dplyr::group_by(nexus, entrypoint, indicator, scenario, year) |> 
  summarise(mean = mean(mean)) |> dplyr::ungroup()

# Relabel indicator names climate
o <- o |> mutate(entrypoint = case_when(entrypoint == "Climate adaptation" ~ "Climate adaptation (risk)",
                                   entrypoint == "Impact" ~ "Climate impacts (risk)",
                                   TRUE ~ entrypoint))

# Rename nexus climate to climate change
o$nexus[o$nexus=="Climate"] <- "Climate change"

# --------------------------------- #
### Build time series plot ####
# Split label
o$indicator[o$indicator=="Length of \r\ninfectious transmission season"] <- "Length of infectious \ntransmission season"
o$indicator[o$indicator=="Total actual water withdrawal"] <- "Total actual \nwater withdrawal"
o$indicator[o$indicator=="Total potential water supply"] <- "Total potential \nwater supply"

gt <- ggplot(data = o, aes(x = year, y = mean, group = scenario, colour = scenario)) + 
  # Theming
  theme_lightgrey(base_size = 18) +
  # Add filled background grid
    # geom_rect(aes(fill = entrypoint), alpha = 0.05,show.legend = FALSE,
              # ymin = -Inf, ymax = Inf, xmin = -Inf, xmax= Inf) +
  geom_hline(yintercept = 0, linetype = "dotted",linewidth = .75) +
  geom_line(linewidth = 1.5) +
  # Nature colour scale
  scale_color_npg() +
    guides(colour = guide_legend(title = "Policy ambition")) +
    theme(legend.position = "bottom") +
  # facet_wrap(nexus~entrypoint,scales = "free_y",ncol = 4) +
  facet_wrap(nexus~indicator,scales = "free",ncol = 5) +
    theme(strip.text = element_text(size = 18)) +
  #Remove y-axis labels
  theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank(),
        axis.text.x.bottom = element_text(size = 14)) +
  # theme(strip.switch.pad.grid = unit(0.1,"cm")) +
  # Axis labels
  ylab(label = expression("Low values" %->% "High values")) + xlab(label = "")

ggsave(filename = paste0(path_figures, "IndicatorTrends.png"), plot = gt,
       width = 16,height = 14,dpi = 400)

# --------------------------------- #
### Add Full indicator breakdown ####
# Figure idea:
# Don't aggregate across Realms and indicators
# Instead showing simply the change relative to baseline.

# Get actual indicator names
inds <- readxl::read_xlsx("C:/Users/tuete/United Nations/Teamkanal - Chapter 3/Explorative Figure/FigureData.xlsx") |> 
  dplyr::filter(!is.na(Indicator), !is.na(Scenario)) |>
  dplyr::select(Nexus:Realm, Scenario.Group, Indicator) |> distinct()
names(inds) <- tolower(names(inds))
# Recode scenarios 
inds$scenario <- ifelse(inds$scenario.group == "Low ambition", "low", "high")
inds <- inds |> dplyr::select(-scenario.group)

# Format data for barplot
o <- df |> filter(year %in% c(2020,2050)) |> 
  dplyr::group_by(nexus, entrypoint, scenario, realm) |> 
  # Add relative change to the mean
  dplyr::mutate(mean = relChange(mean)) |> dplyr::ungroup() |> 
  dplyr::select(nexus, entrypoint, scenario, realm, year, mean) |> 
  # Convert to uppercase
  dplyr::mutate(realm = stringr::str_to_title(realm)) |> distinct()
# Add indicator name to it
o <- o |> dplyr::left_join(inds)

# Relabel indicator names climate
o <- o |> mutate(entrypoint = case_when(entrypoint == "Climate adaptation" ~ "Climate adaptation (risk)",
                                        entrypoint == "Impact" ~ "Climate impacts (risk)",
                                        TRUE ~ entrypoint))


gb <- ggplot(data = o |> dplyr::filter(year == 2050),
             aes(x = entrypoint, y = mean, group = scenario, fill = scenario)) + 
  # Theming
  theme_lightgrey(base_size = 18) +
  geom_hline(yintercept = 0, linetype = "dotted",linewidth = .75) +
  geom_point() +
  geom_line(linewidth = 1.5) +
  # Colours set
  scale_fill_manual(values = cols) +
  facet_wrap(~realm)

  guides(colour = guide_legend(title = "Policy ambition")) +
  theme(legend.position = "bottom") +
  facet_wrap(nexus~entrypoint,scales = "free_y",ncol = 4) +
  #Remove y-axis labels
  theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
  # Axis labels
  ylab(label = expression("Low values" %->% "High values")) + xlab(label = "")
gb

ggsave(filename = paste0(path_figures, "IndicatorTrends.png"), plot = gt,
       width = 10,height = 9,dpi = 400)

