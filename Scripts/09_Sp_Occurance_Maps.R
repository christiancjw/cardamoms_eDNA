# ----------------------------------------------------------------
# 08_species_occurrence_map.R
#
# Plots occurrence of selected species on the rotated river map.
# Change species_to_plot to show different species.
# Map rotation matches 05_alpha_diversity.R (60 degrees clockwise).
#
# Run after 07_alpha_diversity.R — rivers, ctr, angle in memory.
# ----------------------------------------------------------------

library(sf)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(tidyr)


# CONFIG — change species here ----------------------------------------
# Single species: one name
# Multiple species: add more names to the vector

species_to_plot <- c(
  "Poropuntius normani",
   "Schistura sp1",
   "Homaloptera confuzona",
   "Schistura sp2",
  "Glyptothorax laosensis",
  "Channa gachua"
)


# 1: Build occurrence data ----------------------------------------
# pa and meta are in memory from 05_alpha_diversity.R
# If not, load them:
# pa   <- read_csv("Data/Clean/fish_community_pa.csv",   show_col_types = FALSE)
# meta <- read_csv("Data/Clean/fish_community_meta.csv", show_col_types = FALSE)

# Pull presence/absence for selected species per event
occ <- pa %>%
  select(Event, any_of(species_to_plot)) %>%
  pivot_longer(-Event, names_to = "Species", values_to = "presence") %>%
  left_join(meta %>% select(Event, Site, Strahler_f, Season), by = "Event") %>%
  group_by(Site, Species) %>%
  summarise(
    detected    = any(presence == 1),   # TRUE if detected in any event at this site
    n_events    = sum(presence == 1),   # number of events detected
    pct_events  = round(100 * mean(presence), 0),
    .groups = "drop"
  ) %>%
  left_join(meta %>% distinct(Site, Long, Lat, Code, Strahler_f), by = "Site") %>%
  filter(!is.na(Long))

occ  # check before plotting


# 2: Convert and rotate site coordinates ----------------------------------------
# Uses same angle and ctr as the river rotation in 05_alpha_diversity.R

occ_sf <- occ %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

geom_occ     <- st_geometry(occ_sf)
geom_occ_rot <- (geom_occ - ctr) * matrix(c(cos(angle), sin(angle),
                                            -sin(angle), cos(angle)), 2, 2) + ctr
st_geometry(occ_sf) <- geom_occ_rot
st_crs(occ_sf) <- 4326


# 3: Plot — single species ----------------------------------------
# Clean version — detected sites filled, absent sites open

if (length(species_to_plot) == 1) {
  
  fig_occ <- ggplot() +
    geom_sf(data = rivers, colour = "grey60", linewidth = 0.4, alpha = 0.8) +
    # Absent sites — open grey circles
    geom_sf(data   = occ_sf %>% filter(!detected),
            shape  = 21, fill = "white", colour = "grey60",
            size   = 4, stroke = 0.8, alpha = 0.7) +
    # Present sites — filled, sized by % events detected
    geom_sf(data   = occ_sf %>% filter(detected),
            aes(size = pct_events, fill = Strahler_f),
            shape  = 21, colour = "black", stroke = 0.8, alpha = 0.95) +
    geom_sf_text(data     = occ_sf,
                 aes(label = Code),
                 size     = 2.8, fontface = "bold", colour = "grey20",
                 nudge_y  = 0.018) +
    scale_fill_manual(values = strahler_f_colors, name = "Strahler order") +
    scale_size_continuous(range = c(3, 10), name = "% events detected") +
    labs(title = species_to_plot) +
    theme_void(base_size = 11) +
    theme(
      plot.title       = element_text(face = "italic", size = 12, hjust = 0.5),
      legend.position  = "right",
      legend.direction = "vertical"
    )
  
  fig_occ
  
}


# 4: Plot — multiple species (faceted) ----------------------------------------
# One panel per species, same map layout

if (length(species_to_plot) > 1) {
  
  fig_occ_multi <- ggplot() +
    geom_sf(data = rivers, colour = "grey60", linewidth = 0.3, alpha = 0.7) +
    # Absent
    geom_sf(data   = occ_sf %>% filter(!detected),
            shape  = 21, fill = "white", colour = "grey60",
            size   = 3, stroke = 0.6, alpha = 0.6) +
    # Present
    geom_sf(data   = occ_sf %>% filter(detected),
            aes(fill = Strahler_f, size = pct_events),
            shape  = 21, colour = "black", stroke = 0.7, alpha = 0.95) +
    geom_sf_text(data     = occ_sf,
                 aes(label = Code),
                 size     = 2.2, fontface = "bold", colour = "grey20",
                 nudge_y  = 0.018) +
    scale_fill_manual(values = strahler_f_colors, name = "Strahler order") +
    scale_size_continuous(range = c(2, 8), name = "% events detected") +
    facet_wrap(~Species, ncol = 2) +  # change ncol to adjust layout
    theme_void(base_size = 10) +
    theme(
      strip.text       = element_text(face = "italic", size = 9),
      legend.position  = "bottom",
      legend.direction = "horizontal"
    )
  
  fig_occ_multi
  
}

