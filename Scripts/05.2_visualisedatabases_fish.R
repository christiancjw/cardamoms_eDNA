# =============================================================================
# 03_visualise_comparison.R
#
# Produces publication-quality figures comparing MIDORI-only vs combined
# database performance.
#
# FIGURES:
#   Fig 1: Stacked bar — assignment resolution by database
#   Fig 2: Bar — OTU assignment changes (improved/same/degraded)
#   Fig 3: Dot plot — species gained from local library
#
# INPUTS:  Data/Clean/ outputs from 02_compare_databases.R
# OUTPUTS: Figures/db_comparison/
# =============================================================================

# ----------------------------------------
# Setup
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(forcats)

input_dir   <- "Data/Clean/"
figures_dir <- "Figures/db_comparison/"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Consistent colours
col_midori   <- "#4E79A7"
col_combined <- "#F28E2B"

resolution_colours <- c(
  "Species"    = "#2C7BB6",
  "Genus"      = "#74C476",
  "Family"     = "#FD8D3C",
  "Order"      = "#C6DBEF",
  "Unresolved" = "#D9D9D9"
)

resolution_levels <- c("Species", "Genus", "Family", "Order", "Unresolved")


# 1: Load data ----------------------------------------
mid     <- read_csv(file.path(input_dir, "fish_midori.csv"),          show_col_types = FALSE)
comb    <- read_csv(file.path(input_dir, "fish_combined.csv"),         show_col_types = FALSE)
compare <- read_csv(file.path(input_dir, "db_comparison_otus.csv"),    show_col_types = FALSE)
gained  <- read_csv(file.path(input_dir, "db_species_gained.csv"),     show_col_types = FALSE)
summary <- read_csv(file.path(input_dir, "db_comparison_summary.csv"), show_col_types = FALSE)


# 2: Figure 1 — Resolution stacked bar ----------------------------------------

res_data <- bind_rows(
  mid  %>% dplyr::count(resolution) %>% mutate(Database = "MIDORI only"),
  comb %>% dplyr::count(resolution) %>% mutate(Database = "MIDORI + Local")
) %>%
  mutate(
    resolution = factor(resolution, levels = resolution_levels),
    Database   = factor(Database, levels = c("MIDORI only", "MIDORI + Local"))
  )

fig1 <- ggplot(res_data, aes(x = Database, y = n, fill = resolution)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = ifelse(n > 1, n, "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = resolution_colours, name = "Assignment level") +
  labs(
    title = "Taxonomic assignment resolution by database",
    x     = NULL,
    y     = "Number of OTUs"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

fig1

ggsave(file.path(figures_dir, "fig1_resolution.pdf"), fig1, width = 6, height = 5)
ggsave(file.path(figures_dir, "fig1_resolution.png"), fig1, width = 6, height = 5, dpi = 300)


# 3: Figure 2 — OTU changes bar ----------------------------------------

change_colours <- c(
  "Improved"           = "#59A14F",
  "Same"               = "#BDBDBD",
  "Degraded"           = "#E15759",
  "Gained by combined" = "#F28E2B",
  "Lost in combined"   = "#4E79A7"
)

change_data <- compare %>%
  dplyr::count(change) %>%
  filter(!is.na(change)) %>%
  mutate(change = fct_reorder(change, n))

fig2 <- ggplot(change_data, aes(x = change, y = n, fill = change)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.3, size = 4) +
  scale_fill_manual(values = change_colours) +
  coord_flip() +
  labs(
    title = "Per-OTU assignment change: MIDORI only → MIDORI + Local",
    x     = NULL,
    y     = "Number of OTUs"
  ) +
  theme_classic(base_size = 12) +
  expand_limits(y = max(change_data$n) * 1.2)

fig2

ggsave(file.path(figures_dir, "fig2_otu_changes.pdf"), fig2, width = 7, height = 4)
ggsave(file.path(figures_dir, "fig2_otu_changes.png"), fig2, width = 7, height = 4, dpi = 300)


# 4: Figure 3 — Species gained dot plot ----------------------------------------
# Only run if there are species gained

nrow(gained)  # check before plotting

if (nrow(gained) > 0) {
  
  # Join similarity from combined table for the gained species
  gained_sim <- comb %>%
    filter(Species %in% gained$Species) %>%
    group_by(Species, Family) %>%
    summarise(max_similarity = max(Similarity, na.rm = TRUE), .groups = "drop")
  
  fig3 <- ggplot(gained_sim,
                 aes(x = max_similarity,
                     y = fct_reorder(Species, max_similarity),
                     colour = Family)) +
    geom_point(size = 4) +
    geom_vline(xintercept = 97, linetype = "dashed", colour = "grey60") +
    annotate("text", x = 97.2, y = 0.7, label = "97% threshold",
             hjust = 0, colour = "grey50", size = 3.2) +
    labs(
      title    = "Species detected only by combined database",
      subtitle = paste(nrow(gained), "species recovered via local reference library"),
      x        = "Best BLAST similarity (%)",
      y        = NULL,
      colour   = "Family"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right") +
    xlim(95, 101)
  
  fig3
  
  ggsave(file.path(figures_dir, "fig3_species_gained.pdf"),
         fig3, width = 9, height = max(4, nrow(gained) * 0.4 + 2))
  ggsave(file.path(figures_dir, "fig3_species_gained.png"),
         fig3, width = 9, height = max(4, nrow(gained) * 0.4 + 2), dpi = 300)
  
} else {
  cat("No species gained — fig3 skipped.\n")
}


# 5: Print summary table ----------------------------------------
print(summary)

cat("\nAll figures saved to", figures_dir, "\n")
cat("Next step: run 04_prepare_ttt2.R\n")
