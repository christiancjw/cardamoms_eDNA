# =============================================================================
# 02_compare_databases.R
#
# Compares OTU recovery and taxonomic assignment between the MIDORI-only
# and combined (MIDORI + local barcodes) database BLAST runs.
#
# INPUTS:  Data/Clean/fish_midori.csv and fish_combined.csv
# OUTPUTS: Data/Clean/db_comparison_otus.csv
#          Data/Clean/db_comparison_summary.csv
#          Data/Clean/db_species_gained.csv
# =============================================================================

# ----------------------------------------
# Setup
library(dplyr)
library(readr)
library(tidyr)

input_dir  <- "Data/Clean/"
output_dir <- "Data/Clean/"


# 1: Load fish tables ----------------------------------------
mid  <- read_csv(file.path(input_dir, "fish_midori.csv"),   show_col_types = FALSE)
comb <- read_csv(file.path(input_dir, "fish_combined.csv"), show_col_types = FALSE)

nrow(mid)   # MIDORI fish OTUs
nrow(comb)  # Combined fish OTUs


# 2: Add resolution level to each table ----------------------------------------
# What is the finest taxonomic rank that was assigned for each OTU?

add_resolution <- function(df) {
  df %>%
    mutate(resolution = case_when(
      !is.na(Species) & !Species %in% c("NoMatch", "Missing-Accession", "") ~ "Species",
      !is.na(Genus)   & !Genus   %in% c("NoMatch", "Missing-Accession", "") ~ "Genus",
      !is.na(Family)  & !Family  %in% c("NoMatch", "Missing-Accession", "") ~ "Family",
      !is.na(Order)   & !Order   %in% c("NoMatch", "Missing-Accession", "") ~ "Order",
      TRUE ~ "Unresolved"
    ))
}

mid  <- add_resolution(mid)
comb <- add_resolution(comb)

# Check resolution breakdown for each database
mid  %>% dplyr::count(resolution, sort = TRUE)
comb %>% dplyr::count(resolution, sort = TRUE)


# 3: OTU-level comparison ----------------------------------------
# Join both tables on OTU_ID to compare assignment per OTU

otu_compare <- full_join(
  mid  %>% select(OTU_ID, Species_mid  = Species,  Genus_mid  = Genus,
                  Family_mid  = Family,  Similarity_mid  = Similarity,
                  resolution_mid  = resolution),
  comb %>% select(OTU_ID, Species_comb = Species,  Genus_comb = Genus,
                  Family_comb = Family,  Similarity_comb = Similarity,
                  resolution_comb = resolution),
  by = "OTU_ID"
)

# Classify what changed between databases for each OTU
resolution_levels <- c("Unresolved", "Order", "Family", "Genus", "Species")

otu_compare <- otu_compare %>%
  mutate(
    rank_mid  = match(resolution_mid,  resolution_levels),
    rank_comb = match(resolution_comb, resolution_levels),
    change = case_when(
      is.na(rank_mid)  & !is.na(rank_comb) ~ "Gained by combined",
      !is.na(rank_mid) &  is.na(rank_comb) ~ "Lost in combined",
      rank_comb > rank_mid                  ~ "Improved",
      rank_comb == rank_mid                 ~ "Same",
      rank_comb < rank_mid                  ~ "Degraded",
      TRUE                                  ~ "Same"
    )
  )

# Summary of changes
otu_compare %>% dplyr::count(change, sort = TRUE)


# 4: Species-level comparison ----------------------------------------
# Collapse OTUs to unique species for each database

species_mid  <- mid  %>%
  filter(!is.na(Species), !Species %in% c("NoMatch", "Missing-Accession", "")) %>%
  distinct(Species, Genus, Family, Order) %>%
  arrange(Family, Species)

species_comb <- comb %>%
  filter(!is.na(Species), !Species %in% c("NoMatch", "Missing-Accession", "")) %>%
  distinct(Species, Genus, Family, Order) %>%
  arrange(Family, Species)

nrow(species_mid)   # unique species in MIDORI
nrow(species_comb)  # unique species in combined

# Species unique to each database
species_gained <- species_comb %>% filter(!Species %in% species_mid$Species)
species_lost   <- species_mid  %>% filter(!Species %in% species_comb$Species)

cat("Species gained by combined database:\n")
print(species_gained)

cat("\nSpecies in MIDORI only (not in combined):\n")
print(species_lost)


# 5: Key finding — Rasbora correction ----------------------------------------
# OTUs where combined changed the species assignment (not just improved rank)

species_corrections <- otu_compare %>%
  filter(
    !is.na(Species_mid), !is.na(Species_comb),
    !Species_mid %in% c("NoMatch", "Missing-Accession"),
    !Species_comb %in% c("NoMatch", "Missing-Accession"),
    Species_mid != Species_comb
  ) %>%
  select(OTU_ID, Species_mid, Species_comb, Similarity_mid, Similarity_comb)

cat("\nOTUs where species assignment changed between databases:\n")
print(species_corrections)


# 6: Summary table for paper ----------------------------------------

summary_table <- tibble(
  Metric = c(
    "Fish OTUs total",
    "Resolved to species",
    "Resolved to genus",
    "Resolved to family",
    "Unresolved",
    "Unique species detected",
    "Species-level rate (%)"
  ),
  MIDORI_only = c(
    nrow(mid),
    sum(mid$resolution == "Species"),
    sum(mid$resolution == "Genus"),
    sum(mid$resolution == "Family"),
    sum(mid$resolution == "Unresolved"),
    nrow(species_mid),
    round(100 * sum(mid$resolution == "Species") / nrow(mid), 1)
  ),
  Combined = c(
    nrow(comb),
    sum(comb$resolution == "Species"),
    sum(comb$resolution == "Genus"),
    sum(comb$resolution == "Family"),
    sum(comb$resolution == "Unresolved"),
    nrow(species_comb),
    round(100 * sum(comb$resolution == "Species") / nrow(comb), 1)
  )
)

print(summary_table)


# 7: Write outputs ----------------------------------------
write_csv(otu_compare,    file.path(output_dir, "db_comparison_otus.csv"))
write_csv(summary_table,  file.path(output_dir, "db_comparison_summary.csv"))
write_csv(species_gained, file.path(output_dir, "db_species_gained.csv"))

cat("Done. Outputs written to", output_dir, "\n")
cat("Next step: run 03_visualise_comparison.R\n")
