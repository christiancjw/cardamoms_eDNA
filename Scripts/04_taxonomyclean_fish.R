# =============================================================================
# 01_filter_fish.R
#
# Takes the two APSCALE taxonomy outputs (MIDORI-only and combined database)
# and the read table, filters to fish only (Actinopteri/Actinopterygii),
# removes contaminants and misidentifications, and produces clean OTU tables
# ready for database comparison and TTT2.
#
# NOTES:
#   - APSCALE uses "Actinopteri" (not "Actinopterygii") as the class name
#   - Both are checked to be safe
#   - Negative controls are removed if present
#   - The hash/unique ID columns are matched between taxonomy and read table
#   - Contaminant species are removed after fish filtering
# =============================================================================

# ----------------------------------------
# Setup
library(readxl)
library(dplyr)
library(readr)

taxonomy_midori   <- "Data/Raw/midori_blastn_0_24shaw16spaired_sequences_taxonomy.xlsx"
taxonomy_combined <- "Data/Raw/customdb_blastn_0_24shaw16spaired_sequences_taxonomy.xlsx"
read_table        <- "Data/Raw/0_24shaw16spaired_sequence_read_table_part_0.xlsx"
output_dir        <- "Data/Clean/"

# Negative control sample codes — remove these columns
# Set to c() if already excluded from read table
nc_samples <- c()


# ----------------------------------------
# Contaminant / misidentification species to remove
# These are removed from BOTH databases before any analysis.
#
# REMOVE — confirmed misidentifications (MIDORI database errors):
#   Trinectes inscriptus    — Atlantic flatfish, impossible in Mekong
#   Achiroides leucorhynchos — marine/estuarine sole, not in headwaters
#
# FLAG — introduced/non-native (kept in dataset, noted in paper):
#   Oreochromis niloticus      — African cichlid, aquaculture escapee
#   Hypophthalmichthys nobilis — Chinese bighead carp, stocked fisheries
#
# Add any additional confirmed misidentifications to remove_species below.

remove_species <- c(
  "Trinectes inscriptus",
  "Achiroides leucorhynchos"
)

introduced_species <- c(
  "Oreochromis niloticus",
  "Hypophthalmichthys nobilis"
)


# 2: Load files ----------------------------------------
tax_mid  <- read_excel(taxonomy_midori)
tax_comb <- read_excel(taxonomy_combined)
reads    <- read_excel(read_table)

# Check they loaded correctly
nrow(tax_mid)   # should be 229
nrow(tax_comb)  # should be 229
nrow(reads)     # should be 229


# 3: Standardise the ID column name ----------------------------------------
# Taxonomy uses "unique ID", read table uses "hash" — rename both to OTU_ID
tax_mid  <- tax_mid  %>% rename(OTU_ID = `unique ID`)
tax_comb <- tax_comb %>% rename(OTU_ID = `unique ID`)
reads    <- reads    %>% rename(OTU_ID = hash)

# Confirm the first few IDs match between files
head(tax_mid$OTU_ID)
head(reads$OTU_ID)


# 4: Filter to fish only ----------------------------------------
# APSCALE uses "Actinopteri" as the class name for ray-finned fish
# Check what class names are present first
tax_comb %>% dplyr::count(Class, sort = TRUE)

# Filter
fish_mid  <- tax_mid  %>% filter(Class %in% c("Actinopteri", "Actinopterygii"))
fish_comb <- tax_comb %>% filter(Class %in% c("Actinopteri", "Actinopterygii"))

# Check how many fish OTUs
nrow(fish_mid)
nrow(fish_comb)


# 5: Remove contaminant / misidentified species ----------------------------------------
# Check which contaminants are present before removing
fish_mid  %>% filter(Species %in% remove_species) %>% select(OTU_ID, Species, Similarity)
fish_comb %>% filter(Species %in% remove_species) %>% select(OTU_ID, Species, Similarity)

# Remove
fish_mid  <- fish_mid  %>% filter(!Species %in% remove_species)
fish_comb <- fish_comb %>% filter(!Species %in% remove_species)

# Confirm removal
nrow(fish_mid)
nrow(fish_comb)

# Flag introduced species (add column, keep in dataset)
fish_mid  <- fish_mid  %>% mutate(introduced = Species %in% introduced_species)
fish_comb <- fish_comb %>% mutate(introduced = Species %in% introduced_species)

# Check which introduced species are present
fish_comb %>% filter(introduced) %>% select(OTU_ID, Species, Similarity)


# 6: Remove negative controls if present ----------------------------------------
if (length(nc_samples) > 0) {
  nc_present <- nc_samples[nc_samples %in% colnames(reads)]
  if (length(nc_present) > 0) {
    reads <- reads %>% select(-all_of(nc_present))
    cat("Removed negative control columns:", paste(nc_present, collapse = ", "), "\n")
  }
}


# 7: Join read counts onto taxonomy ----------------------------------------
# Get sample columns from read table (everything except OTU_ID and sequence)
sample_cols <- reads %>% select(-OTU_ID, -sequence) %>% colnames()
sample_cols  # should be A01, A02 ... H08

# Join
fish_mid_reads  <- fish_mid  %>% left_join(reads %>% select(OTU_ID, all_of(sample_cols)), by = "OTU_ID")
fish_comb_reads <- fish_comb %>% left_join(reads %>% select(OTU_ID, all_of(sample_cols)), by = "OTU_ID")

# Add total reads column
fish_mid_reads  <- fish_mid_reads  %>% mutate(total_reads = rowSums(across(all_of(sample_cols)), na.rm = TRUE))
fish_comb_reads <- fish_comb_reads %>% mutate(total_reads = rowSums(across(all_of(sample_cols)), na.rm = TRUE))

# Quick check
fish_comb_reads %>% select(OTU_ID, Species, introduced, total_reads) %>% head(10)


# 8: Add resolution to the dataset -------------------------------------------
# Add resolution level
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

fish_mid_reads  <- add_resolution(fish_mid_reads)
fish_comb_reads <- add_resolution(fish_comb_reads)

# 9: Write outputs ----------------------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write_csv(fish_mid_reads,  file.path(output_dir, "fish_midori.csv"))
write_csv(fish_comb_reads, file.path(output_dir, "fish_combined.csv"))

cat("Done. fish_midori.csv and fish_combined.csv written to", output_dir, "\n")
