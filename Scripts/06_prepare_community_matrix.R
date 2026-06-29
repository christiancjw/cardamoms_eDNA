# 04_prepare_community_matrix.R
#
# Builds a presence/absence community matrix at event level.
# A species is PRESENT at an event if detected in EITHER replicate.
# Switch taxon_group to "amphibian" to rerun for frog dataset.
# ----------------------------------------------------------------

# Setup ----------------------------------------
library(dplyr)
library(tidyr)
library(readr)

# CONFIG ----------------------------------------
taxon_group <- "fish"

# 1: Read in data ----------------------------------------
otu  <- read_csv("Data/Clean/fish_combined.csv",       show_col_types = FALSE)
meta <- read_csv("Data/Raw/riparia2024_metadata.csv",  show_col_types = FALSE)
nrow(otu)
nrow(meta)

# Removing Kravanh and Opreil — single season sites, analysed separately --------
meta_excluded <- meta %>% filter(Site %in% c("Kravanh", "Opreil"))
meta          <- meta %>% filter(!Site %in% c("Kravanh", "Opreil"))
nrow(meta_excluded)
nrow(meta)

# 3: Get valid sample columns ----------------------------------------
sample_cols <- colnames(otu)[colnames(otu) %in% meta$Sample]
sample_cols <- setdiff(sample_cols, c("H06", "E07"))   # remove H06 (no metadata entry) and E07 
length(sample_cols)


# 4: Convert reads to presence/absence (1/0) ----------------------------------------
otu_pa <- otu %>%
  select(OTU_ID, Species, Genus, Family, Order, Class, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~ if_else(. > 0, 1L, 0L)))

# 5: Assign taxon labels ----------------------------------------
has_name <- function(x) !is.na(x) & !x %in% c("NoMatch", "Missing-Accession", "")

otu_pa <- otu_pa %>%
  mutate(taxon_label = case_when(
    has_name(Species) ~ Species,
    has_name(Genus)   ~ paste0(Genus,  " sp"),
    has_name(Family)  ~ paste0(Family, " sp"),
    has_name(Order)   ~ paste0(Order,  " sp"),
    has_name(Class)   ~ paste0(Class,  " sp"),
    TRUE              ~ paste0("Unresolved_", substr(OTU_ID, 1, 8))
  ))

# 6: Collapse OTUs sharing a taxon_label, per SAMPLE ----------------------------------------
# Multiple OTUs can map to the same taxon (e.g. two OTUs both = Glyptothorax laosensis).
# Present in a sample if ANY of its OTUs is present in that sample.
community_long <- otu_pa %>%
  pivot_longer(all_of(sample_cols), names_to = "Sample", values_to = "presence") %>%
  group_by(Sample, taxon_label) %>%
  summarise(presence = as.integer(any(presence == 1)), .groups = "drop")

# Spot check
community_long %>%
  filter(taxon_label == "Poropuntius normani") %>%
  arrange(Sample) %>%
  print(n = Inf)

# 7: Pivot to wide matrix (rows = samples, columns = taxa) ----------------------------------------
community_wide <- community_long %>%
  pivot_wider(names_from = taxon_label, values_from = presence, values_fill = 0)

dim(community_wide)   # Samples x (1 + Taxa)

# 8: Build sample-level metadata ----------------------------------------
sample_meta <- meta %>%
  filter(Sample %in% community_wide$Sample) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  arrange(match(Sample, community_wide$Sample)) %>%
  mutate(
    Strahler_f = factor(Strahler, levels = 1:5,
                        labels = c("1st","2nd","3rd","4th","5th"), ordered = TRUE),
    Season     = factor(Season, levels = c("Dry","Wet"))
  )

# Rows must align between matrix and metadata
identical(sample_meta$Sample, community_wide$Sample)   # must be TRUE

# Richness per sample
community_wide %>%
  mutate(richness = rowSums(across(-Sample))) %>%
  select(Sample, richness) %>%
  left_join(sample_meta %>% select(Sample, Event, Site, Strahler, Season), by = "Sample") %>%
  arrange(Strahler, Site, Sample) %>%
  print(n = Inf)

# 9: Write outputs ----------------------------------------
write_csv(community_wide, "Data/Clean/fish_community_pa_sample.csv")
write_csv(sample_meta,    "Data/Clean/fish_community_meta_sample.csv")

dim(community_wide)


# Event Level Outputs -------
community_long_event <- community_long %>%
  left_join(sample_meta %>% select(Sample, Event), by = "Sample") %>%
  group_by(Event, taxon_label) %>%
  summarise(presence = as.integer(any(presence == 1)), .groups = "drop")

community_wide_event <- community_long_event %>%
  pivot_wider(names_from = taxon_label, values_from = presence, values_fill = 0)

dim(community_wide_event)   # Events x (1 + Taxa)


# 11: Build event-level metadata ----------------------------------------
event_meta <- sample_meta %>%
  distinct(Event, .keep_all = TRUE) %>%
  right_join(community_wide_event %>% select(Event), by = "Event") %>%
  arrange(match(Event, community_wide_event$Event))

# Rows must align between matrix and metadata
identical(event_meta$Event, community_wide_event$Event)   # must be TRUE

# Richness per event
community_wide_event %>%
  mutate(richness = rowSums(across(-Event))) %>%
  select(Event, richness) %>%
  left_join(event_meta %>% select(Event, Site, Strahler, Season), by = "Event") %>%
  arrange(Strahler, Site) %>%
  print(n = Inf)


# 12: Write event-level outputs ----------------------------------------
write_csv(community_wide_event, "Data/Clean/fish_community_pa.csv")
write_csv(event_meta,           "Data/Clean/fish_community_meta.csv")
dim(community_wide_event)
