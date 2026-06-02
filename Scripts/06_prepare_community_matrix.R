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


# 2: Set aside Kravanh and Opreil — single season sites, analysed separately --------
meta_excluded <- meta %>% filter(Site %in% c("Kravanh", "Opreil"))
meta          <- meta %>% filter(!Site %in% c("Kravanh", "Opreil"))

nrow(meta_excluded)  # samples set aside
nrow(meta)           # remaining samples


# 3: Get valid sample columns ----------------------------------------
# Sample columns match metadata samples
sample_cols <- colnames(otu)[colnames(otu) %in% meta$Sample]
sample_cols <- setdiff(sample_cols, "H06") # remove H06 (no entry in metadata anyway)


length(sample_cols)  # confirm count


# 4: Convert reads to presence/absence (1/0) ----------------------------------------
otu_pa <- otu %>%
  select(OTU_ID, Species, Genus, Family, Order, Class, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), ~if_else(. > 0, 1L, 0L)))


# 5: Assign taxon labels ----------------------------------------
# Returns TRUE if a taxonomy column has a real value
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

# 6: Merge replicates to event level ----------------------------------------
# For each event (e.g. TC2), combine replicates A01 + A02
# A taxon is present at the event if it appears in EITHER replicate

sample_event <- meta %>% select(Sample, Event, Site, Strahler, QBR,
                                Season, Deployment, QBR_Class, Strahler_Class,
                                Month, Year, pH, `us/cm`, ppm, Long, Lat, Code)

# Convert to long form
community_long <- otu_pa %>%
  pivot_longer(all_of(sample_cols), names_to = "Sample", values_to = "presence") %>%
  left_join(sample_event %>% select(Sample, Event), by = "Sample") %>%
  group_by(Event, taxon_label) %>%
  summarise(presence = as.integer(any(presence == 1)), .groups = "drop")

# Spot check
community_long %>% 
  filter(taxon_label == "Poropuntius normani") %>% 
  print(n = 35)


# 7: Pivot to wide matrix (rows = events, columns = taxa) ----------------------------------------
community_wide <- community_long %>%
  pivot_wider(names_from = taxon_label, values_from = presence, values_fill = 0)

# Dimensions: Events x Taxa
dim(community_wide)


# 8: Build event-level metadata ----------------------------------------
event_meta <- sample_event %>%
  distinct(Event, .keep_all = TRUE) %>%
  right_join(community_wide %>% select(Event), by = "Event") %>%
  arrange(match(Event, community_wide$Event)) %>%
  mutate(
    Strahler_f = factor(Strahler, levels = 1:5,
                        labels = c("1st","2nd","3rd","4th","5th"), ordered = TRUE),
    Season     = factor(Season, levels = c("Dry","Wet"))
  )

# Rows must align between matrix and metadata
identical(event_meta$Event, community_wide$Event)  # must be TRUE

# Richness per event
community_wide %>%
  mutate(richness = rowSums(across(-Event))) %>%
  select(Event, richness) %>%
  left_join(event_meta %>% select(Event, Site, Strahler, Season), by = "Event") %>%
  arrange(Strahler, Site) %>%
  print(n = 35)


# 9: Write outputs ----------------------------------------
write_csv(community_wide, "Data/Clean/fish_community_pa.csv")
write_csv(event_meta,     "Data/Clean/fish_community_meta.csv")

# Confirms: n events x n taxa
dim(community_wide)
print(community_wide)
