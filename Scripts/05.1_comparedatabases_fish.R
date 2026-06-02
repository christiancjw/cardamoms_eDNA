# ----------------------------------------------------------------
# 02_compare_databases.R
# ----------------------------------------------------------------

# Setup ----------------------------------------
library(dplyr)
library(readr)
library(ggplot2)


# 1: Load ----------------------------------------
mid  <- read_csv("Data/Clean/fish_midori.csv",  show_col_types = FALSE)
comb <- read_csv("Data/Clean/fish_combined.csv", show_col_types = FALSE)

nrow(mid)
nrow(comb)


# 2: Add resolution ----------------------------------------
add_resolution <- function(df) {
  df %>% mutate(resolution = case_when(
    !is.na(Species) & !Species %in% c("NoMatch","Missing-Accession","") ~ "Species",
    !is.na(Genus)   & !Genus   %in% c("NoMatch","Missing-Accession","") ~ "Genus",
    !is.na(Family)  & !Family  %in% c("NoMatch","Missing-Accession","") ~ "Family",
    !is.na(Order)   & !Order   %in% c("NoMatch","Missing-Accession","") ~ "Order",
    TRUE ~ "Unresolved"
  ))
}

mid  <- add_resolution(mid)
comb <- add_resolution(comb)

mid  %>% dplyr::count(resolution, sort = TRUE)
comb %>% dplyr::count(resolution, sort = TRUE)


# 3: Fix species names in combined ----------------------------------------

comb <- comb %>%
  mutate(Species = case_when(
    Species == "Nemacheilus masyai/sp" ~ "Nemacheilus masyai",
    Species == "Clarias sp"            ~ "Clarias cf nieuhofii",
    TRUE                               ~ Species
  ))

# Check species names of these
comb %>% filter(grepl("Nemacheilus|Clarias|Schistura", Species)) %>% dplyr::count(Species)


# 4: Join and classify each OTU ----------------------------------------
otu_compare <- mid %>%
  select(OTU_ID, Species_mid = Species, Similarity_mid = Similarity, resolution_mid = resolution) %>%
  full_join(
    comb %>% select(OTU_ID, Species_comb = Species, Similarity_comb = Similarity, resolution_comb = resolution),
    by = "OTU_ID"
  )

# Type 1 — Species reassigned (wrong taxon corrected)
type1_species <- c(
  "Channa limbata",
  "Rasbora elegans/sumatrana",
  "Glyptothorax ventrolineatus",
  "Clarias macrocephalus",
  "Osteochilus hasseltii",
  "Systomus orphoides",
  "Schistura callichromus/fasciolata"
)

# Type 2 — Same species, better similarity via local voucher
type2_species <- c(
  "Pseudomystus siamensis",
  "Mystacoleucus marginatus",
  "Homaloptera confuzona",
  "Parambassis siamensis"
)

otu_compare <- otu_compare %>%
  mutate(improvement_type = case_when(
    Species_mid %in% type1_species ~ "Type 1 — Reassigned",
    Species_mid %in% type2_species & Similarity_comb > Similarity_mid ~ "Type 2 — Confidence",
    TRUE ~ NA_character_
  ))

# Check
otu_compare %>% filter(!is.na(improvement_type)) %>% dplyr::count(improvement_type)

otu_compare %>%
  filter(improvement_type == "Type 1 — Reassigned") %>%
  select(Species_mid, Species_comb, Similarity_mid, Similarity_comb)

otu_compare %>%
  filter(improvement_type == "Type 2 — Confidence") %>%
  select(Species_mid, Species_comb, Similarity_mid, Similarity_comb)


# 5: Summary table ----------------------------------------
species_mid  <- mid  %>% filter(resolution == "Species") %>% distinct(Species)
species_comb <- comb %>% filter(resolution == "Species") %>% distinct(Species)

tibble(
  Metric = c("Fish OTUs","Species-level","Genus-level","Family-level",
             "Unresolved","Unique species","Species rate (%)"),
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


# 6: Plot ----------------------------------------
# Colours
  res_colours <- c(
    "Species"           = "#1A6BAD",   # dark blue
    "Species (Type 1)"  = "#5AADD4",   # mid blue
    "Species (Type 2)"  = "#A8D4EC",   # light blue
    "Genus"             = "#66C2A5",   # teal  (Set2)
    "Family"            = "#FC8D62",   # orange (Set2)
    "Order"             = "#8DA0CB",   # lavender (Set2)
    "Unresolved"        = "#D9D9D9"   
)

res_levels <- c("Species","Species (Type 1)","Species (Type 2)",
                "Genus","Family","Order","Unresolved")

n_type1   <- sum(otu_compare$improvement_type == "Type 1 — Reassigned",  na.rm = TRUE)
n_type2   <- sum(otu_compare$improvement_type == "Type 2 — Confidence",   na.rm = TRUE)
n_species <- sum(mid$resolution == "Species") - n_type1 - n_type2

# MIDORI bar — species split into 3 segments
midori_bar <- bind_rows(
  data.frame(resolution = "Species",          n = n_species),
  data.frame(resolution = "Species (Type 1)", n = n_type1),
  data.frame(resolution = "Species (Type 2)", n = n_type2),
  mid %>% dplyr::count(resolution) %>% filter(resolution != "Species")
) %>% mutate(Database = "MIDORI only")

# Combined bar — species as single block
comb_bar <- comb %>%
  dplyr::count(resolution) %>%
  mutate(Database = "MIDORI + Local")

plot_data <- bind_rows(midori_bar, comb_bar) %>%
  mutate(
    resolution = factor(resolution, levels = res_levels),
    Database   = factor(Database, levels = c("MIDORI only", "MIDORI + Local"))
  )

fig_resolution <- ggplot(plot_data, aes(x = Database, y = n, fill = resolution)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = ifelse(n > 1, n, "")),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = res_colours, name = NULL,
                    labels = c(
                      "Species"          = "Species",
                      "Species (Type 1)" = "Sp. Incorrect Match",
                      "Species (Type 2)" = "Sp. Confidence Increase",
                      "Genus"            = "Genus",
                      "Family"           = "Family",
                      "Order"            = "Order",
                      "Unresolved"       = "Unresolved"
                    )) +
  labs(x = NULL, y = "Number of OTUs") +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "right",
    legend.direction = "vertical"
  )

fig_resolution

