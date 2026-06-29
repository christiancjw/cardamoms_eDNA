# ----------------------------------------------------------------
# 02_compare_databases.R
# ----------------------------------------------------------------

# Setup ----------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

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
  Metric = c("Fish OTUs","Species-level","Genus-level","Family-level", "Order-level",
             "Unresolved","Unique species","Species rate (%)"),
  MIDORI_only = c(
    nrow(mid),
    sum(mid$resolution == "Species"),
    sum(mid$resolution == "Genus"),
    sum(mid$resolution == "Family"),
    sum(comb$resolution == "Order"),
    sum(mid$resolution == "Unresolved"),
    nrow(species_mid),
    round(100 * sum(mid$resolution == "Species") / nrow(mid), 1)
  ),
  Combined = c(
    nrow(comb),
    sum(comb$resolution == "Species"),
    sum(comb$resolution == "Genus"),
    sum(comb$resolution == "Family"),
    sum(comb$resolution == "Order"),
    sum(comb$resolution == "Unresolved"),
    nrow(species_comb),
    round(100 * sum(comb$resolution == "Species") / nrow(comb), 1)
  )
)


# 6: Plot OTUs ----------------------------------------
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

res_colours <- c(
  "Species"           = "#3F007D",   # dark purple   ┐
  "Species (Type 1)"  = "#6A51A3",   # mid purple    ├ one hue, 3 shades
  "Species (Type 2)"  = "#9E9AC8",   # light purple  ┘
  "Genus"             = "#E7298A",   # magenta-pink  — distinct hue
  "Family"            = "#E6AB02",   # mustard/gold  — distinct hue
  "Order"             = "#1B9E77",   # teal-green    — distinct hue
  "Unresolved"        = "#D9D9D9"    # neutral grey  — shared "no call"
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
  labs(x = NULL, y = "Number of OTUs", subtitle = "BLASTn Detections by Database") +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "right",
    legend.direction = "vertical", plot.subtitle    = element_text(face = "bold") 
  )

fig_resolution

## Species collapsed plot:
# 7: Collapse to species ----------------------------------------
# Species-rank OTUs -> one row per distinct species name.
# Genus / Family / Order / Unresolved -> kept as-is (one row per OTU).
# NOTE: relies on clean Species strings. Fix any concatenation
# artifacts (e.g. "Mystus atrifasciatus/multiradiatusMystus
# atrifasciatus/multiradiatus") before this step or distinct()
# will treat the doubled string as its own species.

collapse_to_species <- function(df, db_label) {
  species_n <- df %>%
    filter(resolution == "Species") %>%
    distinct(Species) %>%
    nrow()
  
  other_n <- df %>%
    filter(resolution != "Species") %>%
    dplyr::count(resolution)
  
  bind_rows(
    tibble(resolution = "Species", n = species_n),
    other_n
  ) %>%
    mutate(Database = db_label)
}

species_collapsed <- bind_rows(
  collapse_to_species(mid,  "MIDORI only"),
  collapse_to_species(comb, "MIDORI + Local")
) %>%
  mutate(
    resolution = factor(resolution,
                        levels = c("Species","Genus","Family","Order","Unresolved")),
    Database   = factor(Database,
                        levels = c("MIDORI only","MIDORI + Local"))
  )

species_collapsed

# Quick check: OTUs at species rank vs distinct species (the collapse)
mid  %>% summarise(otu_species = sum(resolution == "Species"),
                   distinct_species = n_distinct(Species[resolution == "Species"]))
comb %>% summarise(otu_species = sum(resolution == "Species"),
                   distinct_species = n_distinct(Species[resolution == "Species"]))


# 7: Collapse to species (with improvement type preserved) ------------------
# Type is resolved per collapsed species: Type 1 if ANY of its OTUs was
# reassigned, else Type 2 if ANY was a confidence gain, else plain Species.
# Requires clean Species strings first (fix concatenation artifacts).

# Per-OTU improvement type, keyed to the MIDORI species name
otu_type <- otu_compare %>%
  select(Species_mid, improvement_type)

# MIDORI species, each tagged with its strongest improvement category
mid_species_typed <- mid %>%
  filter(resolution == "Species") %>%
  left_join(otu_type, by = c("Species" = "Species_mid")) %>%
  group_by(Species) %>%
  summarise(
    seg = case_when(
      any(improvement_type == "Type 1 — Reassigned", na.rm = TRUE) ~ "Species (Type 1)",
      any(improvement_type == "Type 2 — Confidence", na.rm = TRUE) ~ "Species (Type 2)",
      TRUE ~ "Species"
    ),
    .groups = "drop"
  )

midori_species_bar <- bind_rows(
  mid_species_typed %>% dplyr::count(seg, name = "n") %>% rename(resolution = seg),
  mid %>% filter(resolution != "Species") %>% dplyr::count(resolution)
) %>% mutate(Database = "MIDORI only")

# Combined bar — species collapsed, single block (Types are a MIDORI-side story)
comb_species_bar <- bind_rows(
  tibble(resolution = "Species",
         n = n_distinct(comb$Species[comb$resolution == "Species"])),
  comb %>% filter(resolution != "Species") %>% dplyr::count(resolution)
) %>% mutate(Database = "MIDORI + Local")

species_collapsed <- bind_rows(midori_species_bar, comb_species_bar) %>%
  mutate(
    resolution = factor(resolution, levels = res_levels),
    Database   = factor(Database, levels = c("MIDORI only","MIDORI + Local"))
  )

# Sanity: typed species segments must sum to distinct MIDORI species
stopifnot(
  sum(mid_species_typed$seg %in% c("Species","Species (Type 1)","Species (Type 2)")) ==
    n_distinct(mid$Species[mid$resolution == "Species"])
)

species_collapsed


# 8: Plot species ----------------------------------------
fig_species <- ggplot(species_collapsed, aes(x = Database, y = n, fill = resolution)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = ifelse(n > 3, n, "")),
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
  labs(x = NULL, y = "Number of taxa", subtitle = "Blast Detections by Database") +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "right",
    legend.direction = "vertical", plot.subtitle    = element_text(face = "bold") 
  )

fig_species



# 9: Number of species ----------------------------------------
# Distinct named species = species-rank OTUs collapsed to unique names.
# (Genus/family/order OTUs are not species and are excluded here.)

n_species_mid  <- mid  %>% filter(resolution == "Species") %>% distinct(Species) %>% nrow()
n_species_comb <- comb %>% filter(resolution == "Species") %>% distinct(Species) %>% nrow()

cat("Distinct named species (MIDORI only)   :", n_species_mid,  "\n")
cat("Distinct named species (MIDORI + Local):", n_species_comb, "\n")

# Total taxa / MOTUs = distinct named species + every non-species OTU kept as its own unit
# (this is the unit count that should feed richness if you count MOTUs, not just names)
n_taxa_comb <- n_species_comb + sum(comb$resolution != "Species")
cat("Total taxa / MOTUs (MIDORI + Local)    :", n_taxa_comb, "\n")

# Breakdown by rank, for the methods/results sentence
comb %>% dplyr::count(resolution, sort = TRUE)

# 10: Whole-basin composition (stacked, combined database) --------------------
# Publication figure. Unit toggle: "OTU" (rows) or "species" (collapsed).
#   - "OTU"     : each OTU counted once at its rank (detection representation)
#   - "species" : species-rank OTUs collapsed to distinct names; non-species
#                 OTUs retained as individual MOTUs at their rank
# Unresolved is always forced to the end of the stack and greyed.

library(dplyr)
library(ggplot2)
library(forcats)

unit <- "species"   # <- set "OTU" or "species"
blank_vals <- c("NoMatch","Missing-Accession","")

clean_rank <- function(x) ifelse(is.na(x) | x %in% blank_vals, "Unresolved", x)

# Build the counting table for a given rank
basin_comp <- function(df, rank_col, unit) {
  df <- df %>% mutate(.rank = clean_rank(.data[[rank_col]]))
  
  if (unit == "species") {
    # species-rank rows -> one per distinct species; others -> one per OTU
    sp  <- df %>% filter(resolution == "Species") %>%
      distinct(Species, .rank)
    oth <- df %>% filter(resolution != "Species")
    df  <- bind_rows(sp, oth)
  }
  
  df %>%
    dplyr::count(.rank, name = "n") %>%
    arrange(.rank == "Unresolved", desc(n)) %>%      # Unresolved last
    mutate(.rank = fct_inorder(.rank))
}

# ---- Order palette (resolved orders get colour, Unresolved greyed) ----------
# Colourblind-safe, distinct hues. Extend if you have >8 orders.
order_levels <- basin_comp(comb, "Order", unit)$.rank
order_resolved <- setdiff(as.character(order_levels), "Unresolved")

order_pal_base <- c(
  "#1F6FB2",  # blue
  "#2CA25F",  # green
  "#D95F0E",  # orange
  "#6A51A3",  # purple
  "#C51B7D",  # magenta
  "#41B6C4",  # teal
  "#B2912F",  # gold
  "#762A83",
  "#E7298A",  
  "#66A61E"   
)
order_pal <- setNames(order_pal_base[seq_along(order_resolved)], order_resolved)
order_pal["Unresolved"] <- "#D9D9D9"

# ---- Plot helper ------------------------------------------------------------
y_lab <- if (unit == "species") "Number of taxa" else "Number of OTUs"

comp_stacked <- function(df, fill_vals = NULL, subtitle = NULL) {
  p <- ggplot(df, aes(x = "Pursat Basin", y = n, fill = .rank)) +
    geom_col(width = 0.45, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(n >= 4, n, "")),
              position = position_stack(vjust = 0.5),
              size = 3, colour = "white", fontface = "bold") +
    scale_y_continuous(limits = c(0, 120), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = y_lab, fill = NULL, subtitle = subtitle) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x      = element_text(size = 11),
      axis.title.y     = element_text(size = 11),
      legend.text      = element_text(size = 10, face = "italic"),
      legend.key.size  = unit(0.9, "lines"),
      plot.subtitle    = element_text(face = "bold", size = 12),
      legend.position  = "right"
    ) +
    guides(fill = guide_legend(ncol = 1))
  if (!is.null(fill_vals)) p <- p + scale_fill_manual(values = fill_vals)
  p
}

# ---- Order: defined palette -------------------------------------------------
fig_order <- comp_stacked(basin_comp(comb, "Order", unit),
                          fill_vals = order_pal, subtitle = "Pursat Basin Fish Composition by Order")
fig_order

# ---- Family: many levels, use viridis (scales to any count) -----------------
fam_df <- basin_comp(comb, "Family", unit)
fam_resolved <- setdiff(as.character(fam_df$.rank), "Unresolved")
fam_pal <- setNames(
  viridisLite::viridis(length(fam_resolved), option = "D", end = 0.92),
  fam_resolved
)
fam_pal["Unresolved"] <- "#D9D9D9"

fig_family <- comp_stacked(fam_df, fill_vals = fam_pal, subtitle = "Family")
fig_family



# ---- Save (publication) -----------------------------------------------------
# ggsave("Figs/composition_order.png",  fig_order,  width = 4, height = 5, dpi = 600)
# ggsave("Figs/composition_family.png", fig_family, width = 4.5, height = 6, dpi = 600)


patchworkGrob(
  fig_order
  
)