# ================================================================
# SPECIES LIST CONSTRUCTION, ANNOTATION, AND DETECTION PLOT
#
# This script runs in STAGES with MANUAL EDITING between some of them.
# Each stage writes a CSV; you hand-edit in Excel, then the next stage
# reads the edited file back in. Manual steps are flagged clearly.
#
#   STAGE 1  Build species table with site occurrences      -> fish_species_table.csv
#   [MANUAL] Add Corrsp (corrected names), review taxonomy
#   STAGE 2  Add FishBase references + historical status     -> fish_species_table_final.csv
#   [MANUAL] Add voucher rows, edit status categories, Corrsp
#   STAGE 3  Add IUCN Red List status                        -> fish_species_table_complete.csv
#   [MANUAL] Final review
#   STAGE 4  Detection dot-plot (proportion + voucher marks)
# ================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(forcats)


# ================================================================
# STAGE 1 — Build species table with per-site occurrences
# ================================================================

meta <- read.csv("Data/Raw/riparia2024_metadata.csv")
comb <- read.csv("Data/Clean/fish_combined.csv")

# Remove excluded samples
exclude_samples <- c("H06", "E07", "F09")
meta <- meta %>% filter(!Sample %in% exclude_samples)
comb <- comb %>% select(-any_of(exclude_samples))

# Site groups (from site codes)
ta_chey   <- meta$Sample[meta$Code %in% c("TC1", "TC2", "TC3")]
ta_say    <- meta$Sample[meta$Code %in% c("TS1", "TS2", "TS3")]
dam5      <- meta$Sample[meta$Code %in% c("TS4")]
tang_rang <- meta$Sample[meta$Code %in% c("SP1")]
pursat    <- meta$Sample[meta$Code %in% c("PS1")]

# Base taxonomy table + X/blank site occurrence columns
species_table <- comb %>%
  select(Class, Order, Family, Genus, Species, Similarity, introduced, resolution) %>%
  mutate(
    Ta_Chey   = ifelse(rowSums(comb[, ta_chey,   drop = FALSE]) > 0, "X", ""),
    Ta_Say    = ifelse(rowSums(comb[, ta_say,    drop = FALSE]) > 0, "X", ""),
    Dam5      = ifelse(rowSums(comb[, dam5,      drop = FALSE]) > 0, "X", ""),
    Tang_Rang = ifelse(rowSums(comb[, tang_rang, drop = FALSE]) > 0, "X", ""),
    Pursat    = ifelse(rowSums(comb[, pursat,    drop = FALSE]) > 0, "X", "")
  )

write.csv(species_table, "Data/Clean/fish_species_table.csv", row.names = FALSE)

# >>> MANUAL EDIT 1 <<<
# Open fish_species_table.csv, add a "Corrsp" column with corrected
# species names, review taxonomy. Save as fish_species_table_corrected.csv


# ================================================================
# STAGE 2 — FishBase references + historical status
# ================================================================
library(rfishbase)

species_table <- read.csv("Data/Clean/fish_species_table_corrected.csv")

# Collapse duplicate species rows (present if detected in any duplicate)
species_table <- species_table %>%
  filter(resolution == "Species", !is.na(Species)) %>%
  group_by(Species) %>%
  summarise(
    across(c(Class, Order, Family, Genus), first),
    Similarity = max(Similarity),
    introduced = first(introduced),
    across(c(Ta_Chey, Ta_Say, Dam5, Tang_Rang, Pursat),
           ~ ifelse(any(. == "X"), "X", "")),
    .groups = "drop"
  )

# Taxonomic authority from FishBase
fb_lookup <- species(species_table$Species, fields = c("Species", "Author")) %>%
  rename(Reference = Author)

species_table <- species_table %>% left_join(fb_lookup, by = "Species")

# Species FishBase couldn't match (sp./cf./undescribed/synonyms)
species_table %>% filter(is.na(Reference)) %>% pull(Species)

# Cross-reference against historical Cardamom records
historical <- read.csv("Data/historical_cardamom_fish_records.csv")
historical_species <- historical %>%
  distinct(Updated_Taxonomy) %>%
  rename(Species = Updated_Taxonomy) %>%
  mutate(historically_recorded = TRUE)

species_table <- species_table %>%
  left_join(historical_species, by = "Species") %>%
  mutate(
    historically_recorded   = ifelse(is.na(historically_recorded), FALSE, historically_recorded),
    detected_relevant_sites = (Ta_Chey == "X" | Ta_Say == "X" | Dam5 == "X"),
    status = case_when(
      detected_relevant_sites &  historically_recorded ~ "Range Extension / confirmation",
      detected_relevant_sites & !historically_recorded ~ "New Record",
      TRUE                                             ~ NA_character_
    )
  ) %>%
  select(Class, Order, Family, Genus, Species, Reference, Similarity, introduced,
         Ta_Chey, Ta_Say, Dam5, Tang_Rang, Pursat, status)

write.csv(species_table, "Data/Clean/fish_species_table_final.csv", row.names = FALSE)

# >>> MANUAL EDIT 2 <<<
# In fish_species_table_final.csv:
#  - Add Corrsp column (corrected display names)
#  - Append voucher rows: species with vouchers but NOT recovered by eDNA,
#    marked "V" in the site columns where the voucher was collected,
#    status = "Not recovered - voucher confirmed"
#  - Review/edit status categories
# Save as fish_species_table_EDITED.csv


# ================================================================
# STAGE 3 — IUCN Red List status
# ================================================================
library(rredlist)

df <- read.csv("Data/Clean/fish_species_table_EDITED.csv")

# Clean reference — remove parentheses
df <- df %>% mutate(Reference = str_replace_all(Reference, "[()]", ""))

# API key from https://apiv3.iucnredlist.org/api/v3/token
Sys.setenv(IUCN_REDLIST_KEY = "your_key_here")
iucn_key <- Sys.getenv("IUCN_REDLIST_KEY")

# Look up one species' Red List category (v4 API)
get_redlist <- function(sp) {
  if (is.na(sp) || str_detect(sp, "\\bsp\\b|\\bsp[0-9]|\\bcf\\b|/|\\bsp\\.") || sp == "") {
    return(NA_character_)
  }
  parts <- str_split(str_squish(sp), " ")[[1]]
  if (length(parts) < 2) return(NA_character_)
  res <- tryCatch(
    rl_species(genus = parts[1], species = parts[2], key = iucn_key),
    error = function(e) NULL
  )
  if (is.null(res) || length(res$assessments) == 0) return(NA_character_)
  res$assessments$red_list_category_code[1]
}

df <- df %>% mutate(redlist_status = vapply(Corrsp, get_redlist, character(1)))

df %>% select(Corrsp, Reference, redlist_status) %>% as.data.frame()

write.csv(df, "Data/Clean/fish_species_table_complete.csv", row.names = FALSE)

# >>> MANUAL EDIT 3 (optional) <<<
# Final review of fish_species_table_complete.csv


# ================================================================
# STAGE 4 — Detection dot-plot
#   Detected species: sized by proportion of samples detected
#   Voucher species (not eDNA-recovered): X marker where voucher collected
# ================================================================
# ================================================================
# STAGE 4 — Detection dot-plot
#   Detected species: sized by proportion of samples detected
#   Voucher species (not eDNA-recovered): X where voucher collected
# ================================================================

# --- 4a: Read the complete table for names + status ---
sp <- read_csv("Data/Clean/fish_species_table_complete.csv", show_col_types = FALSE) %>%
  mutate(Reference = str_replace_all(Reference, "[()]", ""))

site_cols  <- c("Ta_Chey", "Ta_Say", "Dam5", "Tang_Rang", "Pursat")
site_order <- c("Ta Chey", "Ta Say", "Dam 5", "Tang Rang", "Pursat")

# Label: "Corrsp (Reference)"
sp <- sp %>%
  mutate(sp_label = if_else(is.na(Reference) | Reference == "",
                            Corrsp, paste0(Corrsp, " (", Reference, ")")))


# --- 4b: Detected species — proportion of samples per site ---
comb_species <- comb %>% filter(resolution == "Species")

sample_lookup <- bind_rows(
  tibble(Sample = ta_chey,   Site = "Ta Chey"),
  tibble(Sample = ta_say,    Site = "Ta Say"),
  tibble(Sample = dam5,      Site = "Dam 5"),
  tibble(Sample = tang_rang, Site = "Tang Rang"),
  tibble(Sample = pursat,    Site = "Pursat")
)

all_samples <- sample_lookup$Sample

prop_table <- comb_species %>%
  group_by(Species) %>%
  summarise(across(all_of(all_samples), ~ as.integer(any(. > 0))), .groups = "drop") %>%
  pivot_longer(-Species, names_to = "Sample", values_to = "detected") %>%
  inner_join(sample_lookup, by = "Sample") %>%
  group_by(Species, Site) %>%
  summarise(proportion = mean(detected), .groups = "drop") %>%
  filter(proportion > 0) %>%
  left_join(sp %>% select(Species, Corrsp, sp_label, Order, Family), by = "Species") %>%
  filter(!is.na(sp_label))


# --- 4c: Voucher species (not eDNA-recovered) — X where voucher collected ---
voucher_long <- sp %>%
  filter(status == "Not recovered - voucher confirmed") %>%
  select(Order, Family, Corrsp, sp_label, all_of(site_cols)) %>%
  pivot_longer(all_of(site_cols), names_to = "Site_raw", values_to = "mark") %>%
  filter(mark == "V") %>%
  mutate(Site = recode(Site_raw,
                       Ta_Chey = "Ta Chey", Ta_Say = "Ta Say", Dam5 = "Dam 5",
                       Tang_Rang = "Tang Rang", Pursat = "Pursat")) %>%
  select(Order, Family, Corrsp, sp_label, Site)


# --- 4d: One canonical ordering for ALL species (detected + voucher) ---
# Combine BEFORE arranging so each family forms a single contiguous block
all_species <- bind_rows(
  prop_table   %>% distinct(Order, Family, Corrsp, sp_label),
  voucher_long %>% distinct(Order, Family, Corrsp, sp_label)
) %>%
  distinct(sp_label, .keep_all = TRUE) %>%
  arrange(Order, Family, Corrsp)

# Axis order: reverse so first species sits at the top of the plot
species_levels <- rev(all_species$sp_label)


# --- 4e: Apply the shared factor axis ---
prop_table <- prop_table %>%
  mutate(sp_label = factor(sp_label, levels = species_levels),
         Site     = factor(Site, levels = site_order),
         marker   = "Detected")

voucher_long <- voucher_long %>%
  mutate(sp_label   = factor(sp_label, levels = species_levels),
         Site       = factor(Site, levels = site_order),
         proportion = NA_real_,
         marker     = "Voucher collected, not eDNA-detected")

plot_data <- bind_rows(prop_table, voucher_long)


# --- 4f: Family divider lines and labels (numeric positions, one source) ---
family_lookup <- all_species %>%
  mutate(y_pos = match(sp_label, species_levels)) %>%
  arrange(y_pos)

# Divider lines between families (top edge of each block)
family_breaks <- family_lookup %>%
  mutate(change = Family != lag(Family)) %>%
  filter(change) %>%
  pull(y_pos) - 0.5

# One label per family at the top species position of its block
family_labels <- family_lookup %>%
  group_by(Family) %>%
  summarise(y_pos = max(y_pos), .groups = "drop")

# Sanity check — families should be contiguous (returns 0 rows if all good)
family_lookup %>%
  group_by(Family) %>%
  summarise(min = min(y_pos), max = max(y_pos), n = n(), .groups = "drop") %>%
  mutate(contiguous = (max - min + 1) == n) %>%
  filter(!contiguous)


# --- 4g: Plot ---
site_colours <- c(
  "Ta Chey"   = "#4FA873",
  "Ta Say"    = "#5E9BD4",
  "Dam 5"     = "#98319B",
  "Tang Rang" = "#F2A35C",
  "Pursat"    = "#C0392B"
)

fig_occurrence <- ggplot(plot_data, aes(x = Site, y = sp_label)) +
  # Family divider lines — numeric intercepts on the discrete y-axis
  geom_hline(yintercept = family_breaks, linetype = "dotted",
             colour = "grey60", linewidth = 0.3) +
  # Detected species bubbles
  geom_point(data = plot_data %>% filter(marker == "Detected"),
             aes(size = proportion, fill = Site),
             shape = 21, colour = "black", stroke = 0.3) +
  # Voucher X marks
  geom_point(data = plot_data %>% filter(marker != "Detected"),
             aes(shape = marker),
             size = 1.8, colour = "grey20", stroke = 0.5) +
  # Family labels — numeric y_pos matches the numeric hline breaks
  geom_text(data = family_labels,
            aes(y = y_pos, label = Family),
            x = 0.15, inherit.aes = FALSE,
            hjust = 0, vjust = 0.5,
            size = 2.2, fontface = "bold", colour = "grey30") +
  scale_size(range = c(0.2, 2), name = "Proportion of samples detected") +
  scale_fill_manual(values = site_colours, guide = "none") +
  scale_shape_manual(values = c("Voucher collected, not eDNA-detected" = 4),
                     name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.9, 0.4))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Species") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y          = element_text(size = 7, face = "italic"),
    axis.text.x          = element_text(angle = 45, hjust = 1),
    axis.title.y         = element_text(size = 9),
    axis.title.x         = element_text(size = 9),
    legend.title         = element_text(size = 8),
    legend.text          = element_text(size = 8),
    legend.key.size      = unit(0.4, "cm"),
    legend.position      = "bottom",
    legend.box           = "horizontal",
    legend.justification = "left",
    legend.box.just      = "left",
    legend.box.margin    = margin(l = -85)
  )

fig_occurrence

ggsave(fig_occurrence, file="a4_output.pdf", width = 190, height = 297, units = "mm")

