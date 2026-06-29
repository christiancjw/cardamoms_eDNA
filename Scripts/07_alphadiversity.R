# ----------------------------------------------------------------
# 05_alpha_diversity.R
#
# Alpha diversity — species richness across the hydrological gradient.
# Shannon is not reported (mathematically identical to ln(S) on P/A data).
# ----------------------------------------------------------------

# Setup ----------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(vegan)
library(cowplot)
library(tidyr)
library(iNEXT)
library(sf)
library(ggrepel)
library(ggspatial)
library(lme4)
library(lmerTest)   # adds p-values via Satterthwaite

install.packages("ggstar")
library(ggstar)

# Colour palettes ----------------------------------------
strahler_f_colors <- c(
  "1st" = "#C6DCEF",
  "2nd" = "#84B9D9",
  "3rd" = "#4393C3",
  "4th" = "#1E6FA8",
  "5th" = "#084A87"
)

season_colors <- c("Dry" = "#D4A94A", "Wet" = "#4A90D4")

qbr_colors <- c(
  "Natural" = "#006BA6", "Good" = "#22A122",
  "Fair"    = "#DBCB43", "Poor" = "#FF7134", "Bad" = "#AF3245"
)

site_order <- c(
  "Knaong Bat Sa", "Stung Ta Chey", "Stung Ta Say", "Arai Village",
  "Kronomh Village", "Stung Oda", "Dam 5", "Tang Rang Bridge",
  "Pursat River"
)


# 1: Read in data ----------------------------------------
pa   <- read_csv("Data/Clean/fish_community_pa.csv",   show_col_types = FALSE)
meta <- read_csv("Data/Clean/fish_community_meta.csv", show_col_types = FALSE)

identical(pa$Event, meta$Event)  # must be TRUE


# 2: Calculate richness ----------------------------------------
taxa_cols <- colnames(pa)[colnames(pa) != "Event"]
otu_mat   <- pa %>% select(all_of(taxa_cols)) %>% as.matrix()

alpha_div <- meta %>%
  mutate(
    richness   = specnumber(otu_mat),
    Site       = factor(Site, levels = site_order),
    Strahler_f = factor(Strahler, levels = 1:5,
                        labels = c("1st","2nd","3rd","4th","5th"), ordered = TRUE),
    Season     = factor(Season, levels = c("Dry","Wet")),
    QBR_Class  = factor(QBR_Class, levels = c("Natural","Good","Fair","Poor","Bad"))
  )

alpha_div %>% select(Event, Site, Strahler_f, Season, richness)


# 3: Linear models ----------------------------------------



# Richness ~ Strahler (continuous — tests linear gradient)
lm_strahler <- lm(richness ~ Strahler, data = alpha_div)
summary(lm_strahler)

# Richness ~ Season
lm_season <- lm(richness ~ Season, data = alpha_div)
summary(lm_season)

# Richness ~ Strahler + Season (additive)
lm_both <- lm(richness ~ Strahler + Season, data = alpha_div)
summary(lm_both)

# Richness ~ QBR (for confounding comparison)
lm_qbr <- lm(richness ~ QBR, data = alpha_div)
summary(lm_qbr)


# 4: Figure — Richness by Strahler order ----------------------------------------
fig_richness_strahler <- ggplot(alpha_div,
                                aes(x = Strahler_f, y = richness, fill = Strahler_f)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(colour = Season), width = 0.15, size = 2.5, alpha = 0.9) +
  scale_fill_manual(values = strahler_f_colors, guide = "none") +
  scale_colour_manual(values = season_colors, name = NULL) +
  labs(x = "Strahler order", y = "Species richness (S)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

fig_richness_strahler


# 5: Figure — Richness vs QBR (confounding) ----------------------------------------
fig_richness_qbr <- ggplot(alpha_div,
                           aes(x = QBR, y = richness, colour = Strahler_f, shape = Season)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              linetype = "dashed", linewidth = 0.7) +
  scale_colour_manual(values = strahler_f_colors, name = "Strahler order") +
  scale_shape_manual(values = c("Dry" = 16, "Wet" = 17), name = NULL) +
  labs(x = "QBR score", y = "Species richness (S)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

fig_richness_qbr


# 6: Species accumulation curves ----------------------------------------

# Build specaccum for a given Strahler range
make_sacc <- function(min_s, max_s) {
  events <- meta %>% filter(Strahler >= min_s, Strahler <= max_s) %>% pull(Event)
  mat    <- pa %>% filter(Event %in% events) %>% select(all_of(taxa_cols)) %>% as.matrix()
  specaccum(mat, method = "random", permutations = 999)
}

sacc_full <- make_sacc(1, 5)
sacc_head <- make_sacc(1, 3)
sacc_low  <- make_sacc(4, 5)


# Saturation statistics for each curve ----------------------------------------
# Fits a Michaelis-Menten asymptote to estimate total richness and 95% saturation point

sacc_stats <- function(sacc, label) {
  
  # Fit asymptote using nonlinear least squares (Michaelis-Menten)
  df  <- data.frame(x = sacc$sites, y = sacc$richness)
  fit <- tryCatch(
    nls(y ~ Smax * x / (k + x), data = df, start = list(Smax = max(df$y) * 1.2, k = 5)),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    cat(sprintf("\n%s: could not fit asymptote\n", label))
    return(invisible(NULL))
  }
  
  Smax      <- coef(fit)["Smax"]                 # estimated total richness
  k         <- coef(fit)["k"]                    # half-saturation constant
  events_95 <- ceiling(k * 0.95 / (1 - 0.95))   # events needed for 95% of Smax
  
  cat(sprintf(
    "\n%s:\n  Estimated total richness (asymptote): %.1f\n  Events to reach 95%%: %d\n  Observed at n=%d: %d species\n",
    label, Smax, events_95, max(sacc$sites), round(max(sacc$richness))
  ))
  
  list(Smax = Smax, k = k, events_95 = events_95, fit = fit, label = label)
}

stats_full <- sacc_stats(sacc_full, "Full basin (1–5)")
stats_head <- sacc_stats(sacc_head, "Cardamom headwaters (1–3)")
stats_low  <- sacc_stats(sacc_low,  "Lowland (4–5)")


# SOme statistics for ya - average sp per sample / stdv
per_sample <- function(sacc, label) {
  cat(sprintf("%s: %.1f species per event (SD %.1f)\n",
              label, sacc$richness[1], sacc$sd[1]))
}

per_sample(sacc_full, "Full basin (1–5)")
per_sample(sacc_head, "Cardamom headwaters (1–3)")
per_sample(sacc_low,  "Lowland (4–5)")
# Build tidy data frame for plotting ----------------------------------------

# iNEXT takes a list of incidence frequency data
# For P/A event data: convert each subset to a list of species incidence counts
# (number of events each species was detected in)

make_inext <- function(min_s, max_s) {
  events <- meta %>% filter(Strahler >= min_s, Strahler <= max_s) %>% pull(Event)
  mat    <- pa %>% filter(Event %in% events) %>% select(all_of(taxa_cols))
  # First element = number of sampling events (T), rest = species incidence counts
  inc <- c(nrow(mat), colSums(mat))
  inc[inc > 0]  # drop species never detected
}

inc_full <- make_inext(1, 5)
inc_head <- make_inext(1, 3)
inc_low  <- make_inext(4, 5)

# Run iNEXT — datatype "incidence_freq", extrapolate to 2x observed sample size
inext_out <- iNEXT(
  list(
    "Full basin (Strahler 1–5)"          = inc_full,
    "Cardamom headwaters (Strahler 1–3)" = inc_head,
    "Lowland modified (Strahler 4–5)"    = inc_low
  ),
  q          = 0,           # species richness (not Shannon or Simpson)
  datatype   = "incidence_freq",
  endpoint   = max(nrow(pa) * 2, 60)  # extrapolate to 2x full basin sample size
)

# Print summary statistics — asymptote and coverage estimates
print(inext_out$AsyEst)

# Extract tidy data frame for ggplot
sacc_df <- inext_out$iNextEst$size_based %>%
  mutate(
    Assemblage = factor(Assemblage, levels = c(
      "Full basin (Strahler 1–5)",
      "Cardamom headwaters (Strahler 1–3)",
      "Lowland modified (Strahler 4–5)"
    )),
    # Mark interpolated vs extrapolated
    line_type = if_else(Method == "Rarefaction", "Interpolated", "Extrapolated")
  )

sacc_colors <- c(
  "Full Basin (1–5)"              = "darkgreen",
  "Cardamom Headwaters (1–3)"     = "#0073C2FF",
  "Lowland Modified Rivers (4–5)" = "#CD534CFF"
)

sacc_df <- inext_out$iNextEst$size_based %>%
  mutate(
    Assemblage = case_when(
      Assemblage == "Full basin (Strahler 1–5)"          ~ "Full Basin (1–5)",
      Assemblage == "Cardamom headwaters (Strahler 1–3)" ~ "Cardamom Headwaters (1–3)",
      Assemblage == "Lowland modified (Strahler 4–5)"    ~ "Lowland Modified Rivers (4–5)"
    ),
    Assemblage = factor(Assemblage, levels = c(
      "Full Basin (1–5)", "Cardamom Headwaters (1–3)", "Lowland Modified Rivers (4–5)"
    )),
    line_type = factor(
      if_else(Method == "Extrapolation", "Extrapolated", "Interpolated"),
      levels = c("Interpolated", "Extrapolated")   # Interpolated first in legend
    )
  )

fig_accumulation <- ggplot(sacc_df,
                           aes(x = t, y = qD,
                               colour   = Assemblage,
                               fill     = Assemblage,
                               linetype = line_type)) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL),
              alpha = 0.15, colour = NA, linetype = "solid") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = sacc_colors, name = NULL) +
  scale_fill_manual(values   = sacc_colors, name = NULL) +
  scale_linetype_manual(
    values = c("Interpolated" = "solid", "Extrapolated" = "dashed"),
    name   = NULL
  ) +
  labs(x = "Number of sampling events", y = "Cumulative species richness") +
  guides(
    linetype = guide_legend(
      override.aes = list(linewidth = 1, colour = "grey30"),
      order = 1
    ),
    colour = guide_legend(
      override.aes = list(linewidth = 1),
      order = 2
    ),
    fill = "none"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position      = c(0.98, 0.02),     # bottom right corner
    legend.justification = c("right", "bottom"),
    legend.key.width     = unit(2, "cm"),
    legend.text          = element_text(size = 9)
  )

fig_accumulation


# Run iNEXT by number of samples (replicates, not events) ----------------------------------------
# Uses the full replicate-level P/A matrix — each replicate is one sampling unit
meta_raw <- read_csv("Data/Raw/riparia2024_metadata.csv",  show_col_types = FALSE)

make_inext_samples <- function(min_s, max_s) {
  # Get sample codes for this Strahler range
  samps <- meta_raw %>%
    filter(Strahler >= min_s, Strahler <= max_s) %>%
    pull(Sample)
  mat <- otu_pa %>%
    select(all_of(samps[samps %in% colnames(otu_pa)]))
  # First element = number of samples (T), rest = species incidence counts
  inc <- c(ncol(mat), rowSums(mat))
  inc[inc > 0]
}

inext_out_samples <- iNEXT(
  list(
    "Full basin (Strahler 1–5)"          = make_inext_samples(1, 5),
    "Cardamom headwaters (Strahler 1–3)" = make_inext_samples(1, 3),
    "Lowland modified (Strahler 4–5)"    = make_inext_samples(4, 5)
  ),
  q        = 0,
  datatype = "incidence_freq",
  endpoint = max(length(sample_cols) * 2, 60)
)

print(inext_out_samples$AsyEst)



# Extract tidy data frame for ggplot
sacc_df <- inext_out_samples$iNextEst$size_based %>%
  mutate(
    Assemblage = factor(Assemblage, levels = c(
      "Full basin (Strahler 1–5)",
      "Cardamom headwaters (Strahler 1–3)",
      "Lowland modified (Strahler 4–5)"
    )),
    # Mark interpolated vs extrapolated
    line_type = if_else(Method == "Rarefaction", "Interpolated", "Extrapolated")
  )

sacc_colors <- c(
  "Full Basin (1–5)"              = "darkgreen",
  "Cardamom Headwaters (1–3)"     = "#0073C2FF",
  "Lowland Modified Rivers (4–5)" = "#CD534CFF"
)

sacc_df <- inext_out_samples$iNextEst$size_based %>%
  mutate(
    Assemblage = case_when(
      Assemblage == "Full basin (Strahler 1–5)"          ~ "Full Basin (1–5)",
      Assemblage == "Cardamom headwaters (Strahler 1–3)" ~ "Cardamom Headwaters (1–3)",
      Assemblage == "Lowland modified (Strahler 4–5)"    ~ "Lowland Modified Rivers (4–5)"
    ),
    Assemblage = factor(Assemblage, levels = c(
      "Full Basin (1–5)", "Cardamom Headwaters (1–3)", "Lowland Modified Rivers (4–5)"
    )),
    line_type = factor(
      if_else(Method == "Extrapolation", "Extrapolated", "Interpolated"),
      levels = c("Interpolated", "Extrapolated")   # Interpolated first in legend
    )
  )

fig_accumulation <- ggplot(sacc_df,
                           aes(x = t, y = qD,
                               colour   = Assemblage,
                               fill     = Assemblage,
                               linetype = line_type)) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL),
              alpha = 0.15, colour = NA, linetype = "solid") +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = sacc_colors, name = NULL) +
  scale_fill_manual(values   = sacc_colors, name = NULL) +
  scale_linetype_manual(
    values = c("Interpolated" = "solid", "Extrapolated" = "dashed"),
    name   = NULL
  ) +
  labs(x = "Number of samples collected", y = "Cumulative species richness") +
  guides(
    linetype = guide_legend(
      override.aes = list(linewidth = 1, colour = "grey30"),
      order = 1
    ),
    colour = guide_legend(
      override.aes = list(linewidth = 1),
      order = 2
    ),
    fill = "none"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position      = c(0.98, 0.02),     # bottom right corner
    legend.justification = c("right", "bottom"),
    legend.key.width     = unit(2, "cm"),
    legend.text          = element_text(size = 9)
  )

fig_accumulation

## SP Accu Stats -----
# 7: Coleman species accumulation ----------------------------------------
# Random-placement (Coleman) expectation of richness at a given sample size,
# reported with the Coleman SD. Run for the whole basin and separately for
# the headwater (1–3) and lowland (4–5) groups.

coleman_mat <- function(min_s, max_s) {
  samp <- meta_raw %>% filter(Strahler >= min_s, Strahler <= max_s) %>% pull(Sample)
  samp <- samp[samp %in% colnames(otu_pa)]                  # keep samples present in matrix
  mat  <- t(as.matrix(otu_pa[, samp, drop = FALSE]))        # otu_pa is species × sample → transpose
  mat[, colSums(mat) > 0, drop = FALSE]                     # drop species absent from the subset
}

coleman_report <- function(mat, label, n_small = c(1, 3), ref_n = NULL) {
  per_samp <- rowSums(mat > 0)
  sac      <- specaccum(mat, method = "coleman")
  
  if (is.null(ref_n)) ref_n <- max(sac$sites)
  if (ref_n > max(sac$sites)) {
    warning(sprintf("ref_n (%d) exceeds samples available (%d); using endpoint",
                    ref_n, max(sac$sites)))
    ref_n <- max(sac$sites)
  }
  ref_rich <- sac$richness[which(sac$sites == ref_n)]
  
  cat(sprintf("\n%s  (%d samples, %d species)\n", label, nrow(mat), ncol(mat)))
  cat(sprintf("  Observed per sample: %.2f species (±%.2f SD)\n",
              mean(per_samp), sd(per_samp)))
  for (n in n_small) {
    i <- which(sac$sites == n) 
    cat(sprintf("  Coleman n=%d: %.1f species (±%.1f SD) = %.0f%% of %d-sample richness\n",
                n, sac$richness[i], sac$sd[i],
                100 * sac$richness[i] / ref_rich, ref_n))
  }
  invisible(sac)
}

# ref_n = NULL uses each curve's endpoint; set a number (e.g. 50) for a fixed denominator
coleman_report(coleman_mat(1, 5), "Full basin (1–5)")
coleman_report(coleman_mat(1, 3), "Cardamom headwaters (1–3)")
coleman_report(coleman_mat(4, 5), "Lowland (4–5)")


# Richness Map ------------------------


# 1: Read river shapefile ----------------------------------------
sf_use_s2(FALSE)
rivers <- st_read("Data/Spatial/CardamomRivers.shp")
st_crs(rivers)  # confirm WGS84


# 2: Rotate rivers 60 degrees clockwise around centroid ----------------------------------------
angle <- 60 * pi / 180   # negative = clockwise
geom  <- st_geometry(rivers)
ctr   <- st_centroid(st_union(geom))
geom_rot <- (geom - ctr) * matrix(c(cos(angle), sin(angle),
                                    -sin(angle), cos(angle)), 2, 2) + ctr
st_geometry(rivers) <- geom_rot
st_crs(rivers) <- 4326   # reassign CRS after rotation so north arrow works

# 3: Site summary — richness + coordinates from memory ----------------------------------------
# alpha_div and meta are in memory from 05_alpha_diversity.R
# meta must have Long, Lat, Code columns

site_summary <- alpha_div %>%
  group_by(Site, Strahler, Strahler_f) %>%
  summarise(mean_richness = round(mean(richness), 1), .groups = "drop") %>%
  left_join(meta %>% distinct(Site, Long, Lat, Code), by = "Site") %>%
  filter(!is.na(Long))   # excludes Kravanh and Opreil

# 4: Convert sites to sf and rotate to match rivers ----------------------------------------
sites_sf <- site_summary %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

geom_sites     <- st_geometry(sites_sf)
geom_sites_rot <- (geom_sites - ctr) * matrix(c(cos(angle), sin(angle),
                                                -sin(angle), cos(angle)), 2, 2) + ctr
st_geometry(sites_sf) <- geom_sites_rot
st_crs(sites_sf) <- 4326


# 5: Plot ----------------------------------------

# Add individual nudges to sites_sf — adjust per label
sites_sf <- sites_sf %>%
  mutate(
    nudge_x = case_when(
      Code == "TC1" ~  -0.043,
      Code == "TC2" ~  0.00,
      Code == "TC3" ~  0.00,
      Code == "TS1" ~  -0.02,
      Code == "TS2" ~ -0.020,
      Code == "TS3" ~  0.045,
      Code == "TS4" ~  0.045,
      Code == "SP1" ~  0.00,
      Code == "PS1" ~  0.02,
      TRUE          ~  0.02
    ),
    nudge_y = case_when(
      Code == "TC1" ~  0.0025,
      Code == "TC2" ~ -0.02,
      Code == "TC3" ~  0.02,
      Code == "TS1" ~  0.016,
      Code == "TS2" ~  0.021,
      Code == "TS3" ~  0.00,
      Code == "TS4" ~  0.00,
      Code == "SP1" ~  0.023,
      Code == "PS1" ~  0.02,
      TRUE          ~  0.01
    )
  )



fig_map <- ggplot() +
  geom_sf(data = rivers, aes(colour = "Waterways"), linewidth = 0.4, alpha = 0.8,
          key_glyph = "path") +
  scale_colour_manual(values = c("Waterways" = "black"), name = NULL) +
  
  # Points — shape = Strahler order, fill = richness (continuous colour scale)
  geom_sf(data   = sites_sf,
          aes(shape = Strahler_f, fill = mean_richness),
          size   = 7,
          colour = "black",
          stroke = 0.8,
          alpha  = 1) +
  
  geom_sf_text(
    data = sites_sf,
    aes(label = paste0(Code, " (", round(mean_richness, 1), ")")),
    size = 3,
    fontface = "bold",
    colour = "grey10",
    nudge_x = sites_sf$nudge_x,
    nudge_y = sites_sf$nudge_y
  ) +
  
  # Shapes need to be filled types (21-25) to show fill colour
  scale_shape_manual(values = c("1st" = 24,   # triangle
                                "2nd" = 21,   # circle
                                "3rd" = 22,   # square
                                "4th" = 23,   # diamond (pentagon not available as filled shape)
                                "5th" = 25),  # inverted triangle (octagon not available)
                     name = "Strahler order") +
  scale_fill_viridis_c(name = "Mean species richness", option = "D") +
  
  annotation_north_arrow(
    location    = "br",
    pad_x       = unit(0.8, "cm"),
    pad_y       = unit(1.8, "cm"),
    which_north = "true",
    rotation    = -300,
    style       = north_arrow_fancy_orienteering(
      fill = c("grey40", "white"), line_col = "grey20",
      text_col = "grey20", text_size = 8
    ),
    height = unit(1.5, "cm"), width = unit(1.5, "cm")
  ) +
  
  theme_minimal(base_size = 11) +
  guides(
    shape  = guide_legend(title = "Strahler order", order = 1,
                          override.aes = list(fill = "grey60", size = 4)),
    fill   = guide_colorbar(title = "Mean species richness", order = 2),
    colour = guide_legend(title = NULL, order = 3,
                          override.aes = list(
                            linetype  = 1, linewidth = 0.8,
                            shape = NA, fill = NA
                          ))
  ) +
  theme(
    legend.box            = "horizontal",
    legend.position       = c(0.3, 0.02),
    legend.justification  = c("left", "bottom"),
    legend.title          = element_text(size = 10, face = "bold"),
    legend.text           = element_text(size = 10),
    legend.key.size       = unit(0.4, "cm"),
    panel.background      = element_rect(fill = "white", colour = NA)
  )

fig_map

