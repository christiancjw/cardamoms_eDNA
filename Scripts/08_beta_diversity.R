# ----------------------------------------------------------------
# 08_beta_diversity.R
#
# Beta diversity analysis — community composition across the
# hydrological gradient and seasons.
#
# ANALYSES:
#   1. Jaccard dissimilarity matrix (P/A data)
#   2. NMDS ordination (k=2)
#   3. PERMANOVA: Strahler + Season + QBR
#   4. PERMDISP: test homogeneity of dispersion
#   5. Seasonal beta per site (Jaccard wet vs dry)
#   6. Family/order composition bar charts by site
#
# FIGURES:
#   fig_nmds_strahler    — NMDS coloured by Strahler order
#   fig_nmds_season      — NMDS coloured by season
#   fig_nmds_combined    — cowplot panel A+B
#   fig_seasonal_beta    — per-site seasonal Jaccard by Strahler
#   fig_composition      — stacked bar family composition by site
#


#  Setup ----------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(vegan)
library(tidyr)
library(cowplot)
library(cluster)

# Palettes — consistent with other scripts
strahler_f_colors <- c(
  "1st" = "#C6DCEF", "2nd" = "#84B9D9",
  "3rd" = "#4393C3", "4th" = "#1E6FA8", "5th" = "#084A87"
)

qbr_colors <- c(
  "Natural" = "#006BA6", "Good" = "#22A122",
  "Fair"    = "#DBCB43", "Poor" = "#FF7134", "Bad" = "#AF3245"
)

season_colors <- c("Dry" = "#D4A94A", "Wet" = "#4A90D4")

site_order <- c(
  "Knaong Bat Sa", "Stung Ta Chey", "Stung Ta Say", "Arai Village",
  "Kronomh Village", "Stung Oda", "Dam 5", "Tang Rang Bridge", "Pursat River"
)

code_order <- c(
  "TC1", "TC2", "TC3", "TS1", "TS2", "TS3", "TS4", "SP1", "PS1"
)

code_colors <- c(
  "TC1" = "#A8D8B9",   # light green
  "TC2" = "#4FA873",   # mid green
  "TC3" = "#1B6B3A",   # dark green
  
  "TS1" = "#A8C9E8",   # light blue
  "TS2" = "#5E9BD4",   # mid blue
  "TS3" = "#2C5F9E",   # darker blue
  
  "TS4" = "#98319B",   # deepest blue
  
  "SP1" = "#F2A35C",   # orange
  "PS1" = "#C0392B"    # red-orange
)
## scale_fill_manual(values = site_colors, name = "Site")

# Add the highland vs lowlands 


# 1: Read in data (EVENTS) ----------------------------------------
pa   <- read_csv("Data/Clean/fish_community_pa.csv",   show_col_types = FALSE)
meta <- read_csv("Data/Clean/fish_community_meta.csv", show_col_types = FALSE) %>%
  mutate(
    Site       = factor(Site, levels = site_order),
    Strahler_f = factor(Strahler, levels = 1:5,
                        labels = c("1st","2nd","3rd","4th","5th"), ordered = TRUE),
    Season     = factor(Season,    levels = c("Dry","Wet")),
    QBR_Class  = factor(QBR_Class, levels = c("Natural","Good","Fair","Poor","Bad"))
  )

identical(pa$Event, meta$Event)  # must be TRUE

taxa_cols <- colnames(pa)[colnames(pa) != "Event"]
otu_mat   <- pa %>% select(all_of(taxa_cols)) %>% as.matrix()
rownames(otu_mat) <- pa$Event

# Read in SAMPLES
pa_sample   <- read_csv("Data/Clean/fish_community_pa_sample.csv",   show_col_types = FALSE)
meta_sample <- read_csv("Data/Clean/fish_community_meta_sample.csv", show_col_types = FALSE) %>%
  mutate(
    Site       = factor(Site, levels = site_order),
    Strahler_f = factor(Strahler, levels = 1:5,
                        labels = c("1st","2nd","3rd","4th","5th"), ordered = TRUE),
    Season     = factor(Season,    levels = c("Dry","Wet")),
    QBR_Class  = factor(QBR_Class, levels = c("Natural","Good","Fair","Poor","Bad"))
  )

identical(pa_sample$Sample, meta_sample$Sample)  # must be TRUE

taxa_cols_sample <- colnames(pa_sample)[colnames(pa_sample) != "Sample"]
otu_mat_sample   <- pa_sample %>% select(all_of(taxa_cols_sample)) %>% as.matrix()
rownames(otu_mat_sample) <- pa_sample$Sample

# 2: Family composition bar chart ----------------------------------------
# Shows what communities look like before collapsing to distance matrix
# Sites ordered by Strahler, families < 10% of events collapsed to Other

species_list <- read_csv("Data/Clean/fish_species_list.csv", show_col_types = FALSE)

family_site <- pa %>%
  pivot_longer(-Event, names_to = "taxon_label", values_to = "presence") %>%
  filter(presence == 1) %>%
  left_join(species_list %>% select(taxon_label, Family), by = "taxon_label") %>%
  left_join(meta %>% select(Event, Site, Strahler_f), by = "Event") %>%
  filter(!is.na(Family)) %>%
  group_by(Site, Strahler_f, Family) %>%
  summarise(n_events = n_distinct(Event), .groups = "drop") %>%
  left_join(meta %>% dplyr::count(Site, name = "total_events"), by = "Site") %>%
  mutate(
    prop   = n_events / total_events,
    Family = if_else(prop < 0.1, "Other", Family)
  ) %>%
  group_by(Site, Strahler_f, Family) %>%
  summarise(prop = sum(prop), .groups = "drop") %>%
  mutate(Site = factor(Site, levels = site_order))

fig_composition <- ggplot(family_site,
                          aes(x = Site, y = prop, fill = Family)) +
  geom_col(position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Proportion of events detected") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank()
  )

fig_composition


# 3: Jaccard dissimilarity matrix ----------------------------------------
jac_dist <- vegdist(otu_mat, method = "jaccard", binary = TRUE)

mean(jac_dist)
range(jac_dist)

jac_dist_sample <- vegdist(otu_mat_sample, method = "jaccard", binary = TRUE)

mean(jac_dist_sample)
range(jac_dist_sample)

# 4: NMDS ----------------------------------------
set.seed(123)
nmds <- metaMDS(jac_dist, k = 2, trymax = 100, trace = FALSE)
nmds$stress   # < 0.10 good, < 0.20 acceptable

# Create NMDS Dataframe
nmds_df <- as.data.frame(scores(nmds, display = "sites")) %>%
  tibble::rownames_to_column("Event") %>%
  left_join(meta, by = "Event")

# Stress label for annotation
stress_label <- sprintf("Stress = %.3f", nmds$stress)

# Panel A — Strahler order
fig_nmds_strahler <- ggplot(nmds_df,
                            aes(x = NMDS1, y = NMDS2,
                                colour = Code, fill = Code)) +
  stat_ellipse(geom = "polygon", alpha = 0.08, linewidth = 0.4, level = 0.95) +
  geom_point(size = 3, alpha = 0.9) +
  annotate("text", x = Inf, y = -Inf, label = stress_label,
           hjust = 1.1, vjust = -0.5, size = 3, colour = "grey50") +
  scale_colour_manual(values = code_colors, name = NULL) +
  scale_fill_manual(values = code_colors, name = "Code") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

fig_nmds_strahler

`# Panel B — Season
fig_nmds_season <- ggplot(nmds_df,
                          aes(x = NMDS1, y = NMDS2,
                              colour = Season, fill = Season)) +
  stat_ellipse(geom = "polygon", alpha = 0.08, linewidth = 0.4,
               level = 0.95, aes(group = Season)) +
  geom_point(aes(shape = Strahler_f), size = 3, alpha = 0.9) +
  scale_colour_manual(values = season_colors, name = NULL) +
  scale_fill_manual(values   = season_colors, guide = "none") +
  scale_shape_manual(values  = c("1st"=16,"2nd"=17,"3rd"=18,"4th"=15,"5th"=8),
                     name = NULL) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Panel C — QBR class
fig_nmds_qbr <- ggplot(nmds_df,
                       aes(x = NMDS1, y = NMDS2,
                           colour = QBR_Class, fill = QBR_Class)) +
  stat_ellipse(geom = "polygon", alpha = 0.08, linewidth = 0.4, level = 0.95) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = qbr_colors, name = NULL) +
  scale_fill_manual(values   = qbr_colors, guide = "none") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# View individually
fig_nmds_strahler
fig_nmds_season
fig_nmds_qbr

# Combined panel
fig_nmds_combined <- plot_grid(
  fig_nmds_strahler + theme(legend.position = "none"),
  fig_nmds_season   + theme(legend.position = "none"),
  fig_nmds_qbr      + theme(legend.position = "none"),
  ncol = 3, labels = c("A","B","C"), label_size = 12
)

fig_nmds_combined

### NMDS SAMPLE ----------------

jac_dist_sample <- vegdist(otu_mat_sample, method = "jaccard", binary = TRUE)

mean(jac_dist_sample)
range(jac_dist_sample)

set.seed(123)
nmds_sample <- metaMDS(jac_dist_sample, k = 2, trymax = 100, trace = FALSE)
nmds_sample$stress   # < 0.10 good, < 0.20 acceptable

# Create NMDS Dataframe
nmds_sample_df <- as.data.frame(scores(nmds_sample, display = "sites")) %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(meta_sample, by = "Sample")

# Stress label for annotation
sample_stress_label <- sprintf("Stress = %.3f", nmds_sample$stress)

# Panel A — Strahler order
fig_nmds_strahler <- ggplot(nmds_sample_df,
                            aes(x = NMDS1, y = NMDS2,
                                colour = Code, fill = Code)) +
  stat_ellipse(geom = "polygon", alpha = 0.08, linewidth = 0.4, level = 0.95) +
  geom_point(size = 3, alpha = 0.9) +
  annotate("text", x = Inf, y = -Inf, label = sample_stress_label,
           hjust = 1.1, vjust = -0.5, size = 3, colour = "grey50") +
  scale_colour_manual(values = code_colors, name = NULL) +
  scale_fill_manual(values = code_colors, name = "Code") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

fig_nmds_strahler

`# Panel B — Season
fig_nmds_season <- ggplot(nmds_df,
                          aes(x = NMDS1, y = NMDS2,
                              colour = Season, fill = Season)) +
  stat_ellipse(geom = "polygon", alpha = 0.08, linewidth = 0.4,
               level = 0.95, aes(group = Season)) +
  geom_point(aes(shape = Strahler_f), size = 3, alpha = 0.9) +
  scale_colour_manual(values = season_colors, name = NULL) +
  scale_fill_manual(values   = season_colors, guide = "none") +
  scale_shape_manual(values  = c("1st"=16,"2nd"=17,"3rd"=18,"4th"=15,"5th"=8),
                     name = NULL) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Panel C — QBR class
fig_nmds_qbr <- ggplot(nmds_df,
                       aes(x = NMDS1, y = NMDS2,
                           colour = QBR_Class, fill = QBR_Class)) +
  stat_ellipse(geom = "polygon", alpha = 0.08, linewidth = 0.4, level = 0.95) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = qbr_colors, name = NULL) +
  scale_fill_manual(values   = qbr_colors, guide = "none") +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# View individually
fig_nmds_strahler
fig_nmds_season
fig_nmds_qbr

# Combined panel
fig_nmds_combined <- plot_grid(
  fig_nmds_strahler + theme(legend.position = "none"),
  fig_nmds_season   + theme(legend.position = "none"),
  fig_nmds_qbr      + theme(legend.position = "none"),
  ncol = 3, labels = c("A","B","C"), label_size = 12
)

fig_nmds_combined
# 5: Seasonal beta per site ----------------------------------------
# Jaccard between dry and wet season centroid per site
# One value per site — shows which sites shift most seasonally

seasonal_beta <- meta %>%
  group_by(Site, Strahler_f) %>%
  filter(n_distinct(Season) == 2) %>%   # only sites with both seasons
  group_split() %>%
  lapply(function(s) {
    dry <- s %>% filter(Season == "Dry") %>% pull(Event)
    wet <- s %>% filter(Season == "Wet") %>% pull(Event)
    dry_c <- colMeans(otu_mat[dry, , drop = FALSE])
    wet_c <- colMeans(otu_mat[wet, , drop = FALSE])
    jac   <- as.numeric(vegdist(rbind(dry_c, wet_c), method = "jaccard"))
    data.frame(Site = s$Site[1], Strahler_f = s$Strahler_f[1], seasonal_jaccard = jac)
  }) %>%
  bind_rows() %>%
  mutate(Site = factor(Site, levels = site_order))

seasonal_beta

fig_seasonal_beta <- ggplot(seasonal_beta,
                            aes(x = Strahler_f, y = seasonal_jaccard,
                                colour = Strahler_f, label = Site)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(nudge_y = 0.02, size = 3, colour = "grey30") +
  scale_colour_manual(values = strahler_f_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Strahler order", y = "Seasonal Jaccard (dry vs wet)") +
  theme_classic(base_size = 12)

fig_seasonal_beta


# 6: Clustering ----------------------------------------

# Silhouette width to choose optimal number of clusters (k)
# Higher = better defined clusters
hclust_ward <- hclust(jac_dist, method = "ward.D2")

sil_widths <- sapply(2:6, function(k) {
  mean(silhouette(cutree(hclust_ward, k), jac_dist)[, 3])
})

plot(2:6, sil_widths, type = "b",
     xlab = "Number of clusters (k)", ylab = "Mean silhouette width",
     pch = 16)

# Choose k based on silhouette peak — change k here
k <- 3
meta$cluster <- factor(cutree(hclust_ward, k))

# Dendrogram with cluster boxes
plot(hclust_ward, labels = meta$Site, hang = -1, cex = 0.7,
     main = paste("Ward clustering — k =", k))
rect.hclust(hclust_ward, k = k, border = "red")

# Overlay clusters on NMDS
nmds_df <- nmds_df %>% left_join(meta %>% select(Event, cluster), by = "Event")

fig_nmds_cluster <- ggplot(nmds_df,
                           aes(x = NMDS1, y = NMDS2,
                               colour = cluster, shape = Strahler_f)) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_shape_manual(values = c("1st"=16,"2nd"=17,"3rd"=18,"4th"=15,"5th"=8),
                     name = "Strahler order") +
  labs(x = "NMDS1", y = "NMDS2", colour = "Cluster") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

fig_nmds_cluster

# Contingency tables — do clusters align with gradients?
# Rows = cluster, columns = grouping variable
table(meta$cluster, meta$Strahler_f)   # vs Strahler
table(meta$cluster, meta$Season)       # vs Season
table(meta$cluster, meta$QBR_Class)    # vs QBR
table(meta$cluster, meta$Branch)


# Clustering Sp?
# For each cluster, what proportion of events detected each species?
cluster_species <- otu_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Event") %>%
  left_join(meta %>% select(Event, cluster), by = "Event") %>%
  group_by(cluster) %>%
  summarise(across(all_of(taxa_cols), mean)) %>%
  pivot_longer(-cluster, names_to = "Species", values_to = "prop_detected")

# Top 10 species per cluster
cluster_species %>%
  group_by(cluster) %>%
  slice_max(prop_detected, n = 10) %>%
  print(n = 30)

# Species detected in > 50% of events in one cluster but < 20% in others
cluster_wide <- cluster_species %>%
  pivot_wider(names_from = cluster, values_from = prop_detected,
              names_prefix = "cluster_")

# Species strongly associated with Cluster 1 (headwaters)
cluster_wide %>% filter(cluster_1 > 0.5, cluster_2 < 0.2, cluster_3 < 0.2) %>%
  arrange(desc(cluster_1))

# Species strongly associated with Cluster 3 (lowland)
cluster_wide %>% filter(cluster_3 > 0.5, cluster_1 < 0.2, cluster_2 < 0.2) %>%
  arrange(desc(cluster_3))

## Sp occurance ranking------
data.frame(
  species   = taxa_cols,
  n_events  = colSums(otu_mat),
  pct_events = round(100 * colSums(otu_mat) / nrow(otu_mat), 1)
) %>% arrange(desc(n_events))

# Betapart decomposition -------
library(betapart)
beta_pair <- beta.pair(otu_mat, index.family = "jaccard")
# beta_pair$beta.jtu = turnover
# beta_pair$beta.jne = nestedness

# Indicator Sp
library(indicspecies)
multipatt(otu_mat, meta$Strahler_f, func = "IndVal.g", control = how(nperm = 999))
multipatt(otu_mat, meta$Branch, func = "IndVal.g", control = how(nperm = 999))

##  interactive NMDS --------
library(plotly)

# Add a hover label combining all useful info
nmds_df <- nmds_df %>%
  mutate(hover = paste0(
    "Event: ",     Event,     "<br>",
    "Site: ",      Site,      "<br>",
    "Strahler: ",  Strahler,  "<br>",
    "Season: ",    Season,    "<br>",
    "Cluster: ",   cluster,   "<br>",
    "QBR: ",       QBR_Class
  ))

p_interactive <- plot_ly(
  data   = nmds_df,
  x      = ~NMDS1,
  y      = ~NMDS2,
  color  = ~cluster,
  symbol = ~Strahler_f,
  text   = ~hover,
  hoverinfo = "text",
  type   = "scatter",
  mode   = "markers",
  marker = list(size = 10, opacity = 0.9)
) %>%
  layout(
    title  = paste("NMDS — k =", k, "clusters"),
    xaxis  = list(title = "NMDS1"),
    yaxis  = list(title = "NMDS2"),
    legend = list(orientation = "v")
  )

p_interactive

