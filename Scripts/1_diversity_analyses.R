# Load Libraries
library(tidyverse)
library(vegan)
library(reshape2)
library(ggplot2)

### Setup ---------------------------------------------------------------------------------
# Load Data
read_table <- read.csv("Data/shaw16S_final_read_table.csv")
metadata   <- read.csv("Data/metadata.csv")

# Check Data
head(read_table)
str(read_table)

head(metadata)
str(metadata)

# Prepare Read_Table columns
taxonomy_cols <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)

# Everything after 'Seq' are sample columns
# Find the column index of 'Seq'
seq_col <- which(names(read_table) == "Seq")

# Sample columns are everything AFTER 'Seq'
sample_cols <- names(read_table)[(seq_col + 1):ncol(read_table)]


# Collapse ASVs into OTUs 
otu_table <- read_table %>%
  # Keep only taxonomy + sample columns
  select(all_of(taxonomy_cols), all_of(sample_cols)) %>%
  # Group by taxonomy (Species defines the OTU)
  group_by(across(all_of(taxonomy_cols))) %>%
  # Sum reads across ASVs for each species
  summarise(
    across(
      all_of(sample_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Inspect result
head(otu_table)
dim(otu_table)


# To join metadata:
taxonomy_cols <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)

# Convert OTU table to long-form
otu_long <- otu_table %>%
  pivot_longer(
    cols = -all_of(taxonomy_cols),
    names_to = "Sample",
    values_to = "Reads"
  )

# Attach Metadata
otu_meta <- otu_long %>%
  left_join(metadata, by = "Sample")
# Remove 0 read rows
otu_meta <- otu_meta %>%
  filter(Reads > 0)

head(otu_meta)
str(otu_meta)


