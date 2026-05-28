# =============================================================================
# format_local_refs_apscale.R
#
# Converts local Sanger reference sequences (exported from Geneious) into
# APSCALE/MIDORI-compatible FASTA format for eDNA metabarcoding pipelines.
#
# WORKFLOW:
#   1. Export sequences from Geneious as FASTA
#      Headers must look like: >CC04.Lipinia_cf_vittigera
#   2. Set paths in CONFIG below
#   3. Run whole script — queries NCBI once, caches results, writes review table
#   4. Check review_table.csv, fix any flagged rows, re-run from STEP 4 only
#
# OUTPUT HEADER FORMAT:
#   >CC04.Lipinia_cf_vittigera###root_1;Eukaryota_2759;Chordata_7711;
#     Lepidosauria_8504;Squamata_8509;Scincidae_66056;Lipinia_71016;Lipinia_sp_71016
#
# REQUIREMENTS:
#   R packages: ape, stringr, taxize, dplyr, seqinr, readr
#   NCBI API key recommended — get one free at https://www.ncbi.nlm.nih.gov/account/
#     Run usethis::edit_r_environ() and add: ENTREZ_KEY="your_key_here"
#
# AUTHORS: Christian Ching — Cardamoms eDNA Metabarcoding Project
# =============================================================================

# ----------------------------------------
# Setup
library(ape)
library(stringr)
library(taxize)
library(dplyr)
library(seqinr)
library(readr)

# CONFIG ----------------------------------------
input_fasta  <- "databases/sanger_databases/16S_database.fasta"
cache_file   <- "databases/sanger_databases/ncbi_cache.csv"
review_file  <- "databases/sanger_databases/review_table.csv"
output_fasta <- "databases/sanger_databases/2025cardamoms16S_database_midorimode.fasta"

# How to handle cf. identifications (e.g. Lipinia_cf_vittigera):
#   "species" = strip cf. and assign to species — Lipinia vittigera
#   "genus"   = treat as genus-level only — Lipinia sp.
cf_handling <- "species"


# STEP 1: Read and validate input FASTA ----------------------------------------

cat("\nReading input FASTA...\n")

if (!file.exists(input_fasta)) {
  stop("Input file not found:\n  ", input_fasta,
       "\nCheck your path in CONFIG.")
}

fasta_seqs <- read.fasta(
  file           = input_fasta,
  seqtype        = "DNA",
  as.string      = TRUE,
  set.attributes = FALSE
)

cat(sprintf("  %d sequences loaded.\n", length(fasta_seqs)))

# Check all headers have the expected . delimiter
bad_headers <- names(fasta_seqs)[!str_detect(names(fasta_seqs), "\\.")]
if (length(bad_headers) > 0) {
  warning(sprintf(
    "%d headers missing . delimiter — these will be skipped:\n%s",
    length(bad_headers), paste(" -", bad_headers, collapse = "\n")
  ))
}

# Parse sample codes and species names from headers
header_df <- data.frame(original_header = names(fasta_seqs),
                        stringsAsFactors = FALSE) %>%
  filter(str_detect(original_header, "\\.")) %>%
  mutate(
    sample_code  = str_extract(original_header, "^[^.]+"),
    species_raw  = str_extract(original_header, "(?<=\\.).+"),
    is_cf        = str_detect(species_raw, "_cf_"),
    is_sp        = str_detect(species_raw, "_sp$"),
    lookup_name  = case_when(
      is_cf & cf_handling == "genus"    ~ str_extract(species_raw, "^[^_]+"),
      is_cf & cf_handling == "species"  ~ str_replace_all(str_replace(species_raw, "_cf_", "_"), "_", " "),
      is_sp                             ~ str_extract(species_raw, "^[^_]+"),
      TRUE                              ~ str_replace_all(species_raw, "_", " ")
    ),
    target_level = case_when(
      is_sp                           ~ "genus",
      is_cf & cf_handling == "genus"  ~ "genus",
      TRUE                            ~ "species"
    )
  )

cat(sprintf("  Parsed: %d sequences (%d cf., %d sp., %d species-level)\n",
            nrow(header_df), sum(header_df$is_cf),
            sum(header_df$is_sp), sum(header_df$target_level == "species")))


# STEP 2: Query NCBI (with caching) ----------------------------------------
# Results saved to cache_file after first run.
# Re-running is always safe — only new sequences are queried.
# If NCBI times out mid-run, just re-run: cache preserves progress.
#
# NOTE: taxize 0.10.0 returns Eukaryota with rank "domain" not "superkingdom".
# rank_preferences handles this by accepting multiple rank names per position.

rank_preferences <- list(
  superkingdom = c("domain", "superkingdom", "kingdom"),
  phylum       = c("phylum", "subphylum"),
  class        = c("class", "superclass", "subclass"),
  order        = c("order", "superorder", "suborder"),
  family       = c("family", "superfamily", "subfamily"),
  genus        = c("genus", "subgenus"),
  species      = c("species")
)

pick_ranks <- function(lin, rank_preferences) {
  result <- lapply(names(rank_preferences), function(position) {
    prefs <- rank_preferences[[position]]
    for (pref in prefs) {
      row <- lin[lin$rank == pref, ]
      if (nrow(row) > 0) {
        return(data.frame(name = row$name[1], rank = position,
                          id   = row$id[1], stringsAsFactors = FALSE))
      }
    }
    return(NULL)
  })
  bind_rows(result)
}

fetch_ncbi <- function(lookup_name, target_level) {
  
  safe_get_uid <- function(name) {
    for (attempt in 1:3) {
      result <- tryCatch(
        suppressWarnings(get_uid(name, messages = FALSE)[[1]]),
        error = function(e) { Sys.sleep(2 * attempt); NA }
      )
      if (!is.na(result)) return(result)
      Sys.sleep(1)
    }
    return(NA)
  }
  
  safe_classification <- function(uid) {
    for (attempt in 1:3) {
      result <- tryCatch({
        res <- classification(uid, db = "ncbi")
        if (is.list(res) && length(res) > 0 && is.data.frame(res[[1]])) {
          res[[1]]
        } else {
          NULL
        }
      },
      error   = function(e) { Sys.sleep(2 * attempt); NULL },
      warning = function(w) { Sys.sleep(1); NULL }
      )
      if (!is.null(result) && is.data.frame(result) && nrow(result) > 0) {
        return(result)
      }
      Sys.sleep(1)
    }
    return(NULL)
  }
  
  uid    <- safe_get_uid(lookup_name)
  status <- "ok"
  notes  <- ""
  
  if (is.na(uid)) {
    genus <- str_split(lookup_name, " ")[[1]][1]
    uid   <- safe_get_uid(genus)
    if (!is.na(uid)) {
      status <- "genus_fallback"
      notes  <- sprintf("'%s' not in NCBI; fell back to genus '%s'", lookup_name, genus)
    } else {
      return(list(taxonomy_string = NA, status = "failed",
                  notes = sprintf("Could not find '%s' or its genus in NCBI", lookup_name)))
    }
  }
  
  lin <- safe_classification(uid)
  
  if (is.null(lin) || !is.data.frame(lin) || nrow(lin) == 0) {
    return(list(taxonomy_string = NA, status = "failed",
                notes = sprintf("Empty classification returned for uid %s", uid)))
  }
  
  lin_filtered <- pick_ranks(lin, rank_preferences)
  
  sk <- lin_filtered %>% filter(rank == "superkingdom") %>% pull(name)
  if (length(sk) == 0 || !sk[1] %in% c("Eukaryota", "Eukaryotes")) {
    return(list(taxonomy_string = NA, status = "failed",
                notes = "Non-eukaryote or missing superkingdom — check NCBI hit"))
  }
  
  if (target_level == "genus" || !("species" %in% lin_filtered$rank)) {
    genus_row <- lin_filtered %>% filter(rank == "genus")
    if (nrow(genus_row) == 0) {
      return(list(taxonomy_string = NA, status = "failed",
                  notes = "No genus node found in NCBI lineage"))
    }
    sp_row <- data.frame(
      name = paste0(genus_row$name[1], "_sp"),
      rank = "species",
      id   = genus_row$id[1],
      stringsAsFactors = FALSE
    )
    lin_filtered <- bind_rows(lin_filtered, sp_row)
    if (status == "ok") notes <- "Genus-level: synthetic species node added"
  }
  
  lin_ordered <- lin_filtered %>%
    mutate(rank = factor(rank, levels = names(rank_preferences))) %>%
    arrange(rank)
  
  if (nrow(lin_ordered) != 7) {
    missing_r <- setdiff(names(rank_preferences), as.character(lin_ordered$rank))
    notes <- paste0(notes, sprintf(" | WARNING: %d/7 nodes — missing: %s",
                                   nrow(lin_ordered), paste(missing_r, collapse = ", ")))
  }
  
  nodes      <- lin_ordered %>% mutate(label = paste0(name, "_", id)) %>% pull(label)
  tax_string <- paste(c("root_1", nodes), collapse = ";")
  
  list(taxonomy_string = tax_string, status = status, notes = str_trim(notes))
}

# Load or initialise cache
if (file.exists(cache_file)) {
  cache <- read_csv(cache_file, show_col_types = FALSE)
  cat(sprintf("\nCache loaded: %d previously queried names.\n", nrow(cache)))
} else {
  cache <- data.frame(lookup_name = character(), target_level = character(),
                      taxonomy_string = character(), status = character(),
                      notes = character(), stringsAsFactors = FALSE)
  cat("\nNo cache found — querying all sequences from NCBI.\n")
}

to_query <- header_df %>%
  select(lookup_name, target_level) %>%
  distinct() %>%
  anti_join(cache, by = c("lookup_name", "target_level"))

if (nrow(to_query) == 0) {
  cat("All sequences already cached — skipping NCBI queries.\n")
} else {
  cat(sprintf("Querying NCBI for %d unique names", nrow(to_query)))
  if (Sys.getenv("ENTREZ_KEY") == "") cat(" (set ENTREZ_KEY in .Renviron for faster queries)")
  cat("...\n")
  
  new_results <- vector("list", nrow(to_query))
  for (i in seq_len(nrow(to_query))) {
    nm  <- to_query$lookup_name[i]
    lvl <- to_query$target_level[i]
    cat(sprintf("  [%d/%d] %-40s", i, nrow(to_query), paste0(nm, " (", lvl, ")")))
    res <- fetch_ncbi(nm, lvl)
    new_results[[i]] <- data.frame(
      lookup_name     = nm,
      target_level    = lvl,
      taxonomy_string = ifelse(is.null(res$taxonomy_string), NA, res$taxonomy_string),
      status          = res$status,
      notes           = res$notes,
      stringsAsFactors = FALSE
    )
    cat(sprintf("[%s]\n", res$status))
    Sys.sleep(0.34)
  }
  
  cache <- bind_rows(cache, bind_rows(new_results))
  write_csv(cache, cache_file)
  cat(sprintf("Cache saved to: %s\n", cache_file))
}


# STEP 3: Build and save review table ----------------------------------------
# Open review_table.csv in Excel and check anything flagged other than "OK".
#
# FLAG meanings:
#   OK             — No action needed
#   GENUS_FALLBACK — Species not in NCBI; genus used instead. Check spelling
#                    or whether species is too recently described for NCBI.
#   NODE_WARNING   — Taxonomy has fewer than 8 nodes. See notes column.
#   FAILED         — Lookup failed. Either:
#                    a) Fix the name, delete from cache_file, re-run; OR
#                    b) Manually enter taxonomy_string in the CSV,
#                       change flag to MANUAL, re-run from STEP 4.
#
# MANUAL TAXONOMY STRING FORMAT (8 nodes, semicolon-separated):
#   root_1;Eukaryota_2759;Chordata_7711;Amphibia_8292;Anura_8342;
#   Ranidae_8397;Sylvirana_1659747;Sylvirana_mortenseni_1659781
#   Find taxon IDs at: https://www.ncbi.nlm.nih.gov/taxonomy

review <- header_df %>%
  left_join(cache, by = c("lookup_name", "target_level")) %>%
  mutate(
    apscale_header = case_when(
      !is.na(taxonomy_string) ~ paste0(sample_code, ".", species_raw, "###", taxonomy_string),
      TRUE ~ NA_character_
    ),
    flag = case_when(
      is.na(taxonomy_string)       ~ "FAILED",
      status == "genus_fallback"   ~ "GENUS_FALLBACK",
      str_detect(notes, "WARNING") ~ "NODE_WARNING",
      TRUE                         ~ "OK"
    )
  ) %>%
  select(sample_code, species_raw, is_cf, is_sp, target_level,
         lookup_name, status, flag, notes, taxonomy_string, apscale_header, original_header)

write_csv(review, review_file)

cat(sprintf("\n--- Review table summary ---\n"))
cat(sprintf("  OK:              %d\n", sum(review$flag == "OK")))
cat(sprintf("  Genus fallback:  %d\n", sum(review$flag == "GENUS_FALLBACK")))
cat(sprintf("  Node warning:    %d\n", sum(review$flag == "NODE_WARNING")))
cat(sprintf("  FAILED:          %d\n", sum(review$flag == "FAILED")))
cat(sprintf("\nReview table saved to: %s\n", review_file))

if (sum(review$flag == "FAILED") > 0) {
  cat(sprintf("\n*** ACTION REQUIRED: %d sequences failed. ***\n", sum(review$flag == "FAILED")))
  cat("Fix in review_table.csv then re-run from STEP 4.\n\nFailed sequences:\n")
  review %>% filter(flag == "FAILED") %>% select(sample_code, species_raw, notes) %>% print()
}


# STEP 4: Write output FASTA ----------------------------------------
# Re-reads review_table.csv so any manual edits are picked up.
# Re-run from here after editing review_table.csv — no need to re-query NCBI.

cat("\nReading review table (picks up any manual edits)...\n")
review_final <- read_csv(review_file, show_col_types = FALSE)

to_write <- review_final %>% filter(!is.na(taxonomy_string), !is.na(apscale_header))
excluded <- review_final %>% filter(is.na(taxonomy_string) | is.na(apscale_header))

cat(sprintf("  Writing:  %d sequences\n", nrow(to_write)))
if (nrow(excluded) > 0) {
  cat(sprintf("  Excluded: %d sequences:\n", nrow(excluded)))
  cat(paste(" -", excluded$sample_code, excluded$species_raw, collapse = "\n"), "\n")
}

out_seqs    <- fasta_seqs[to_write$original_header]
out_headers <- str_remove(to_write$apscale_header, "^>")

# Write uppercase sequences wrapped at 60 chars per line — matching MIDORI format
wrap_sequence <- function(seq_str, width = 60) {
  seq_upper <- toupper(seq_str)
  starts    <- seq(1, nchar(seq_upper), by = width)
  chunks    <- substring(seq_upper, starts, starts + width - 1)
  paste(chunks, collapse = "\n")
}

con <- file(output_fasta, open = "w")
for (i in seq_along(out_seqs)) {
  writeLines(paste0(">", out_headers[i]), con)
  writeLines(wrap_sequence(out_seqs[[i]]), con)
}
close(con)

# Validate node counts
node_counts <- str_count(out_headers, ";") + 1
if (all(node_counts == 8)) {
  cat(sprintf("Node count check PASSED: all %d headers have exactly 8 nodes.\n", nrow(to_write)))
} else {
  bad <- out_headers[node_counts != 8]
  cat(sprintf("WARNING: %d headers do not have 8 nodes:\n", length(bad)))
  cat(paste(" -", bad, collapse = "\n"), "\n")
}

# Spot check output format
cat("\nOutput format check (first sequence):\n")
cat(paste(readLines(output_fasta, n = 3), collapse = "\n"), "\n")

cat(sprintf("\nDone. %d sequences written to:\n  %s\n", nrow(to_write), output_fasta))
if (nrow(excluded) > 0) cat(sprintf("  (%d excluded)\n", nrow(excluded)))
cat("\nNext step: run 0_2_merge_fastas_build_apscale_db.py\n")