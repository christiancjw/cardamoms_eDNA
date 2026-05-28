"""
fix_nomatch.py

Post-processes APSCALE blast taxonomy output to override NoMatch assignments
for OTUs that have a high-similarity BLAST hit but were assigned NoMatch due
to LCA conflict resolution.

LOGIC:
    For each NoMatch OTU with similarity > 0:
        1. Find all BLAST hits for that OTU across all subset CSVs
        2. Take the top hit (highest similarity)
        3. If top hit similarity >= override_threshold (default 97%):
           - Parse taxonomy directly from the top hit sseqid
           - Override the NoMatch with that taxonomy
           - Flag the row as "Manual override - top hit"

WHEN IS THIS JUSTIFIED:
    A 97-100% match to a voucher barcode from your own fieldwork is an
    unambiguous identification. APSCALE's NoMatch in these cases is an
    artefact of the LCA algorithm when hits to distantly-related taxa are
    included in the conflict set. The override is scientifically defensible
    and should be documented in your methods.

USAGE:
    Edit the CONFIG section below, then:
        conda activate apscale4
        python fix_nomatch.py

AUTHORS: Christian Ching — Cardamoms eDNA Metabarcoding Project
"""

import os
import re
import glob
import pandas as pd
from pathlib import Path


# =============================================================================
# CONFIG
# =============================================================================

# Final taxonomy file from APSCALE blast
taxonomy_file = Path("/Users/christianching/anaconda3/envs/apscale4/APSCALE_projects/24shaw16spaired_apscale/11_read_table/data/blastn_0_24shaw16spaired_sequences/blastn_0_24shaw16spaired_sequences_taxonomy.xlsx")

# Folder containing all subset_X_blastn.csv files
blast_folder = Path("/Users/christianching/anaconda3/envs/apscale4/APSCALE_projects/24shaw16spaired_apscale/11_read_table/data/blastn_0_24shaw16spaired_sequences/subsets")

# Output file
output_file = taxonomy_file.parent / (taxonomy_file.stem + "_fixed.xlsx")

# Similarity threshold for override (%)
# OTUs with top hit >= this value will be overridden if currently NoMatch
override_threshold = 97.0

# If True, only override when the top hit is a LOCAL sequence (See local_prefixes below)
# If False, override using any top hit >= threshold
local_only = False

# Local sequence prefixes (used when local_only = True)
local_prefixes = ("CC", "CTS", "RUA")


# =============================================================================
# FUNCTIONS
# =============================================================================

def parse_taxonomy_from_sseqid(sseqid):
    """
    Parses taxonomy from APSCALE sseqid string.
    Format: accession###root_1;Eukaryota_2759;Chordata_7711;...
    Returns dict with keys: Kingdom, Phylum, Class, Order, Family, Genus, Species
    """
    if "###" not in sseqid:
        return None

    tax_string = sseqid.split("###")[1]
    nodes = tax_string.split(";")

    # nodes[0] = root_1, nodes[1:] = taxonomy ranks
    if len(nodes) < 8:
        return None

    def clean_node(node):
        # Remove trailing taxon ID (_numbers) and underscores
        name = re.sub(r"_\d+$", "", node)
        return name.replace("_", " ").strip()

    return {
        "Kingdom": clean_node(nodes[1]),
        "Phylum":  clean_node(nodes[2]),
        "Class":   clean_node(nodes[3]),
        "Order":   clean_node(nodes[4]),
        "Family":  clean_node(nodes[5]),
        "Genus":   clean_node(nodes[6]),
        "Species": clean_node(nodes[7]),
    }


def load_all_blast_hits(blast_folder):
    """
    Loads all subset_X_blastn.csv files and returns a dict:
        { qseqid: [ (sseqid, pident, evalue), ... ] sorted by pident desc }
    """
    csv_files = sorted(glob.glob(str(blast_folder / "subset_*_blastn.csv")))
    if not csv_files:
        raise FileNotFoundError(
            f"No subset_*_blastn.csv files found in:\n  {blast_folder}\n"
            "Check your blast_folder path in CONFIG."
        )

    print(f"Loading {len(csv_files)} subset BLAST files...")
    all_hits = {}

    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file, header=None, sep=";;", engine="python")
            df.columns = ["qseqid", "sseqid", "pident", "evalue"]
            for _, row in df.iterrows():
                qid = row["qseqid"]
                if qid not in all_hits:
                    all_hits[qid] = []
                all_hits[qid].append((row["sseqid"], float(row["pident"]), float(row["evalue"])))
        except Exception as e:
            print(f"  Warning: could not read {csv_file}: {e}")

    # Sort each OTU's hits by similarity descending
    for qid in all_hits:
        all_hits[qid].sort(key=lambda x: x[1], reverse=True)

    print(f"Loaded hits for {len(all_hits)} unique OTUs.")
    return all_hits


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("Loading taxonomy file...")
    tax = pd.read_excel(taxonomy_file)
    print(f"  {len(tax)} OTUs loaded.")

    # Load all BLAST hits
    all_hits = load_all_blast_hits(blast_folder)

    # Find NoMatch OTUs with real similarity (not 0.0)
    nomatch_mask = (tax["Kingdom"].isin(["NoMatch", "Missing-Accession"])) & (tax["Similarity"] > 0)
    nomatch_otus = tax[nomatch_mask]
    print(f"\nNoMatch OTUs with similarity > 0: {len(nomatch_otus)}")

    overrides = 0
    skipped   = 0

    for idx, row in nomatch_otus.iterrows():
        otu_id = row["unique ID"]

        if otu_id not in all_hits:
            print(f"  [NOT FOUND] {otu_id[:20]}... — no BLAST hits in CSVs")
            skipped += 1
            continue

        hits = all_hits[otu_id]
        top_sseqid, top_pident, top_evalue = hits[0]

        # Check local_only filter
        if local_only:
            acc = top_sseqid.split("###")[0]
            if not any(acc.startswith(p) for p in local_prefixes):
                # Find best local hit instead
                local_hits = [(s, p, e) for s, p, e in hits
                              if any(s.split("###")[0].startswith(pr) for pr in local_prefixes)]
                if not local_hits:
                    print(f"  [SKIP] {otu_id[:20]}... — no local hits found")
                    skipped += 1
                    continue
                top_sseqid, top_pident, top_evalue = local_hits[0]

        # Check threshold
        if top_pident < override_threshold:
            print(f"  [BELOW THRESHOLD] {otu_id[:20]}... — top hit {top_pident:.1f}% < {override_threshold}%")
            skipped += 1
            continue

        # Parse taxonomy from top hit
        taxonomy = parse_taxonomy_from_sseqid(top_sseqid)
        if taxonomy is None:
            print(f"  [PARSE FAILED] {otu_id[:20]}... — could not parse taxonomy from sseqid")
            skipped += 1
            continue

        acc = top_sseqid.split("###")[0]
        is_local = any(acc.startswith(p) for p in local_prefixes)
        source = "local barcode" if is_local else "MIDORI"

        print(f"  [OVERRIDE] {otu_id[:20]}... {top_pident:.1f}% → {taxonomy['Species']} ({source})")

        # Apply override
        for col in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]:
            tax.at[idx, col] = taxonomy[col]

        tax.at[idx, "Similarity"] = top_pident
        tax.at[idx, "evalue"]     = top_evalue
        tax.at[idx, "Flag"]       = f"Manual override - top hit ({top_pident:.1f}%, {source})"
        tax.at[idx, "Status"]     = "overridden"

        overrides += 1

    # Summary
    print(f"\n--- Summary ---")
    print(f"  Overridden: {overrides}")
    print(f"  Skipped:    {skipped}")
    print(f"  Unchanged (similarity=0): {len(tax[tax['Kingdom'].isin(['NoMatch','Missing-Accession']) & (tax['Similarity']==0)])}")
    # Write output
    tax.to_excel(output_file, index=False)
    print(f"\nOutput written to:\n  {output_file}")


if __name__ == "__main__":
    main()
