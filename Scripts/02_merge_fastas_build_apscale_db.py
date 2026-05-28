"""
0_2_merge_fastas_build_apscale_db.py

Builds an APSCALE4-compatible BLAST database from local Sanger reference
sequences and the MIDORI database.

HOW APSCALE4 LOOKUP WORKS (from b_filter.py source):
    1. BLAST runs and returns sseqid — the full sequence ID from makeblastdb
    2. For MIDORI entries the sseqid = full header e.g.:
           JX066741.1.<1.>592###root_1;Eukaryota_2759;Chordata_7711;...
    3. apscale_blast extracts 'accession' from sseqid by taking everything
       before the first '##' separator
    4. It then looks up that accession in db_taxonomy.parquet.snappy:
           taxid_df[taxid_df['Accession'] == accession]
    5. Falls back to startswith(accession[:30]) if exact match fails
    6. Returns 'Missing-Accession' if both fail

WHAT THIS MEANS FOR OUR DATABASE:
    - For MIDORI entries: accession = 'JX066741.1.<1.>592' (before ###)
    - For local entries:  accession = 'CC59.Poropuntius_normani' (before ###)
    - The parquet Accession column must contain these exact strings
    - The FASTA headers must produce these exact strings as the BLAST sseqid
    - makeblastdb uses everything before first whitespace as sseqid
    - Therefore: NO space in header, keep ### delimiter intact

INPUT FASTA FORMAT (from R script output):
    >CC59.Poropuntius_normani###root_1;Eukaryota_2759;Chordata_7711;...
    SEQUENCE

OUTPUT:
    - Combined FASTA (local + MIDORI)
    - Gzipped combined FASTA
    - BLAST database (via makeblastdb)
    - db_taxonomy.parquet.snappy (Accession -> taxonomy lookup)
    - Zipped database folder for APSCALE

REQUIRES:
    conda activate apscale4
    pip install biopython pandas pyarrow openpyxl

USAGE:
    Edit CONFIG section, then:
    python 0_2_merge_fastas_build_apscale_db.py

AUTHORS: Christian Ching — Cardamoms eDNA Metabarcoding Project
"""

import os
import re
import gzip
import shutil
import datetime
import pandas as pd
from pathlib import Path


# =============================================================================
# CONFIG
# =============================================================================

# Local reference FASTA from format_local_refs_apscale.R
# Headers: >CC59.Poropuntius_normani###root_1;Eukaryota_2759;...
LOCAL_FASTA = Path("/Users/christianching/Documents/Projects/cardamoms_eDNA/databases/sanger_databases/2025cardamoms16S_database_midorimode.fasta")

# MIDORI BLAST FASTA (plain or gzipped)
MIDORI_FASTA = Path("/Users/christianching/Documents/Projects/cardamoms_eDNA/databases/blast_databases/MIDORI2_UNIQ_NUC_GB266_lrRNA_BLAST.fasta")

# Output directory
OUTPUT_DIR = Path("/Users/christianching/Documents/Projects/cardamoms_eDNA/databases/blast_databases/")

# Database name
COMBINED_NAME = "16S_MIDORI_2025cardamoms_combined"

# Local sequence prefixes
LOCAL_PREFIXES = ("CC", "CTS", "RUA")


# =============================================================================
# FUNCTIONS
# =============================================================================

def log(msg):
    ts = datetime.datetime.now().strftime("%H:%M:%S")
    print(f"{ts} : {msg}")


def parse_taxonomy_from_header(header):
    """
    Parses accession and taxonomy from a MIDORI-format header.

    Input:  >CC59.Poropuntius_normani###root_1;Eukaryota_2759;Chordata_7711;
              Actinopteri_186623;Cypriniformes_7952;Cyprinidae_7953;
              Poropuntius_497005;Poropuntius_normani_643356

    Returns:
        accession: 'CC59.Poropuntius_normani'
        taxonomy:  ['Eukaryota', 'Chordata', 'Actinopteri', 'Cypriniformes',
                    'Cyprinidae', 'Poropuntius', 'Poropuntius normani']

    The accession is everything before ### — this is what makeblastdb stores
    as the sequence ID (no whitespace in header so full pre-### string is used)
    and what APSCALE extracts for parquet lookup.
    """
    header = header.lstrip(">").strip()

    if "###" not in header:
        return None, None

    # CRITICAL: APSCALE b_filter.py uses the FULL sseqid string as the
    # parquet Accession lookup key — no splitting on ### before lookup.
    # The full header (minus >) must be stored as Accession.
    accession = header  # full string: e.g. CC59.Poropuntius_normani###root_1;Eukaryota_2759;...
    tax_string = header.split("###")[1]
    nodes = tax_string.split(";")

    # nodes[0] = root_1, nodes[1:] = 7 taxonomy ranks
    def clean(node):
        name = re.sub(r"_\d+$", "", node)
        return name.replace("_", " ").strip()

    taxonomy = [clean(n) for n in nodes[1:]]

    # Pad to 7 if short, trim if long
    while len(taxonomy) < 7:
        taxonomy.append("")
    taxonomy = taxonomy[:7]

    return accession, taxonomy


def combine_and_parse(local_fasta, midori_fasta, output_fasta):
    """
    Combines local and MIDORI FASTAs into one file.
    Simultaneously parses taxonomy into a list for the parquet.
    Returns taxonomy records list.
    """
    log("Combining FASTAs and parsing taxonomy...")
    records = []
    n_local = 0
    n_midori = 0
    n_failed = 0

    with open(output_fasta, "w") as out:

        # --- Local sequences: write as-is, headers unchanged ---
        log(f"  Local sequences: {local_fasta.name}")
        with open(local_fasta) as f:
            current_header = None
            for line in f:
                out.write(line)
                if line.startswith(">"):
                    current_header = line.strip()
                    accession, taxonomy = parse_taxonomy_from_header(current_header)
                    if accession is not None:
                        records.append([accession] + taxonomy)
                        n_local += 1
                    else:
                        log(f"  WARNING: could not parse header: {current_header[:80]}")
                        n_failed += 1

        # --- MIDORI sequences ---
        log(f"  MIDORI sequences: {midori_fasta.name}")
        open_fn = gzip.open if str(midori_fasta).endswith(".gz") else open
        mode = "rt" if str(midori_fasta).endswith(".gz") else "r"
        with open_fn(midori_fasta, mode) as f:
            for line in f:
                out.write(line)
                if line.startswith(">"):
                    accession, taxonomy = parse_taxonomy_from_header(line.strip())
                    if accession is not None:
                        records.append([accession] + taxonomy)
                        n_midori += 1

    log(f"  Local:  {n_local} sequences")
    log(f"  MIDORI: {n_midori} sequences")
    log(f"  Failed to parse: {n_failed}")
    log(f"  Total:  {n_local + n_midori}")

    return records


def validate_local_sequences(records):
    """
    Validates that local sequences parsed correctly.
    Checks accession format, node count, Eukaryota presence.
    """
    log("Validating local sequence taxonomy...")
    issues = []

    local_records = [r for r in records
                     if any(r[0].startswith(p) for p in LOCAL_PREFIXES)]

    for r in local_records:
        accession = r[0]
        taxonomy = r[1:]

        if taxonomy[0] != "Eukaryota":
            issues.append(f"  BAD superkingdom '{taxonomy[0]}': {accession}")
        if "" in taxonomy:
            issues.append(f"  EMPTY ranks {[i for i,t in enumerate(taxonomy) if t=='']}: {accession}")

    if issues:
        log(f"  WARNING: {len(issues)} local sequences have issues:")
        for i in issues[:10]:
            log(i)
    else:
        log(f"  OK: all {len(local_records)} local sequences validated")

    return len(issues) == 0


def build_parquet(records, db_folder):
    """
    Writes db_taxonomy.parquet.snappy — the lookup table APSCALE uses.

    Column: Accession = everything before ### in the FASTA header
              = what makeblastdb stores as sequence ID
              = what APSCALE extracts from sseqid for lookup

    Columns: superkingdom, phylum, class, order, family, genus, species
    """
    df = pd.DataFrame(
        records,
        columns=["Accession", "superkingdom", "phylum", "class",
                 "order", "family", "genus", "species"]
    )

    # Remove duplicates (same accession appearing in both local and MIDORI)
    n_before = len(df)
    df = df.drop_duplicates(subset="Accession", keep="first")
    n_after = len(df)
    if n_before != n_after:
        log(f"  Removed {n_before - n_after} duplicate accessions (local takes priority)")

    parquet_path = db_folder / "db_taxonomy.parquet.snappy"
    df.to_parquet(parquet_path, compression="snappy")
    log(f"  Parquet saved: {parquet_path}")
    log(f"  Total accessions: {len(df)}")

    # Spot check local entries in parquet
    local_in_parquet = df[df["Accession"].str.startswith(LOCAL_PREFIXES, na=False)]
    log(f"  Local sequences in parquet: {len(local_in_parquet)}")

    # Show a sample
    if len(local_in_parquet) > 0:
        sample = local_in_parquet.iloc[0]
        log(f"  Sample local entry:")
        log(f"    Accession:    '{sample['Accession']}'")
        log(f"    superkingdom: '{sample['superkingdom']}'")
        log(f"    species:      '{sample['species']}'")

    return df


def gzip_fasta(fasta_path):
    gz_path = fasta_path.with_suffix(".fasta.gz")
    log(f"Gzipping to {gz_path.name}...")
    with open(fasta_path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return gz_path


def build_blast_db(fasta_gz, db_path):
    log("Building BLAST database with makeblastdb...")
    os.makedirs(db_path.parent, exist_ok=True)
    command = f"zcat < {fasta_gz} | makeblastdb -in - -title db -dbtype nucl -out {db_path}"
    log(f"  {command}")
    ret = os.system(command)
    if ret != 0:
        raise RuntimeError(
            "makeblastdb failed. Ensure you are in the apscale4 conda environment "
            "and makeblastdb is on PATH."
        )
    log("  BLAST database built successfully.")


def zip_db_folder(output_dir, db_name):
    archive = output_dir / db_name
    shutil.make_archive(str(archive), "zip", str(output_dir / db_name))
    log(f"  Database zipped: {archive}.zip")


# =============================================================================
# MAIN
# =============================================================================

def main():
    log("=" * 60)
    log("APSCALE4 database build pipeline")
    log(f"  Local FASTA:  {LOCAL_FASTA}")
    log(f"  MIDORI FASTA: {MIDORI_FASTA}")
    log(f"  Output dir:   {OUTPUT_DIR}")
    log("=" * 60)

    if not LOCAL_FASTA.exists():
        raise FileNotFoundError(f"Local FASTA not found: {LOCAL_FASTA}")
    if not MIDORI_FASTA.exists():
        raise FileNotFoundError(f"MIDORI FASTA not found: {MIDORI_FASTA}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Combine FASTAs and parse taxonomy simultaneously
    combined_fasta = OUTPUT_DIR / f"{COMBINED_NAME}.fasta"
    records = combine_and_parse(LOCAL_FASTA, MIDORI_FASTA, combined_fasta)

    # Step 2: Validate local sequences
    validate_local_sequences(records)

    # Step 3: Gzip combined FASTA
    combined_gz = gzip_fasta(combined_fasta)

    # Step 4: Build BLAST database
    db_name = f"db_{COMBINED_NAME}"
    db_blast_path = OUTPUT_DIR / db_name / "db"
    build_blast_db(combined_gz, db_blast_path)

    # Step 5: Build taxonomy parquet
    log("Building taxonomy parquet...")
    build_parquet(records, OUTPUT_DIR / db_name)

    # Step 6: Zip
    log("Zipping database folder...")
    zip_db_folder(OUTPUT_DIR, db_name)

    log("=" * 60)
    log("Pipeline complete.")
    log(f"Database: {OUTPUT_DIR / db_name}.zip")
    log(f"Place the unzipped folder in your APSCALE_databases directory.")
    log("=" * 60)


if __name__ == "__main__":
    main()