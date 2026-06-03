#!/bin/bash

# =========================================================
# eggNOG-mapper v2.1.12 (Lasioglossum malachurum)
# =========================================================

set -euxo pipefail
set -o errtrace
# ----------------------------
# USER INPUT
# ----------------------------

FASTA_FILE="/path/to/LMAR_OGS_v2.1.1_pep.fasta"
WORKDIR="/path/to/eggnog_project"

THREADS=4
PREFIX="Lmalachurum_eggnog"

# ----------------------------
# AUTO PATHS
# ----------------------------

DB_DIR="$WORKDIR/eggnog_db"
OUTPUT_DIR="$WORKDIR/output"

mkdir -p "$WORKDIR" "$DB_DIR" "$OUTPUT_DIR"

# ----------------------------
# STEP 1: DOWNLOAD eggNOG DATABASE (Hymenoptera)
# Tax ID: 7399 Taxonomy id for Hymenoptra
# ----------------------------

echo "Downloading eggNOG database (Hymenoptera)..."

python download_eggnog_data.py \
    -H \
    -d 7399 \
    --dbname "Hymenoptera" \
    -y \
    --data_dir "$DB_DIR"

# ----------------------------
# STEP 2: SET DIAMOND DATABASE PATH
# ----------------------------

DIAMOND_DB="$DB_DIR/eggnog_proteins.dmnd"

# ----------------------------
# STEP 3: RUN eggNOG-mapper
# ----------------------------

echo "Running eggNOG-mapper annotation..."

emapper.py \
    -i "$FASTA_FILE" \
    -m diamond \
    --dmnd_db "$DIAMOND_DB" \
    --output_dir "$OUTPUT_DIR" \
    -o "$PREFIX" \
    --cpu "$THREADS" \
    --evalue 0.001

# ----------------------------
# DONE
# ----------------------------

echo "Annotation completed successfully!"
echo "Results saved in: $OUTPUT_DIR"
