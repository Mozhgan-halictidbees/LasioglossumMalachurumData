#!/bin/bash
set -euxo pipefail
set -o errtrace

# ==========================================
# Bee protein downloader + longest isoforms
# Final output: TAXID.faa
# Keeps only longest isoform per gene
# Handles nested unzip directories
# Produces log file
# Keeps add.sh at end
# ==========================================

WORKDIR="/media/....../new_phylostratigraphy/faa"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

LOGFILE="$WORKDIR/add.log"

echo "=======================================" | tee "$LOGFILE"
echo "STARTING RUN: $(date)" | tee -a "$LOGFILE"
echo "Working directory: $WORKDIR" | tee -a "$LOGFILE"
echo "=======================================" | tee -a "$LOGFILE"

# -------------------------------
# Download URLs + TaxIDs
# -------------------------------

declare -A URLS
declare -A TAXIDS

URLS[NMEL]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/NMEL.zip"
TAXIDS[NMEL]="2448451"

URLS[AVIR]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/AVIR.zip"
TAXIDS[AVIR]="115084"

URLS[Dnov]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/Dnov.zip"
TAXIDS[Dnov]="178035"

URLS[HQUA]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/HQUA.zip"
TAXIDS[HQUA]="115107"

URLS[LALB]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/LALB.zip"
TAXIDS[LALB]="88501"

URLS[LLEU]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/LLEU.zip"
TAXIDS[LLEU]="88532"

URLS[LPAU]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/LPAU.zip"
TAXIDS[LPAU]="88516"

URLS[LZEP]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/LZEP.zip"
TAXIDS[LZEP]="88500"

URLS[Mgen]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/Mgen.zip"
TAXIDS[Mgen]="115081"

URLS[LFIG]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/LFIG.zip"
TAXIDS[LFIG]="160208"

URLS[LMAL]="https://beenomes.princeton.edu/wp-content/uploads/2021/04/LMAL.zip"
TAXIDS[LMAL]="88512"

# -------------------------------
# Main loop
# -------------------------------

for PREFIX in "${!URLS[@]}"; do

    echo "" | tee -a "$LOGFILE"
    echo "=======================================" | tee -a "$LOGFILE"
    echo "Processing $PREFIX" | tee -a "$LOGFILE"
    echo "=======================================" | tee -a "$LOGFILE"

    # -------------------------------
    # Remove leftovers from failed runs
    # -------------------------------

    rm -rf "$PREFIX"
    rm -f "${PREFIX}.zip"
    rm -f "${PREFIX}.zip.1"

    # -------------------------------
    # Download
    # -------------------------------

    echo "Downloading ${PREFIX}.zip ..." | tee -a "$LOGFILE"

    wget -c \
        --tries=20 \
        --timeout=30 \
        --waitretry=5 \
        --retry-connrefused \
        "${URLS[$PREFIX]}" \
        >> "$LOGFILE" 2>&1

    if [[ $? -ne 0 ]]; then
        echo "ERROR: download failed for $PREFIX" | tee -a "$LOGFILE"
        continue
    fi

    # -------------------------------
    # Unzip
    # -------------------------------

    echo "Unzipping ${PREFIX}.zip ..." | tee -a "$LOGFILE"

    unzip -o "${PREFIX}.zip" >> "$LOGFILE" 2>&1

    # -------------------------------
    # Find peptide fasta
    # Handles:
    #   ./PREFIX_pep.fasta
    #   ./PREFIX/....pep.fasta
    # -------------------------------

    PEP=$(find . -type f \( -name "*_pep.fasta" -o -name "*pep*.fa*" \) | grep "$PREFIX" | head -n 1)

    if [[ -z "$PEP" ]]; then
        echo "ERROR: peptide fasta not found for $PREFIX" | tee -a "$LOGFILE"
        continue
    fi

    echo "Protein file found: $PEP" | tee -a "$LOGFILE"

    # -------------------------------
    # Keep longest isoforms only
    # -------------------------------

    python3 <<EOF >> "$LOGFILE" 2>&1
from Bio import SeqIO

input_fasta = r"$PEP"
output_fasta = r"${TAXIDS[$PREFIX]}.faa"

best_isoforms = {}

for record in SeqIO.parse(input_fasta, "fasta"):

    prot_id = record.id

    # Example:
    # LMAL_12294-RA -> LMAL_12294
    gene_id = prot_id.split("-")[0]

    if gene_id not in best_isoforms:
        best_isoforms[gene_id] = record
    else:
        if len(record.seq) > len(best_isoforms[gene_id].seq):
            best_isoforms[gene_id] = record

with open(output_fasta, "w") as out_handle:
    SeqIO.write(best_isoforms.values(), out_handle, "fasta")

print(f"Saved: {output_fasta}")
print(f"Genes kept: {len(best_isoforms)}")
EOF

    if [[ $? -ne 0 ]]; then
        echo "ERROR: python processing failed for $PREFIX" | tee -a "$LOGFILE"
        continue
    fi

    echo "Finished $PREFIX" | tee -a "$LOGFILE"

    # -------------------------------
    # Cleanup
    # Keep:
    #   *.faa
    #   add.sh
    #   add.log
    # -------------------------------

    find . -mindepth 1 \
        ! -name "*.faa" \
        ! -name "add.sh" \
        ! -name "add.log" \
        -exec rm -rf {} + 2>/dev/null

done

echo "" | tee -a "$LOGFILE"
echo "=======================================" | tee -a "$LOGFILE"
echo "ALL DONE" | tee -a "$LOGFILE"
echo "=======================================" | tee -a "$LOGFILE"

echo "" | tee -a "$LOGFILE"
echo "Final FAA files:" | tee -a "$LOGFILE"
ls -lh *.faa | tee -a "$LOGFILE"

echo "" | tee -a "$LOGFILE"
echo "Log saved to:" | tee -a "$LOGFILE"
echo "$LOGFILE" | tee -a "$LOGFILE"

echo "=======================================" | tee -a "$LOGFILE"
