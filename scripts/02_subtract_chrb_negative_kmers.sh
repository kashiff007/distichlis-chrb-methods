#!/bin/bash
set -euo pipefail

# ============================================================
# Remove any k-mer found in ANY chrB-negative sample
#
# Input:
#   WGS_marker_candidates_shared_Bpositive.k21.txt
#   chrB_no.samples
#   kmers/SAMPLE.k21.txt
#
# Output:
#   WGS_marker_set_chrB_specific.k21.txt
# ============================================================

COMMON_BPLUS="WGS_marker_candidates_shared_Bpositive.k21.txt"
BMINUS_LIST="chrB_no.samples"
KMER_DIR="kmers"
TMPDIR="chrB_filter_tmp"
OUTPUT="WGS_marker_set_chrB_specific.k21.txt"

mkdir -p "$TMPDIR"

echo "========================================"
echo "chrB-specific k-mer filtering"
echo "========================================"

# Check input files
if [ ! -f "$COMMON_BPLUS" ]; then
    echo "ERROR: $COMMON_BPLUS not found"
    exit 1
fi

if [ ! -f "$BMINUS_LIST" ]; then
    echo "ERROR: $BMINUS_LIST not found"
    exit 1
fi


# ------------------------------------------------------------
# Sort starting B+ common k-mer set
# ------------------------------------------------------------

echo
echo "Preparing common B+ k-mers..."

LC_ALL=C sort -u "$COMMON_BPLUS" \
    > "$TMPDIR/current.kmers"

START=$(wc -l < "$TMPDIR/current.kmers")

echo "Starting common B+ k-mers: $START"


# ------------------------------------------------------------
# Count B- samples
# ------------------------------------------------------------

TOTAL=$(grep -cv '^[[:space:]]*$' "$BMINUS_LIST")

echo "B- samples to process: $TOTAL"
echo


# ------------------------------------------------------------
# Progressively remove k-mers occurring in each B- sample
# ------------------------------------------------------------

COUNT=0

while read -r sample
do
    # Remove Windows carriage return if present
    sample=$(echo "$sample" | tr -d '\r')

    [ -z "$sample" ] && continue

    COUNT=$((COUNT + 1))

    FILE="${KMER_DIR}/${sample}.k21.txt"

    echo "[$COUNT/$TOTAL] Processing B- sample: $sample"

    if [ ! -f "$FILE" ]; then
        echo "WARNING: Missing file: $FILE"
        continue
    fi


    # Extract k-mer column and sort
    awk '{print $1}' "$FILE" |
        LC_ALL=C sort -u \
        > "$TMPDIR/Bminus.current.sorted"


    # Keep only B+ k-mers NOT present in this B- sample
    LC_ALL=C comm -23 \
        "$TMPDIR/current.kmers" \
        "$TMPDIR/Bminus.current.sorted" \
        > "$TMPDIR/next.kmers"


    mv "$TMPDIR/next.kmers" \
       "$TMPDIR/current.kmers"


    REMAIN=$(wc -l < "$TMPDIR/current.kmers")

    echo "    Remaining chrB-specific candidates: $REMAIN"


    # Stop early if nothing remains
    if [ "$REMAIN" -eq 0 ]; then
        echo
        echo "WARNING: No candidate k-mers remain."
        break
    fi

done < "$BMINUS_LIST"


# ------------------------------------------------------------
# Final output
# ------------------------------------------------------------

cp "$TMPDIR/current.kmers" "$OUTPUT"

FINAL=$(wc -l < "$OUTPUT")
REMOVED=$((START - FINAL))


echo
echo "========================================"
echo "FINISHED"
echo "========================================"
echo "Starting common B+ k-mers : $START"
echo "Removed by B- samples     : $REMOVED"
echo "Final chrB-specific kmers : $FINAL"
echo
echo "Output:"
echo "$OUTPUT"
echo "========================================"


# Clean temporary per-sample file
rm -f "$TMPDIR/Bminus.current.sorted"
