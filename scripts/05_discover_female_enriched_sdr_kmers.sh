#!/bin/bash
set -euo pipefail

# ============================================================
# Female-enriched k-mer discovery
#
# Inputs:
#   female.samples
#   male.samples
#   kmers/SAMPLE.k21.txt
#
# Each k-mer file:
#   KMER [COUNT]
#
# Outputs:
#   all_kmer_prevalence.tsv
#   female_SDR_95pct_1pct.tsv
#   female_SDR_90pct_5pct.tsv
#   female_SDR_80pct_10pct.tsv
#   corresponding kmer-only files
#
# Designed for:
#   30 CPUs
#   ~300 GB RAM
# ============================================================


# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

FEMALE_LIST="female.samples"
MALE_LIST="male.samples"
KMER_DIR="kmers"

THREADS="${SLURM_CPUS_PER_TASK:-30}"

# Keep below your total 300 GB SLURM allocation
SORT_MEM="180G"

# IMPORTANT:
# Use current filesystem, NOT /tmp
TMP="$PWD/female_enriched_tmp_${SLURM_JOB_ID:-$$}"

mkdir -p "$TMP"


# ------------------------------------------------------------
# OUTPUTS
# ------------------------------------------------------------

ALL_PREVALENCE="all_kmer_prevalence.tsv"

OUT95="female_SDR_95pct_1pct.tsv"
OUT90="female_SDR_90pct_5pct.tsv"
OUT80="female_SDR_80pct_10pct.tsv"

KMER95="female_SDR_95pct_1pct.k21.txt"
KMER90="female_SDR_90pct_5pct.k21.txt"
KMER80="female_SDR_80pct_10pct.k21.txt"


# ------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------

die() {
    echo "ERROR: $*" >&2
    exit 1
}


# ------------------------------------------------------------
# CHECK INPUTS
# ------------------------------------------------------------

[ -f "$FEMALE_LIST" ] || die "$FEMALE_LIST not found"
[ -f "$MALE_LIST" ]   || die "$MALE_LIST not found"
[ -d "$KMER_DIR" ]    || die "$KMER_DIR directory not found"


NF=$(grep -cv '^[[:space:]]*$' "$FEMALE_LIST")
NM=$(grep -cv '^[[:space:]]*$' "$MALE_LIST")


echo "============================================================"
echo "Female-enriched k-mer discovery"
echo "============================================================"
echo "Female samples : $NF"
echo "Male samples   : $NM"
echo "Threads        : $THREADS"
echo "Sort memory    : $SORT_MEM"
echo "Temporary dir  : $TMP"
echo
echo "Disk space:"
df -h "$PWD"
echo "============================================================"
echo


# ------------------------------------------------------------
# CHECK ALL SAMPLE FILES
# ------------------------------------------------------------

echo "Checking female and male k-mer files..."

missing=0

while read -r sample
do
    sample=$(printf "%s" "$sample" | tr -d '\r')
    [ -z "$sample" ] && continue

    file="${KMER_DIR}/${sample}.k21.txt"

    if [ ! -f "$file" ]; then
        echo "MISSING: $file"
        missing=$((missing + 1))
    fi

done < "$FEMALE_LIST"


while read -r sample
do
    sample=$(printf "%s" "$sample" | tr -d '\r')
    [ -z "$sample" ] && continue

    file="${KMER_DIR}/${sample}.k21.txt"

    if [ ! -f "$file" ]; then
        echo "MISSING: $file"
        missing=$((missing + 1))
    fi

done < "$MALE_LIST"


if [ "$missing" -gt 0 ]; then
    die "$missing k-mer files are missing"
fi

echo "All sample files found."
echo


# ============================================================
# STEP 1
# FEMALE PREVALENCE
# ============================================================

echo "============================================================"
echo "STEP 1: Collecting female k-mers"
echo "============================================================"

FEMALE_ALL="$TMP/female_all.kmers"
FEMALE_COUNTS="$TMP/female_counts.tsv"

: > "$FEMALE_ALL"

COUNT=0

while read -r sample
do
    sample=$(printf "%s" "$sample" | tr -d '\r')
    [ -z "$sample" ] && continue

    COUNT=$((COUNT + 1))

    FILE="${KMER_DIR}/${sample}.k21.txt"

    echo "[$COUNT/$NF] Female: $sample"

    # Only first column = k-mer
    awk '{print $1}' "$FILE" >> "$FEMALE_ALL"

done < "$FEMALE_LIST"


echo
echo "Female concatenation complete."

echo -n "Female observations: "
wc -l < "$FEMALE_ALL"

echo
echo "Sorting female k-mers and counting prevalence..."

LC_ALL=C sort \
    --parallel="$THREADS" \
    -S "$SORT_MEM" \
    -T "$TMP" \
    "$FEMALE_ALL" \
| uniq -c \
| awk 'BEGIN{OFS="\t"} {print $2,$1}' \
> "$FEMALE_COUNTS"


echo -n "Unique female-observed k-mers: "
wc -l < "$FEMALE_COUNTS"

echo


# Remove huge concatenated file once no longer needed
rm -f "$FEMALE_ALL"


# ============================================================
# STEP 2
# MALE PREVALENCE
# ============================================================

echo "============================================================"
echo "STEP 2: Collecting male k-mers"
echo "============================================================"

MALE_ALL="$TMP/male_all.kmers"
MALE_COUNTS="$TMP/male_counts.tsv"

: > "$MALE_ALL"

COUNT=0

while read -r sample
do
    sample=$(printf "%s" "$sample" | tr -d '\r')
    [ -z "$sample" ] && continue

    COUNT=$((COUNT + 1))

    FILE="${KMER_DIR}/${sample}.k21.txt"

    echo "[$COUNT/$NM] Male: $sample"

    awk '{print $1}' "$FILE" >> "$MALE_ALL"

done < "$MALE_LIST"


echo
echo "Male concatenation complete."

echo -n "Male observations: "
wc -l < "$MALE_ALL"

echo
echo "Sorting male k-mers and counting prevalence..."

LC_ALL=C sort \
    --parallel="$THREADS" \
    -S "$SORT_MEM" \
    -T "$TMP" \
    "$MALE_ALL" \
| uniq -c \
| awk 'BEGIN{OFS="\t"} {print $2,$1}' \
> "$MALE_COUNTS"


echo -n "Unique male-observed k-mers: "
wc -l < "$MALE_COUNTS"

echo


rm -f "$MALE_ALL"


# ============================================================
# STEP 3
# SORT COUNT TABLES FOR JOIN
# ============================================================

echo "============================================================"
echo "STEP 3: Sorting count tables for join"
echo "============================================================"

FEMALE_SORTED="$TMP/female_counts.sorted.tsv"
MALE_SORTED="$TMP/male_counts.sorted.tsv"


LC_ALL=C sort \
    --parallel="$THREADS" \
    -S "$SORT_MEM" \
    -T "$TMP" \
    -k1,1 \
    "$FEMALE_COUNTS" \
> "$FEMALE_SORTED"


LC_ALL=C sort \
    --parallel="$THREADS" \
    -S "$SORT_MEM" \
    -T "$TMP" \
    -k1,1 \
    "$MALE_COUNTS" \
> "$MALE_SORTED"


# ============================================================
# STEP 4
# JOIN FEMALE + MALE COUNTS
# ============================================================

echo
echo "============================================================"
echo "STEP 4: Joining female and male k-mer prevalence"
echo "============================================================"

COMBINED="$TMP/combined.tsv"


join \
    -t $'\t' \
    -a 1 \
    -e 0 \
    -o '1.1,1.2,2.2' \
    "$FEMALE_SORTED" \
    "$MALE_SORTED" \
> "$COMBINED"


echo -n "Female k-mers evaluated: "
wc -l < "$COMBINED"


# ============================================================
# STEP 5
# CALCULATE PREVALENCE
# ============================================================

echo
echo "============================================================"
echo "STEP 5: Calculating prevalence and enrichment"
echo "============================================================"


awk \
-v NF="$NF" \
-v NM="$NM" '
BEGIN {
    OFS="\t"

    print \
    "kmer", \
    "female_count", \
    "male_count", \
    "female_prevalence", \
    "male_prevalence", \
    "female_percent", \
    "male_percent", \
    "prevalence_difference", \
    "enrichment_ratio"
}

{
    fp=$2/NF
    mp=$3/NM

    diff=fp-mp

    enrichment=(fp+0.001)/(mp+0.001)

    print \
    $1, \
    $2, \
    $3, \
    fp, \
    mp, \
    fp*100, \
    mp*100, \
    diff, \
    enrichment
}
' "$COMBINED" \
> "$ALL_PREVALENCE"


echo -n "Rows in prevalence table: "
tail -n +2 "$ALL_PREVALENCE" | wc -l


# ============================================================
# STEP 6
# SELECT CANDIDATE SETS
# ============================================================

echo
echo "============================================================"
echo "STEP 6: Selecting female-enriched candidate sets"
echo "============================================================"


# ------------------------------------------------------------
# Female >=95%, Male <=1%
# ------------------------------------------------------------

awk '
BEGIN {FS=OFS="\t"}

NR==1 {
    print
    next
}

$4 >= 0.95 && $5 <= 0.01 {
    print
}
' "$ALL_PREVALENCE" \
> "$OUT95"


awk '
BEGIN {FS=OFS="\t"}
NR>1 {print $1}
' "$OUT95" \
> "$KMER95"


# ------------------------------------------------------------
# Female >=90%, Male <=5%
# ------------------------------------------------------------

awk '
BEGIN {FS=OFS="\t"}

NR==1 {
    print
    next
}

$4 >= 0.90 && $5 <= 0.05 {
    print
}
' "$ALL_PREVALENCE" \
> "$OUT90"


awk '
BEGIN {FS=OFS="\t"}
NR>1 {print $1}
' "$OUT90" \
> "$KMER90"


# ------------------------------------------------------------
# Female >=80%, Male <=10%
# ------------------------------------------------------------

awk '
BEGIN {FS=OFS="\t"}

NR==1 {
    print
    next
}

$4 >= 0.80 && $5 <= 0.10 {
    print
}
' "$ALL_PREVALENCE" \
> "$OUT80"


awk '
BEGIN {FS=OFS="\t"}
NR>1 {print $1}
' "$OUT80" \
> "$KMER80"


# ============================================================
# STEP 7
# SUMMARY
# ============================================================

N95=$(tail -n +2 "$OUT95" | wc -l)
N90=$(tail -n +2 "$OUT90" | wc -l)
N80=$(tail -n +2 "$OUT80" | wc -l)


echo
echo "============================================================"
echo "ANALYSIS FINISHED"
echo "============================================================"
echo
echo "Female samples : $NF"
echo "Male samples   : $NM"
echo
echo "Candidate counts:"
echo
echo "Female >=95%, Male <=1%  : $N95"
echo "Female >=90%, Male <=5%  : $N90"
echo "Female >=80%, Male <=10% : $N80"
echo
echo "Main prevalence table:"
echo "  $ALL_PREVALENCE"
echo
echo "Candidate tables:"
echo "  $OUT95"
echo "  $OUT90"
echo "  $OUT80"
echo
echo "K-mer-only files:"
echo "  $KMER95"
echo "  $KMER90"
echo "  $KMER80"
echo
echo "Temporary working directory:"
echo "  $TMP"
echo
echo "Disk usage of temporary directory:"
du -sh "$TMP" || true
echo
echo "============================================================"
