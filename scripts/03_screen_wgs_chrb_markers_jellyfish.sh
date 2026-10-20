#!/usr/bin/env bash
set -euo pipefail

# Screen paired-end WGS libraries with a retained ChrB-specific 21-mer panel.
#
# Usage, one sample per Slurm array task:
#   sbatch --array=1-$(wc -l < samples.tsv) --cpus-per-task=30 --mem=50G \
#     03_screen_wgs_chrb_markers_jellyfish.sh samples.tsv WGS_marker_set_chrB_specific.k21.txt results
#
# samples.tsv must have three tab-separated columns and no header:
#   sample_id<TAB>/path/to/R1.fastq.gz<TAB>/path/to/R2.fastq.gz
#
# Output: one result file per sample: results/SAMPLE.tsv. Concatenate these
# after the array finishes. Independent outputs avoid concurrent appends.

SAMPLE_SHEET=${1:?"Provide samples.tsv"}
MARKERS=${2:?"Provide chrB-specific 21-mer list"}
OUTDIR=${3:?"Provide output directory"}

THREADS=${SLURM_CPUS_PER_TASK:-30}
HASH_SIZE=${HASH_SIZE:-500M}
K=${K:-21}
READ_PAIR_BP=${READ_PAIR_BP:-300}
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

mkdir -p "$OUTDIR"
if [[ ! -s "$SAMPLE_SHEET" || ! -s "$MARKERS" ]]; then
  echo "ERROR: sample sheet or marker list is missing/empty" >&2
  exit 1
fi

line=$(sed -n "${TASK_ID}p" "$SAMPLE_SHEET")
if [[ -z "$line" ]]; then
  echo "ERROR: no sample at array index ${TASK_ID}" >&2
  exit 1
fi

IFS=$'\t' read -r sample r1 r2 <<< "$line"
if [[ -z "${sample:-}" || ! -r "${r1:-}" || ! -r "${r2:-}" ]]; then
  echo "ERROR: expected sample_id, readable R1, and readable R2 in row ${TASK_ID}" >&2
  exit 1
fi

module load jellyfish
tmp=$(mktemp -d "${TMPDIR:-/tmp}/${sample}.chrB.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
jf="$tmp/${sample}.k${K}.jf"
query="$tmp/${sample}.marker_query.tsv"

# -C merges every k-mer with its reverse complement (canonical counting).
jellyfish count -C -m "$K" -s "$HASH_SIZE" -t "$THREADS" -o "$jf" \
  <(zcat -- "$r1") <(zcat -- "$r2")

# jellyfish query reports kmer<TAB>count for every marker, including zeroes.
jellyfish query -i "$jf" < "$MARKERS" > "$query"

marker_n=$(wc -l < "$MARKERS")
detected_n=$(awk '$2 > 0 {n++} END {print n + 0}' "$query")
total_marker_instance_count=$(awk '{s += $2} END {print s + 0}' "$query")
read_pairs=$(zcat -- "$r1" | awk 'END {print NR / 4}')
if [[ "$read_pairs" == "0" ]]; then
  echo "ERROR: zero reads detected in ${r1}" >&2
  exit 1
fi

# DBWGS normalization: distinct marker k-mers detected per million sequenced
# bases. This matches the original set-intersection implementation, where a
# marker is counted once when observed at least once in the library.
normalized=$(awk -v h="$detected_n" -v n="$read_pairs" -v bp="$READ_PAIR_BP" \
  'BEGIN {printf "%.10f", (h * 1000000) / (n * bp)}')
predicted_copy=$(awk -v a="$normalized" 'BEGIN {printf "%.10f", a / 12}')
proportion=$(awk -v d="$detected_n" -v m="$marker_n" 'BEGIN {printf "%.8f", d / m}')

{
  printf 'sample\tmarker_panel_size\tchrB_specific_kmer_hits\tdetected_marker_fraction\ttotal_marker_instance_count\tpaired_end_read_pairs\tnormalized_abundance_per_million_bases\tpredicted_chrB_copy_number\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$sample" "$marker_n" "$detected_n" "$proportion" "$total_marker_instance_count" \
    "$read_pairs" "$normalized" "$predicted_copy"
} > "$OUTDIR/${sample}.tsv"

echo "Finished ${sample}: ${detected_n}/${marker_n} markers detected"
