#!/usr/bin/env bash
set -euo pipefail

# Paired-end WGS downsampling for one Slurm array task.
#
# Usage:
#   sbatch --array=1-$(wc -l < downsample_samples.tsv) --cpus-per-task=1 \
#     06_downsample_paired_end_wgs.sh downsample_samples.tsv downsampled
#
# Input columns: sample_id<TAB>R1.fastq.gz<TAB>R2.fastq.gz
# The identical seed is used for R1 and R2, preserving pair membership.

SAMPLES=${1:?"Provide downsample_samples.tsv"}
OUTDIR=${2:?"Provide output directory"}
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
SEED=${SEED:-100}

module load seqtk
line=$(sed -n "${TASK_ID}p" "$SAMPLES")
IFS=$'\t' read -r sample r1 r2 <<< "$line"
[[ -n "${sample:-}" && -r "${r1:-}" && -r "${r2:-}" ]] || { echo "ERROR: malformed or unreadable sample row" >&2; exit 1; }

# Fractions correspond to approximately 15x, 10x, 5x, 3x, 2x, 1x, 0.5x,
# 0.05x, 0.01x, and 0.001x from the original approximately 30x libraries.
fractions=(0.50 0.333333 0.166667 0.10 0.066667 0.033333 0.016667 0.001667 0.000333 0.0000333)
labels=(15x 10x 5x 3x 2x 1x 0.5x 0.05x 0.01x 0.001x)

for i in "${!fractions[@]}"; do
  fraction=${fractions[$i]}
  label=${labels[$i]}
  destination="$OUTDIR/$sample/$label"
  mkdir -p "$destination"
  seqtk sample -s "$SEED" "$r1" "$fraction" | gzip > "$destination/${sample}_R1.fastq.gz"
  seqtk sample -s "$SEED" "$r2" "$fraction" | gzip > "$destination/${sample}_R2.fastq.gz"
done
