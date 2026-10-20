#!/usr/bin/env bash
set -euo pipefail

# Build the assembly-derived (reference) ChrB 21-mer marker set described in
# the manuscript. The script retains canonical 21-mers shared by every listed
# ChrB assembly, then removes any k-mer detected in a standard-chromosome or
# organellar assembly.
#
# Usage:
#   bash 04_build_reference_chrb_marker_set.sh chrB_assemblies.list \
#        nonchrB_background_assemblies.list reference_marker_output
#
# Each list contains one readable FASTA path per line. The first list must
# include the three ChrB assemblies. The background list should contain the
# A-chromosome, chloroplast, and mitochondrial sequences used for subtraction.

B_LIST=${1:?"Provide chrB_assemblies.list"}
BACKGROUND_LIST=${2:?"Provide nonchrB_background_assemblies.list"}
OUTDIR=${3:?"Provide output directory"}

K=${K:-21}
THREADS=${SLURM_CPUS_PER_TASK:-30}
HASH_SIZE=${HASH_SIZE:-1G}

module load jellyfish
mkdir -p "$OUTDIR/tmp"

mapfile -t b_fastas < <(grep -vE '^[[:space:]]*(#|$)' "$B_LIST")
mapfile -t background_fastas < <(grep -vE '^[[:space:]]*(#|$)' "$BACKGROUND_LIST")
[[ ${#b_fastas[@]} -ge 3 ]] || { echo "ERROR: provide at least three ChrB assemblies" >&2; exit 1; }
[[ ${#background_fastas[@]} -ge 1 ]] || { echo "ERROR: provide at least one non-ChrB background assembly" >&2; exit 1; }

for fa in "${b_fastas[@]}" "${background_fastas[@]}"; do
  [[ -r "$fa" ]] || { echo "ERROR: unreadable FASTA: $fa" >&2; exit 1; }
done

first=1
for fa in "${b_fastas[@]}"; do
  label=$(basename "$fa" | sed 's/\.[^.]*$//')
  jf="$OUTDIR/tmp/${label}.k${K}.jf"
  kmers="$OUTDIR/tmp/${label}.canonical.k${K}.txt"
  jellyfish count -C -m "$K" -s "$HASH_SIZE" -t "$THREADS" -o "$jf" "$fa"
  jellyfish dump -c "$jf" | awk '{print $1}' | LC_ALL=C sort -u > "$kmers"
  if [[ "$first" == 1 ]]; then
    cp "$kmers" "$OUTDIR/tmp/shared.current.txt"
    first=0
  else
    LC_ALL=C comm -12 "$OUTDIR/tmp/shared.current.txt" "$kmers" > "$OUTDIR/tmp/shared.next.txt"
    mv "$OUTDIR/tmp/shared.next.txt" "$OUTDIR/tmp/shared.current.txt"
  fi
done

background_jf="$OUTDIR/tmp/nonchrB_background.k${K}.jf"
jellyfish count -C -m "$K" -s "$HASH_SIZE" -t "$THREADS" -o "$background_jf" "${background_fastas[@]}"
jellyfish query -i "$background_jf" < "$OUTDIR/tmp/shared.current.txt" \
  | awk '$2 == 0 {print $1}' \
  | LC_ALL=C sort -u > "$OUTDIR/reference_marker_set_chrB_specific.k${K}.txt"

printf 'Shared canonical ChrB k-mers before background subtraction: %s\n' "$(wc -l < "$OUTDIR/tmp/shared.current.txt")"
printf 'Reference-marker set after A-chromosome and organellar subtraction: %s\n' "$(wc -l < "$OUTDIR/reference_marker_set_chrB_specific.k${K}.txt")"
