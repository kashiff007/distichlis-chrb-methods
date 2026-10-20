#!/usr/bin/env bash
set -euo pipefail

# Run on a Linux compute node with PLINK 2 and GEMMA available on PATH.
# Example: module load plink2 gemma; bash run_missing10_chrB_GLM_GEMMA_LMM.sh

THREADS="${THREADS:-30}"
PLINK2="${PLINK2:-plink2}"
GEMMA="${GEMMA:-gemma}"
OUT="outputs/spicata_stricta_missing10_PCA"
PGEN="$OUT/spicata_stricta_missing10_QC"
PRUNE="$OUT/spicata_stricta_missing10_LD.prune.in"

mkdir -p "$OUT/gemma_output"
Rscript scripts_prepare_missing10_GEMMA_inputs.R

# Export the two tested marker sets to GEMMA-compatible PLINK 1 binary files.
"$PLINK2" --pfile "$PGEN" --make-bed --threads "$THREADS" --out "$OUT/spicata_stricta_missing10_allQC"
"$PLINK2" --pfile "$PGEN" --extract "$PRUNE" --make-bed --threads "$THREADS" --out "$OUT/spicata_stricta_missing10_LDpruned"

# PC- and coverage-adjusted logistic GLM for the binary ChrB phenotype.
"$PLINK2" --pfile "$PGEN" \
  --pheno "$OUT/chrB_status_unpruned.tsv" --pheno-name ChrB_status \
  --covar "$OUT/covariates_PC1_10_coverage_unpruned.tsv" \
  --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Coverage \
  --covar-variance-standardize --glm firth-fallback hide-covar --threads "$THREADS" \
  --out "$OUT/chrB_status_unpruned_GLM"
"$PLINK2" --pfile "$PGEN" --extract "$PRUNE" \
  --pheno "$OUT/chrB_status_LDpruned.tsv" --pheno-name ChrB_status \
  --covar "$OUT/covariates_PC1_10_coverage_LDpruned.tsv" \
  --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,Coverage \
  --covar-variance-standardize --glm firth-fallback hide-covar --threads "$THREADS" \
  --out "$OUT/chrB_status_LDpruned_GLM"

# Estimate one relationship matrix from LD-pruned markers.  Reusing it for both
# scans avoids including the full tested SNP set in the kinship estimate.
"$GEMMA" -bfile "$OUT/spicata_stricta_missing10_LDpruned" -gk 1 \
  -o chrB_LDpruned_kinship -outdir "$OUT/gemma_output"
KINSHIP="$OUT/gemma_output/chrB_LDpruned_kinship.cXX.txt"

# GEMMA LMM: 0/1 ChrB status, PC1-PC10, and standardized coverage as covariates.
# This is the conventional plant-GWAS LMM sensitivity analysis. For final binary
# effect sizes, verify lead loci with a logistic mixed model (e.g. SAIGE/GMMAT).
"$GEMMA" -bfile "$OUT/spicata_stricta_missing10_allQC" -k "$KINSHIP" \
  -p "$OUT/chrB_status_unpruned.gemma.pheno" -c "$OUT/chrB_status_unpruned.gemma.covar" \
  -lmm 4 -n 1 -o chrB_status_unpruned_GEMMA_LMM -outdir "$OUT/gemma_output"
"$GEMMA" -bfile "$OUT/spicata_stricta_missing10_LDpruned" -k "$KINSHIP" \
  -p "$OUT/chrB_status_LDpruned.gemma.pheno" -c "$OUT/chrB_status_LDpruned.gemma.covar" \
  -lmm 4 -n 1 -o chrB_status_LDpruned_GEMMA_LMM -outdir "$OUT/gemma_output"

echo "Finished. GEMMA results: $OUT/gemma_output/*.assoc.txt"

