# Manuscript-methods and script audit

This audit compares the repository workflow against the finalized manuscript. The manuscript was read only; it was not edited.

| Manuscript method | Publication-aligned repository script | Status | Notes |
| --- | --- | --- | --- |
| Canonical 21-mer counting with Jellyfish v2.3.0 | `01_intersect_chrb_positive_kmers.sh` input requirement; `03_screen_wgs_chrb_markers_jellyfish.sh` | Aligned | Inputs must be generated with Jellyfish `-C -m 21`; the screening script performs canonical counting directly. |
| Exact intersection across 26 cytology-confirmed ChrB-positive genets | `01_intersect_chrb_positive_kmers.sh` | Aligned | Produces `WGS_marker_candidates_shared_Bpositive.k21.txt`. |
| Remove all k-mers observed among 75 ChrB-negative genets | `02_subtract_chrb_negative_kmers.sh` | Aligned | Produces `WGS_marker_set_chrB_specific.k21.txt`. |
| Screen 757 paired-end WGS libraries and normalize by sequenced bases | `03_screen_wgs_chrb_markers_jellyfish.sh` | Aligned | Writes per-library output. Normalization uses distinct detected marker k-mers per million sequenced bases, consistent with the original marker-set intersection metric. |
| Reference-marker panel shared by three ChrB assemblies and absent from A/organellar sequences | `04_build_reference_chrb_marker_set.sh` | Aligned | Requires all three ChrB FASTAs and a background FASTA list. Uses canonical 21-mers and requires presence in every ChrB assembly. |
| Female-associated marker panel: at least 95% of females and at most 1% of males | `05_discover_female_enriched_sdr_kmers.sh` | Aligned | The 95%/1% output is the manuscript-aligned SDR candidate set. |
| SeqTK downsampling to 15x, 10x, 5x, 3x, 2x, 1x, 0.5x, 0.05x, 0.01x, 0.001x | `06_downsample_paired_end_wgs.sh` | Aligned | Uses the same seed for R1 and R2 to preserve read pairing. |
| Eastern *D. spicata* GWAS: GMMAT binomial logistic mixed model, five PCs, genomic relationship matrix, no LD pruning | Not yet provided | Pending | The final GMMAT preparation, model-fitting, annotation, and enrichment scripts are still required for repository completeness. |

## Archived scripts

The following supplied scripts are retained under `scripts/archive/` for provenance but are not linked as primary manuscript workflows:

- `reference_marker_prototype.sh`: did not require a k-mer to be present in each of the three ChrB assemblies and did not use canonical Jellyfish counting; therefore it does not implement the final manuscript reference-marker method exactly.
- `strict_female_kmer_intersection_exploratory.sh` and `remove_male_observed_female_kmers_exploratory.sh`: define an all-female/shared and all-male-absent marker set, whereas the manuscript uses a prevalence threshold of at least 95% in females and at most 1% in males.
- `GLM_GEMMA_sensitivity_analysis.sh`: implements PLINK GLM and GEMMA LMM sensitivity analyses for a different sample configuration. It is not the manuscript's final eastern *D. spicata* GMMAT logistic mixed-model pipeline.
