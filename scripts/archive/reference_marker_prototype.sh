#!/bin/bash
# ============================================================================
# CHRB-SPECIFIC K-MER DETECTION PIPELINE
# 
# This pipeline identifies chromosome B-specific k-mers and searches for their
# presence in whole genome sequencing data from multiple samples.
# 
# Pipeline overview:
# 1. Create a reference genome excluding chrB from 10 assemblies
# 2. Extract chrB from 3 assemblies and identify unique chrB k-mers
# 3. Filter k-mers to those unique to chrB (not found elsewhere)
# 4. Search for these unique k-mers in WGS data from thousands of samples
# ============================================================================

###############################################################################
# SECTION 1: CREATE COMPOSITE REFERENCE GENOME (chrB-EXCLUDED)
###############################################################################

# Objective: Concatenate a comprehensive reference genome containing all 
# chromosomes except chrB from 10 different assemblies. This creates a 
# "chrB-depleted" background genome for identifying chrB-specific sequences.

# Concatenate all genomic sequences from 10 assemblies, specifically including:
# - chrA chromosomes (autosomes)
# - chrP chromosomes (if present, potentially plastid/genome)
# - chrM chromosomes (mitochondrial genomes)
# - All other chromosomes EXCLUDING chrB
# This creates 'all.fasta' which serves as the background genome

cat all your genomes without chrB > all.fasta  

# Convert the multi-line FASTA to single-line format for easier processing
# -w 0: Set line width to 0 (no wrapping, single line per sequence)
# This creates 'all_linear.fasta' where each chromosome is on a single line

seqkit seq -w 0 all.fasta > all_linear.fasta

###############################################################################
# SECTION 2: EXTRACT AND CHARACTERIZE CHRB-SPECIFIC K-MERS
###############################################################################

# Objective: Extract chromosome B from 3 specific assemblies and identify 
# 21-mers (k=21) that are unique to or characteristic of chrB.

# Concatenate chrB sequences from three representative assemblies:
# 1. W6-54757-006_chrB.fa
# 2. 2021-KS-Q6_chrB.fa
# 3. W6-50019-003_chrB.fa
# This creates a comprehensive chrB reference set

cat W6-54757-006_chrB.fa 2021-KS-Q6_chrB W6-50019-003_chrB > chrB.fa

# Generate k-mer counts from the combined chrB sequences:
# -m 21: Use k-mer length of 21 bases
# -s 1G: Allocate 1GB hash table size for storing k-mers
# -t 30: Use 30 threads for parallel processing
# -o chrB_specific_kmers.jf: Output Jellyfish database file
# This creates a database of all 21-mers present in chrB sequences

jellyfish count -m 21 -s 1G -t 30 -o chrB_specific_kmers.jf chrB.fa

# Convert the binary Jellyfish database to human-readable text format
# -o chrB_specific_kmers.txt: Output text file with k-mer sequences
# Each line contains one 21-mer found in the chrB sequences

jellyfish dump -o chrB_specific_kmers.txt chrB_specific_kmers.jf

###############################################################################
# SECTION 3: FILTER FOR CHRB-UNIQUE K-MERS
###############################################################################

# Objective: Filter the chrB k-mers to only those that are NOT found in the 
# background genome (all chromosomes except chrB). These are truly chrB-specific.

# AWK script that performs the filtering:
# 1. First block (NR==FNR): Load chrB k-mers into array 'g'
# 2. Second block: Process the linearized background genome:
#    - Convert to uppercase
#    - Extract all 21-mers from the background genome
#    - Store them in array 'seen'
# 3. END block: Output chrB k-mers that are NOT in 'seen' array
#    These are k-mers unique to chrB (not found elsewhere in the genome)

awk 'NR==FNR {g[$0]=1; next} 
{
    seq=toupper($0)
    for(i=1;i<=length(seq)-20;i++){
        k=substr(seq,i,21)
        seen[k]=1
    }
}
END{
    for(k in g)
        if(!(k in seen))
            print k
}' chrB_specific_kmers.txt <(grep -v "^>" all_linear.fasta | tr -d '\n') > chrB_unique_only.txt

# Convert the unique k-mer list to FASTA format for use with seqkit
# Each k-mer gets a header like ">1", ">2", etc.

awk '{print ">" NR "\n" $0}'  genome_temp/chrB_unique_only.txt > genome_temp/chrB_unique_only.fasta

###############################################################################
# SECTION 4: PREPARE WGS DATA FILE LIST
###############################################################################

# Objective: Create a comprehensive list of all whole genome sequencing files
# to be analyzed, then split into manageable chunks for parallel processing.

# List all FASTQ files (both .fastq.gz and .fq.gz extensions) from the
# saltgrass WGS project directory. This includes paired-end and single-end data.

ls /ibex/project/c2141/saltgrass/wgs/*/*.fastq.gz /ibex/project/c2141/saltgrass/wgs/*/*.fq.gz > fastq_list.txt

# Split the large file list (~3500+ entries) into 3 manageable chunks for
# parallel processing using SLURM job arrays:
# - fastq_list1.txt: Lines 1-1500
# - fastq_list2.txt: Lines 1501-3000
# - fastq_list3.txt: Lines 3001-4000

sed -n '1,1500p' fastq_list.txt > fastq_list1.txt
sed -n '1501,3000p' fastq_list.txt > fastq_list2.txt
sed -n '3001,4000p' fastq_list.txt > fastq_list3.txt

###############################################################################
# SECTION 5: PARALLEL SEARCH FOR CHRB-SPECIFIC K-MERS IN WGS DATA
###############################################################################

# Objective: Search each WGS file for the presence of chrB-unique k-mers
# using parallel processing on a SLURM cluster.

# Calculate the number of files in the first list for array job sizing

N=$(wc -l < fastq_list1.txt)

# Submit a SLURM job array to process each FASTQ file in parallel:
# -t 1-01:00:00: Time limit of 1 day, 1 hour
# --job-name=count: Job name for tracking
# --output=slurm_out_%A_%a.txt: Output log files with job and array IDs
# --mem=50GB: Memory allocation per job
# -N 1: One node per job
# -c 30: 30 CPUs per job
# --array=1-$N: Array job with indices from 1 to N (number of files)

# Each array task does:
# 1. Extract filename for this array index
# 2. Use seqkit locate to search for chrB-unique k-mers in the file
# 3. Count matches (excluding header line)
# 4. Output filename and count to a TSV file

sbatch -t 1-01:00:00 --job-name=count --output=slurm_out_%A_%a.txt --mem=50GB -N 1 -c 30   --array=1-$N   --wrap 'f=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_list1.txt); c=$(seqkit locate -f chrB_unique_only.txt -j 30 "$f" | awk "NR>1{n++} END{print n+0}"); echo -e "$f\t$c" > chrB_count_${SLURM_ARRAY_TASK_ID}.tsv'

###############################################################################
# SECTION 6: PROCESS REMAINING FILE LISTS
###############################################################################

# Note: Repeat the SLURM job submission for the other two file lists
# (fastq_list2.txt and fastq_list3.txt) with appropriate adjustments:
# 1. Update the input file list in the sed command
# 2. Adjust the array job range based on file count in each list
# 3. Consider changing output file naming to avoid overwrites

# Example for processing the second list:
# N2=$(wc -l < fastq_list2.txt)
# sbatch ... --array=1-$N2 --wrap '
#   f=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_list2.txt)
#   ... (same processing logic) ...
# '

# Example for processing the third list:
# N3=$(wc -l < fastq_list3.txt)
# sbatch ... --array=1-$N3 --wrap '
#   f=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_list3.txt)
#   ... (same processing logic) ...
# '

###############################################################################
# POST-PROCESSING NOTES
###############################################################################

# After all jobs complete, combine results:
# cat chrB_count_*.tsv > all_chrB_counts.tsv

# Analyze results:
# 1. Files with high counts may contain chrB sequences
# 2. Files with zero counts likely lack chrB
# 3. Consider normalizing by read count for fair comparison

# Potential improvements:
# 1. Add quality control checks on input files
# 2. Implement result validation with positive/negative controls
# 3. Consider using k-mer counting tools (jellyfish) for more efficient searching
# 4. Add error handling and logging within each job
