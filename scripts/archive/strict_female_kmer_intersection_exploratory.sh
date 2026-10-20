mkdir -p common_kmer_tmp

# Number of samples
N=$(grep -cv '^[[:space:]]*$' female.samples)
echo "Number of samples: $N"

# Sort unique k-mers from each sample
while read -r sample
do
    [ -z "$sample" ] && continue

    echo "Processing $sample"

    LC_ALL=C cut -d' ' -f1 "kmers/${sample}.k21.txt" \
        | LC_ALL=C sort -u \
        > "common_kmer_tmp/${sample}.kmers.sorted"

done < female.samples


first=$(grep -v '^[[:space:]]*$' female.samples | head -1)

cp "common_kmer_tmp/${first}.kmers.sorted" \
   common_kmer_tmp/common.current

grep -v '^[[:space:]]*$' female.samples | tail -n +2 | \
while read -r sample
do
    echo "Intersecting with $sample"

    LC_ALL=C comm -12 \
        common_kmer_tmp/common.current \
        "common_kmer_tmp/${sample}.kmers.sorted" \
        > common_kmer_tmp/common.next

    mv common_kmer_tmp/common.next \
       common_kmer_tmp/common.current

    echo -n "Remaining common k-mers: "
    wc -l < common_kmer_tmp/common.current

done



cp common_kmer_tmp/common.current female_SDR.all_common.k21.txt




