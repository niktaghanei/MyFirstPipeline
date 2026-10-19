mkdir -p aligned log qc aligned


ls *.fq.gz > file_names.txt
cut -f6,7 Tco2.fqfiles | tr '\t' '\n' > metadata_names.txt
comm file_names.txt metadata_names.txt


fastqc -t 64 -o qc/ *.fq.gz
ls qc/*.zip | wc -l
for z in qc/*.zip; do unzip -q "$z" -d qc/unzipped/; done

for file in qc/unzipped/*/summary.txt; do
    grep -E "FAIL|WARN" "$file" >> error.txt
done

bowtie2-build --threads 32 TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta Tcongo_index
tail -n +2 Tco2.fqfiles | while IFS=$'\t'  read -r sample type rep time treat End1 End2; do
	echo "working on $sample using $End1 and $End2"
	bowtie2 -p 64 -x Tcongo_index -1 "$End1" -2 "$End2" 2> log/${sample}_align.log | 
	samtools sort -@ 32 -o aligned/${sample}.sorted.bam
done

for bam in aligned/*.sorted.bam
do
    sample=$(basename "$bam" .sorted.bam)
    bedtools bamtobed -i "$bam" > aligned/${sample}.aligned.bed
done


mkdir -p counts

for bam in aligned/*.sorted.bam
do
    sample=$(basename "$bam" .sorted.bam)
    bedtools intersect -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b "$bam" -c > counts/${sample}_counts.txt
done




mkdir -p means



tail -n +2 Tco2.fqfiles | cut -f2,4,5 | sort -u | while read type time treat
do
        files=$(grep -P "\t${type}\t.*\t${time}\t${treat}\t" Tco2.fqfiles | cut -f1 | sed "s/^/counts\//; s/$/_counts.txt/")
        n=$(echo "$files" | wc -w)
        outfile="means/${type}_${time}_${treat}_mean.txt"

        echo -e "gene\tdescription\tmean_count" > "$outfile"

        paste $files | awk -v n="$n" -F'\t' 'BEGIN{OFS="\t"} {
                sum = 0;
                gene = $4;
                description = $5;
                for (i = 6; i <= NF; i += 6) {
                        sum += $i
               }
                print gene, description, sum / n
        }' >> "$outfile"
done


mkdir -p comparisons

tail -n +2 Tco2.fqfiles | cut -f2,4 | sort -u | while read cell_line time
do
    induced_file="means/${cell_line}_${time}_Induced_mean.txt"
    uninduced_file="means/${cell_line}_${time}_Uninduced_mean.txt"

    outfile="comparisons/${cell_line}_${time}_Induced_vs_Uninduced.txt"

    echo -e "gene\tdescription\tfold_change" > "$outfile"

    paste "$induced_file" "$uninduced_file" | awk -F'\t' 'BEGIN{OFS="\t"} NR > 1 {
        gene = $1;
        description = $2;
        induced_mean = $3;
        uninduced_mean = $6;

        fc = induced_mean / (uninduced_mean + 0.0001);

        print gene, description, fc
    }' | sort -k3,3nr >> "$outfile"
done



tail -n +2 Tco2.fqfiles | cut -f2,5 | sort -u | while read cell_line induction
do
    file_0="means/${cell_line}_0_${induction}_mean.txt"

    file_24="means/${cell_line}_24_${induction}_mean.txt"
    outfile_24="comparisons/${cell_line}_${induction}_24_vs_0.txt"

    echo -e "gene\tdescription\tfold_change" > "$outfile_24"

    paste "$file_24" "$file_0" | awk -F'\t' 'BEGIN{OFS="\t"} NR > 1 {
        gene = $1;
        description = $2;
        mean_24 = $3;
        mean_0 = $6;

        fc = mean_24 / (mean_0 + 0.0001);

        print gene, description, fc
    }' | sort -k3,3nr >> "$outfile_24"

    file_48="means/${cell_line}_48_${induction}_mean.txt"
    outfile_48="comparisons/${cell_line}_${induction}_48_vs_0.txt"

    echo -e "gene\tdescription\tfold_change" > "$outfile_48"

    paste "$file_48" "$file_0" | awk -F'\t' 'BEGIN{OFS="\t"} NR > 1 {
        gene = $1;
        description = $2;
        mean_48 = $3;
        mean_0 = $6;

        fc = mean_48 / (mean_0 + 0.0001);

        print gene, description, fc
    }' | sort -k3,3nr >> "$outfile_48"
done
