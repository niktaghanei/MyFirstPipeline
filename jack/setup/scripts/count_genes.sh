#!/bin/bash

ALIGN_DIR="$1"

source ./config.txt || exit 1

# Keep counts beside the alignments from this run.
COUNT_DIR="$ALIGN_DIR/counts"
mkdir -p "$COUNT_DIR" || exit 1

total_samples=$(awk 'END {print NR}' "$ALIGN_DIR/samples.tsv")
sample_count=0

while read sample end1 end2; do

    if [[ "$MODE" == "interactive" ]]; then
        printf '\n\033[1;36mCounting gene overlaps for %s' "$sample" >&2
    fi

    start_time=$SECONDS

    # Append an overlap count to each gene annotation.
    bedtools intersect \
        -a "$GENE_BED" \
        -b "$ALIGN_DIR/$sample.sorted.bam" \
        -c \
        > "$COUNT_DIR/$sample.gene_counts.bed" \
        2> "$COUNT_DIR/$sample.counting.log" &

    counting_pid=$!

    if [[ "$MODE" == "interactive" ]]; then
        while kill -0 "$counting_pid" 2>/dev/null; do
            printf '.' >&2
            sleep 1
        done
    fi

    wait "$counting_pid"
    counting_status=$?
    elapsed=$((SECONDS - start_time))

    if [[ "$MODE" == "interactive" ]]; then
        printf ' (%ds)\033[0m\n' "$elapsed" >&2
    fi

    awk -F '\t' '
    BEGIN {
        OFS = "\t"
        print "GeneID", "Description", "ReadAlignmentCount"
    }
    {
        print $4, $5, $NF
    }
    ' "$COUNT_DIR/$sample.gene_counts.bed" \
        > "$COUNT_DIR/$sample.counts.tsv" || exit 1

    # Rank the same data by count, keeping the header at the top.
    printf 'GeneID\tDescription\tReadAlignmentCount\n' \
        > "$COUNT_DIR/$sample.counts_ranked.tsv" || exit 1

    awk 'NR > 1' "$COUNT_DIR/$sample.counts.tsv" |
        sort -t $'\t' -k3,3nr -k1,1 \
        >> "$COUNT_DIR/$sample.counts_ranked.tsv"

    sample_count=$((sample_count + 1))

    if [[ "$MODE" == "interactive" ]]; then
        percent=$((100 * sample_count / total_samples))

        printf '\033[1;32m✅ %s gene counts saved — %d of %d (%d%% complete)\033[0m\n' \
            "$sample" "$sample_count" "$total_samples" "$percent" >&2
    fi

done < "$ALIGN_DIR/samples.tsv"

printf '\nGene counting completed for %d samples.\n' "$sample_count"
printf 'Counts directory: %s\n' "$COUNT_DIR"
