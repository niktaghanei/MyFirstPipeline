#!/bin/bash

ALIGNMENT_SAMPLES_FILE="$1"

source ./config.txt || exit 1

# Allow a thread override in interactive mode.
if [[ "$MODE" == "interactive" ]]; then
    printf '\nAlignment threads are set to %s.\n' "$THREADS" >&2
    printf 'Enter a new value, or press Enter to keep this setting: ' >&2
    read override

    if [[ -n "$override" ]]; then
        THREADS="$override"
    fi
fi

RUN_ID=$(date '+%Y%m%d_%H%M%S')
ALIGN_DIR="$WORK_DIR/alignments/$RUN_ID"
INDEX_DIR="$ALIGN_DIR/bowtie2_index"
INDEX_PREFIX="$INDEX_DIR/Tcongolense_IL3000"

mkdir -p "$INDEX_DIR" || exit 1

awk -F '\t' 'NR > 1 && NF {print $0}' \
    "$ALIGNMENT_SAMPLES_FILE" \
    > "$ALIGN_DIR/samples.tsv" || exit 1

total_samples=$(awk 'END {print NR}' "$ALIGN_DIR/samples.tsv")

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;36mBuilding the Bowtie2 genome index' >&2
fi

start_time=$SECONDS

bowtie2-build "$GENOME_FASTA" "$INDEX_PREFIX" \
    > "$ALIGN_DIR/index.log" 2>&1 &

index_pid=$!

if [[ "$MODE" == "interactive" ]]; then
    while kill -0 "$index_pid" 2>/dev/null; do
        printf '.' >&2
        sleep 1
    done
fi

wait "$index_pid"
index_status=$?
elapsed=$((SECONDS - start_time))

if [[ "$MODE" == "interactive" ]]; then
    printf ' (%ds)\033[0m\n' "$elapsed" >&2
fi

sample_count=0

while read sample end1 end2; do

    if [[ "$MODE" == "interactive" ]]; then
        printf '\n\033[1;33mAligning %s using %s threads' \
            "$sample" "$THREADS" >&2
    fi

    start_time=$SECONDS

    bowtie2 \
        -p "$THREADS" \
        -x "$INDEX_PREFIX" \
        -1 "$FASTQ_DIR/$end1" \
        -2 "$FASTQ_DIR/$end2" \
        -S "$ALIGN_DIR/$sample.sam" \
        > "$ALIGN_DIR/$sample.bowtie2.log" 2>&1 &

    alignment_pid=$!

    if [[ "$MODE" == "interactive" ]]; then
        while kill -0 "$alignment_pid" 2>/dev/null; do
            printf '.' >&2
            sleep 1
        done
    fi

    wait "$alignment_pid"
    alignment_status=$?
    elapsed=$((SECONDS - start_time))

    if [[ "$MODE" == "interactive" ]]; then
        printf ' (%ds)\033[0m\n' "$elapsed" >&2
    fi

    # Convert SAM to BAM and sort by genomic position.
    samtools sort \
        -o "$ALIGN_DIR/$sample.sorted.bam" \
        "$ALIGN_DIR/$sample.sam" \
        > "$ALIGN_DIR/$sample.samtools.log" 2>&1

    # Create the accompanying BAM index.
    samtools index "$ALIGN_DIR/$sample.sorted.bam" \
        >> "$ALIGN_DIR/$sample.samtools.log" 2>&1

    # Save alignment statistics for later assessment.
    samtools flagstat "$ALIGN_DIR/$sample.sorted.bam" \
        > "$ALIGN_DIR/$sample.flagstat.txt" \
        2>> "$ALIGN_DIR/$sample.samtools.log"

    sample_count=$((sample_count + 1))

    if [[ "$MODE" == "interactive" ]]; then
        percent=$((100 * sample_count / total_samples))

        printf '\033[1;32m✅ %s aligned, sorted and indexed — %d of %d (%d%% complete)\033[0m\n' \
            "$sample" "$sample_count" "$total_samples" "$percent" >&2
    fi

done < "$ALIGN_DIR/samples.tsv"

printf '\nAlignment completed for %d samples.\n' "$sample_count"
printf 'Results directory: %s\n' "$ALIGN_DIR"

# Record this successful run for gene counting.
printf '%s\n' "$ALIGN_DIR" \
    > "$WORK_DIR/current_alignment_dir.txt" || exit 1
