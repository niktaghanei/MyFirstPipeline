#!/bin/bash

if [[ $# -ne 5 ]]; then
    printf 'Usage: fastqc.sh FASTQ_DIR WORK_DIR THREADS METADATA_FILE MODE\n' >&2
    exit 1
fi

FASTQ_DIR="$1"
WORK_DIR="$2"
THREADS="$3"
METADATA_FILE="$4"
MODE="$5"

EXPECTED_FILES="$WORK_DIR/prepare_inputs/expected_files.txt"


if [[ "$MODE" == "interactive" ]]; then
    printf '\nFastQC threads are set to %s.\n' "$THREADS" >&2
    printf 'Enter a new value, or press Enter to keep this setting: ' >&2
    read override

    if [[ -n "$override" ]]; then
        THREADS="$override"
    fi
fi

RUN_ID=$(date '+%Y%m%d_%H%M%S')
FASTQC_RUN_DIR="$WORK_DIR/fastqc/$RUN_ID"
FASTQC_REPORT_DIR="$FASTQC_RUN_DIR/reports"
FASTQC_EXTRACT_DIR="$FASTQC_RUN_DIR/extracted"

mkdir -p "$FASTQC_REPORT_DIR" "$FASTQC_EXTRACT_DIR"

fastqc_files=()

while read filename; do
    if [[ -n "$filename" ]]; then
        fastqc_files+=("$FASTQ_DIR/$filename")
    fi
done < "$EXPECTED_FILES"


if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;33mProcessing FASTQ files using FastQC' >&2
fi

fastqc --quiet --threads "$THREADS" \
    --outdir "$FASTQC_REPORT_DIR" "${fastqc_files[@]}" \
    > "$FASTQC_RUN_DIR/fastqc.log" 2>&1 &

fastqc_pid=$!

# Show progress only in interactive mode.
if [[ "$MODE" == "interactive" ]]; then
    while kill -0 "$fastqc_pid" 2>/dev/null; do
        printf '.' >&2
        sleep 1
    done
    printf '\033[0m\n' >&2
fi

wait "$fastqc_pid"
fastqc_status=$?

if [[ "$fastqc_status" -ne 0 ]]; then
    printf 'FastQC failed. See %s/fastqc.log\n' "$FASTQC_RUN_DIR" >&2
    exit 1
fi

# Initialise the combined files
printf '' > "$FASTQC_EXTRACT_DIR/combined_summaries.txt"

printf 'Filename\tTotalReads\tReadLength\tGCPercent\tMinimumMedianQuality\tMinimumLowerQuartileQuality\tModeSequenceQuality\tMaximumNPercentage\n' \
    > "$FASTQC_EXTRACT_DIR/combined_metrics.txt" || exit 1

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;36mExtracting reports and collecting metrics' >&2
fi

file_count=0

while read filename; do
    [[ -z "$filename" ]] && continue

    name="${filename%.fq.gz}"
    report="$FASTQC_REPORT_DIR/${name}_fastqc.zip"
    result_dir="$FASTQC_EXTRACT_DIR/${name}_fastqc"

    if ! unzip -q "$report" -d "$FASTQC_EXTRACT_DIR" \
        >> "$FASTQC_RUN_DIR/extraction.log" 2>&1; then
        printf 'Extraction failed for %s. See %s/extraction.log\n' \
            "$filename" "$FASTQC_RUN_DIR" >&2
        exit 1
    fi

    cat "$result_dir/summary.txt" \
        >> "$FASTQC_EXTRACT_DIR/combined_summaries.txt"

    awk -F '\t' '
    BEGIN { OFS = "\t" }

    # Basic statistics.
    $1 == "Filename" { filename = $2 }
    $1 == "Total Sequences" { total_reads = $2 }
    $1 == "Sequence length" { read_length = $2 }
    $1 == "%GC" { gc_percent = $2 }

    # Identify each relevant module.
    $1 == ">>Per base sequence quality" {
        in_quality = 1
        next
    }
    $1 == ">>Per sequence quality scores" {
        in_read_quality = 1
        next
    }
    $1 == ">>Per base N content" {
        in_N_content = 1
        next
    }

    $1 == ">>END_MODULE" {
        in_quality = 0
        in_read_quality = 0
        in_N_content = 0
    }

    # Minimum median and lower-quartile base quality.
    in_quality == 1 && $1 != "#Base" {
        quality_rows++

        if (quality_rows == 1 || $3 < minimum_median) {
            minimum_median = $3
        }
        if (quality_rows == 1 || $4 < minimum_lower_quartile) {
            minimum_lower_quartile = $4
        }
    }

    # Most frequent mean quality score across reads.
    in_read_quality == 1 && $1 != "#Quality" {
        read_quality_rows++

        if (read_quality_rows == 1 || $2 > largest_count) {
            largest_count = $2
            mode_mean_quality = $1
        }
    }

    # Highest N percentage across positions or position bins.
    in_N_content == 1 && $1 != "#Base" {
        N_content_rows++

        if (N_content_rows == 1 || $2 > largest_N_percentage) {
            largest_N_percentage = $2
        }
    }

    END {
        print filename, total_reads, read_length, gc_percent,
              minimum_median, minimum_lower_quartile,
              mode_mean_quality, largest_N_percentage
    }
    ' "$result_dir/fastqc_data.txt" \
        >> "$FASTQC_EXTRACT_DIR/combined_metrics.txt" || exit 1

    file_count=$((file_count + 1))

    if [[ "$MODE" == "interactive" ]]; then
        printf '.' >&2
    fi
done < "$EXPECTED_FILES"

if [[ "$MODE" == "interactive" ]]; then
    printf '\033[0m\n' >&2
fi

# Join R1 and R2 metrics to their sample names.
bash ./scripts/analyse_fastq.sh \
    "$METADATA_FILE" "$FASTQC_EXTRACT_DIR" || exit 1

printf '\nProcessed and collected FastQC results for %d files.\n' "$file_count"
printf 'Summaries: %s/combined_summaries.txt\n' "$FASTQC_EXTRACT_DIR"
printf 'Metrics: %s/combined_metrics.txt\n' "$FASTQC_EXTRACT_DIR"
printf 'Paired metrics: %s/paired_metrics.tsv\n' "$FASTQC_EXTRACT_DIR"

# Record this successful run for the next pipeline stage.
printf '%s\n' "$FASTQC_EXTRACT_DIR" \
    > "$WORK_DIR/current_fastqc_extract_dir.txt" || exit 1

# Assess the outputs from this FastQC run.
bash ./scripts/assess_fastqc.sh \
    "$METADATA_FILE" "$FASTQC_EXTRACT_DIR" || exit 1












