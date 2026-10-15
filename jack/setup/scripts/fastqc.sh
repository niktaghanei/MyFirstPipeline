#!/bin/bash

FASTQ_DIR="$1"
WORK_DIR="$2"
THREADS="$3"
METADATA_FILE="$4"
MODE="$5"

EXPECTED_FILES="$WORK_DIR/prepare_inputs/expected_files.txt"

# Allow a key interrupt in interactive mode where the user can specify the number of threads during the run.
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

# get a list of the files to pass through fastqc. the expected files are from the metadata file so contain the complete list needed for the analysis
# they have already passed checking to ensure they all exist in the folder so this should not cause any errors
while read filename; do
    if [[ -n "$filename" ]]; then
        fastqc_files+=("$FASTQ_DIR/$filename")
    fi
done < "$EXPECTED_FILES"


if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;33mProcessing FASTQ files using FastQC' >&2 
    #the output >&2 is used to put some printed statements into the error output so that when we call tee on the script to generate the markdown report these statements will not go into the markdown file and will just print to the screen
    # a number of the screen messages are per file which would clutter the markdown report if all outputs mirrored the screen
fi

# pass fastqc the array of expected files so it runs them all together with parallel processing
fastqc --quiet --threads "$THREADS" \
    --outdir "$FASTQC_REPORT_DIR" "${fastqc_files[@]}" \
    # Put the run output into a log so we can see if a particular file failed
    > "$FASTQC_RUN_DIR/fastqc.log" 2>&1 & # use & to run this in the background so we can output something to the screen at the same time so the user knows it is still working

#$! is the process ID of the last command run in the background i.e. the fastqc
fastqc_pid=$!

# output dots while the fastqc is still running
if [[ "$MODE" == "interactive" ]]; then
    # kill -0 tests whether the process is still running
    #2>/dev/null hides the kill error once it is finished
    while kill -0 "$fastqc_pid" 2>/dev/null; do
        printf '.' >&2
        sleep 1
    done
    printf '\033[0m\n' >&2
fi

#Now extract the results from the outputs of the fastqc process
#it produces a zip file for each extract which contains a summary file with the inbuilt quality test pass/fail results and another file with detailed numerical breakdowns

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

    # if it doesn't unzip then add it to the log
    if ! unzip -q "$report" -d "$FASTQC_EXTRACT_DIR" \
        >> "$FASTQC_RUN_DIR/extraction.log" 2>&1; then
        printf 'Extraction failed for %s. See %s/extraction.log\n' \
            "$filename" "$FASTQC_RUN_DIR" >&2
        exit 1
    fi

    # Otherwise put the whole contents into the combined summary file
    cat "$result_dir/summary.txt" \
        >> "$FASTQC_EXTRACT_DIR/combined_summaries.txt"

    # In order to read the numerical metrics we need to determine whether we are in the right 'module' of output
    # modules are started with >>Module Name and ended with >>END_MODULE
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

    # reset when we reach the end of a module
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

# Join R1 and R2 metrics to their sample names so we can make assessments of the quality based on both strands at the same time
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












