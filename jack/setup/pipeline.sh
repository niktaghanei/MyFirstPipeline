#!/bin/bash

source ./config.txt

if [[ "$MODE" != "interactive" && "$MODE" != "quiet" ]]; then
    printf 'MODE is not set to interactive or quiet - defaulting to quiet mode.\n' >&2
    MODE="quiet"
fi

mkdir -p "$WORK_DIR"

RUN_ID=$(date '+%Y%m%d_%H%M%S')
REPORT_DIR="$WORK_DIR/reports/$RUN_ID"
mkdir -p "$REPORT_DIR"
QUALITY_REPORT="$REPORT_DIR/quality_report.md"

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;33mInitialising pipeline'
    for dot in 1 2 3 4 5; do
        sleep 0.2
        printf '.'
    done
    printf '\033[0m\n'

    printf '\n\033[1;36m'
    printf '%s\n' \
        '--------------------------------------------------' \
        '  01 — INPUT CHECKS AND EXPERIMENTAL COMPARISONS' \
        '--------------------------------------------------'
    printf '\033[0m\n'
fi

{
    printf '# Input validation and sequencing quality report\n\n'
    printf 'Run date: %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
    printf '## Input checks and planned comparisons\n\n```\n'
} > "$QUALITY_REPORT"

# Generate groups and comparisons.
if [[ "$MODE" == "interactive" ]]; then
    bash ./scripts/prepare_comparisons.sh \
        "$METADATA_FILE" "$WORK_DIR" "$FASTQ_DIR" |
        tee -a "$QUALITY_REPORT"
else
    bash ./scripts/prepare_comparisons.sh \
        "$METADATA_FILE" "$WORK_DIR" "$FASTQ_DIR" \
        >> "$QUALITY_REPORT"
fi

# Check the expected sequencing filenames.
if [[ "$MODE" == "interactive" ]]; then
    bash ./scripts/check_input_files.sh \
        "$METADATA_FILE" "$WORK_DIR" "$FASTQ_DIR" |
        tee -a "$QUALITY_REPORT"
else
    bash ./scripts/check_input_files.sh \
        "$METADATA_FILE" "$WORK_DIR" "$FASTQ_DIR" \
        >> "$QUALITY_REPORT"
fi

printf '\n```\n\n## FastQC results\n\n```text\n' \
    >> "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;32m✅ Input checks and comparison preparation complete.\033[0m\n'
    printf 'Press Enter to continue to FastQC, or Ctrl+C to stop: '
    read
    printf '\n'
fi

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;36m'
    printf '%s\n' \
        '--------------------------------------------------' \
        '  02 — SEQUENCING QUALITY: FASTQC' \
        '--------------------------------------------------'
    printf '\033[0m\n'
fi

if [[ "$MODE" == "interactive" ]]; then
    bash ./scripts/fastqc.sh \
        "$FASTQ_DIR" "$WORK_DIR" "$THREADS" \
        "$METADATA_FILE" "$MODE" |
        tee -a "$QUALITY_REPORT"
else
    bash ./scripts/fastqc.sh \
        "$FASTQ_DIR" "$WORK_DIR" "$THREADS" \
        "$METADATA_FILE" "$MODE" \
        >> "$QUALITY_REPORT"
fi

# Close the FastQC section in the Markdown report.
printf '\n```\n\n' >> "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;32m✅ FastQC collection and sample assessment complete.\033[0m\n'
    printf 'Press Enter to continue to alignment, or Ctrl+C to stop: '
    read
    printf '\n'
fi

# Read the directory recorded by the successful FastQC stage.
FASTQC_EXTRACT_DIR=$(cat "$WORK_DIR/current_fastqc_extract_dir.txt") || exit 1

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;36m'
    printf '%s\n' \
        '--------------------------------------------------' \
        '  03 — ALIGN READS TO THE GENOME' \
        '--------------------------------------------------'
    printf '\033[0m\n'
fi

printf '## Genome alignment\n\n```text\n' >> "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    bash ./scripts/align_reads.sh \
        "$FASTQC_EXTRACT_DIR/alignment_samples.tsv" |
        tee -a "$QUALITY_REPORT"
else
    bash ./scripts/align_reads.sh \
        "$FASTQC_EXTRACT_DIR/alignment_samples.tsv" \
        >> "$QUALITY_REPORT"
fi

printf '\n```\n\n' >> "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;32m✅ Genome alignment complete.\033[0m\n'
    printf 'Press Enter to continue to gene counting, or Ctrl+C to stop: '
    read
    printf '\n'
fi

# Read the directory recorded by the successful alignment stage.
ALIGN_DIR=$(cat "$WORK_DIR/current_alignment_dir.txt") || exit 1

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;36m'
    printf '%s\n' \
        '--------------------------------------------------' \
        '  04 — COUNT READ ALIGNMENTS OVERLAPPING GENES' \
        '--------------------------------------------------'
    printf '\033[0m\n'
fi

printf '## Gene counts\n\n```text\n' >> "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    bash ./scripts/count_genes.sh "$ALIGN_DIR" |
        tee -a "$QUALITY_REPORT"
else
    bash ./scripts/count_genes.sh "$ALIGN_DIR" \
        >> "$QUALITY_REPORT"
fi

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;32m✅ Labelled genes counted.\033[0m\n'
fi

printf 'Report saved at: %s\n' "$QUALITY_REPORT"
printf 'Gene counts saved at: %s/counts\n' "$ALIGN_DIR"

if [[ "$MODE" == "interactive" ]]; then
    printf '\nPress Enter to continue: '
    read
fi
