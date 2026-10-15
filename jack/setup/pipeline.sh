#!/bin/bash

source ./config.txt


#check which mode it is in - interactive mode prints to the screen and allows override of variables
if [[ "$MODE" != "interactive" && "$MODE" != "quiet" ]]; then
    printf 'MODE is not set to interactive or quiet - defaulting to quiet mode.\n' >&2
    MODE="quiet"
fi

mkdir -p "$WORK_DIR"

# save in datetime stamped folders so different runs do not get mixed up
RUN_ID=$(date '+%Y%m%d_%H%M%S')
REPORT_DIR="$WORK_DIR/reports/$RUN_ID"
mkdir -p "$REPORT_DIR"

#save the outputs of the run to a quality report in markdown for later review
QUALITY_REPORT="$REPORT_DIR/quality_report.md"

if [[ "$MODE" == "interactive" ]]; then
    #change the colour to yellow and print some dots to show the pipeline starting
    printf '\n\033[1;33mInitialising pipeline'
    for dot in 1 2 3 4 5; do
        sleep 0.2
        printf '.'
    done
    #reset the colour
    printf '\033[0m\n'

    #set the colour to cyan and output the header for the first pipeline stage
    printf '\n\033[1;36m'
    printf '%s\n' \
        '--------------------------------------------------' \
        '  01 — INPUT CHECKS AND EXPERIMENTAL COMPARISONS' \
        '--------------------------------------------------'
    printf '\033[0m\n'
fi

# set up the quality report with the title and date
{
    printf '# Input validation and sequencing quality report\n\n'
    printf 'Run date: %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
    printf '## Input checks and planned comparisons\n\n```\n'
} > "$QUALITY_REPORT"

#generate groups and comparisons
if [[ "$MODE" == "interactive" ]]; then
    bash ./scripts/prepare_comparisons.sh \
        "$METADATA_FILE" "$WORK_DIR" "$FASTQ_DIR" |
        # tee will put the results that are sent to the screen also in the quality report
  	tee -a "$QUALITY_REPORT"
else
    bash ./scripts/prepare_comparisons.sh \
	#if in quiet mode then just output the results into the quality report
        "$METADATA_FILE" "$WORK_DIR" "$FASTQ_DIR" \
        >> "$QUALITY_REPORT"
fi

# ensure the expected sequencing filenames are present as specified in the metadata file
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

#run the fastqc script to generate the quality reports
#this script also runs analyse_fastq.sh to extract the metrics from the FastQC reports and assess_fastqc.sh to compare the results to the thresholds specified in config.txt
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

# Close the FastQC section in the .md report
printf '\n```\n\n' >> "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;32m✅ FastQC collection and sample assessment complete.\033[0m\n'
    printf 'Press Enter to continue to alignment, or Ctrl+C to stop: '
    read
    printf '\n'
fi

#store the directory recorded by the successful FastQC stage for this run so it can be passed to later steps in the pathway
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

# Read the directory recorded by the successful alignment stage for passing to later pathway steps
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

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;36m'
    printf '%s\n' \
        '--------------------------------------------------' \
        '  05 — GROUP MEANS AND FOLD CHANGES' \
        '--------------------------------------------------'
    printf '\033[0m\n'
fi

printf '## Group means and fold changes\n\n```text\n' >> "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    bash ./scripts/compare_expressions.sh "$ALIGN_DIR" |
        tee -a "$QUALITY_REPORT" || exit 1
else
    bash ./scripts/compare_expressions.sh "$ALIGN_DIR" \
        >> "$QUALITY_REPORT" || exit 1
fi

printf '\n```\n\n' >> "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    printf '\n\033[1;32m✅ Group means and fold changes calculated.\033[0m\n'
fi

printf 'Expression results: %s/counts/expression\n' "$ALIGN_DIR"
printf 'Report saved at: %s\n' "$QUALITY_REPORT"

if [[ "$MODE" == "interactive" ]]; then
    printf '\nPress Enter to finish: '
    read
fi














