#!/bin/bash

METADATA_FILE="$1"
FASTQC_EXTRACT_DIR="$2"

# Run from setup, as with the other scripts.
source ./config.txt || exit 1

FLAGS_FILE="$FASTQC_EXTRACT_DIR/fastqc_flags.tsv"
ASSESSMENT_FILE="$FASTQC_EXTRACT_DIR/sample_assessment.tsv"

# Match each FastQC warning or failure to its sample and mate.
awk -F '\t' '
BEGIN {
    OFS = "\t"
    print "SampleName", "Mate", "Filename", "Module", "Status"
}

FILENAME == ARGV[1] {
    if (FNR > 1 && NF) {
        sample[$6] = $1
        mate[$6] = "R1"

        sample[$7] = $1
        mate[$7] = "R2"
    }
    next
}

$1 == "WARN" || $1 == "FAIL" {
    if (!($3 in sample)) {
        print "FastQC filename not found in metadata: " $3 > "/dev/stderr"
        unmatched = 1
        next
    }

    print sample[$3], mate[$3], $3, $2, $1
}

END {
    if (unmatched) exit 1
}
' "$METADATA_FILE" "$FASTQC_EXTRACT_DIR/combined_summaries.txt" \
    > "$FLAGS_FILE" || exit 1


# Assess both mates using the configured numerical thresholds.
awk -F '\t' \
    -v median_flag="$QC_MEDIAN_FLAG" \
    -v median_review="$QC_MEDIAN_REVIEW" \
    -v quartile_flag="$QC_LOWER_QUARTILE_FLAG" \
    -v quartile_review="$QC_LOWER_QUARTILE_REVIEW" \
    -v mode_flag="$QC_MODE_FLAG" \
    -v mode_review="$QC_MODE_REVIEW" \
    -v n_flag="$QC_N_FLAG_PERCENT" \
    -v n_review="$QC_N_REVIEW_PERCENT" \
    -v align_review="$QC_ALIGN_REVIEW_SAMPLES" '
BEGIN {
    OFS = "\t"

    # Each metric has a name, column position and two thresholds.
    metric[1] = "MinimumMedianQuality"
    column[1] = 8
    flag[1] = median_flag
    review[1] = median_review

    metric[2] = "MinimumLowerQuartileQuality"
    column[2] = 10
    flag[2] = quartile_flag
    review[2] = quartile_review

    metric[3] = "ModeSequenceQuality"
    column[3] = 12
    flag[3] = mode_flag
    review[3] = mode_review

    metric[4] = "MaximumNPercentage"
    column[4] = 14
    flag[4] = n_flag
    review[4] = n_review

    print "SampleName", "Assessment", "AlignmentAction",
          "FastQCWarnings", "FastQCFailures", "Reasons"
}

# First file: count FastQC module warnings and failures per sample.
FILENAME == ARGV[1] {
    if (FNR > 1 && NF) {
        if ($5 == "WARN") warnings[$1]++
        if ($5 == "FAIL") failures[$1]++
    }
    next
}

# Second file: assess the paired metrics.
FNR == 1 { next }
NF == 0 { next }

{
    assessment = "PASS"
    reasons = ""

    # Check the read counts for both mates.
    for (mate_number = 1; mate_number <= 2; mate_number++) {
        read_column = mate_number + 1
        read_count = $read_column

        if (read_count !~ /^[0-9]+$/) {
            assessment = "UNUSABLE"
            reasons = reasons "R" mate_number ": missing or invalid read count; "
        } else if (read_count + 0 == 0) {
            assessment = "UNUSABLE"
            reasons = reasons "R" mate_number ": zero reads; "
        }
    }

    # Check each numerical quality metric for R1 and R2.
    for (i = 1; i <= 4; i++) {
        for (mate_number = 1; mate_number <= 2; mate_number++) {
            field_number = column[i] + mate_number - 1
            value = $field_number

            if (value !~ /^[0-9]+([.][0-9]+)?$/) {
                assessment = "UNUSABLE"
                reasons = reasons "R" mate_number " " metric[i] ": missing or invalid value; "
                continue
            }

            issue = ""
            threshold = ""

            # Lower quality scores are worse.
            if (i <= 3) {
                if (value + 0 < review[i] + 0) {
                    issue = "REVIEW"
                    threshold = review[i]
                } else if (value + 0 < flag[i] + 0) {
                    issue = "FLAG"
                    threshold = flag[i]
                }
            }

            # Higher N percentages are worse.
            if (i == 4) {
                if (value + 0 > review[i] + 0) {
                    issue = "REVIEW"
                    threshold = review[i]
                } else if (value + 0 > flag[i] + 0) {
                    issue = "FLAG"
                    threshold = flag[i]
                }
            }

            if (issue != "") {
                reasons = reasons "R" mate_number " " metric[i] \
                          "=" value " (" issue " threshold=" threshold "); "

                if (issue == "REVIEW" && assessment != "UNUSABLE") {
                    assessment = "REVIEW"
                } else if (issue == "FLAG" && assessment == "PASS") {
                    assessment = "FLAG"
                }
            }
        }
    }

    if (reasons == "") {
        reasons = "No thresholds from the config settings are crossed"
    } else {
        sub(/; $/, "", reasons)
    }

    # Module flags are recorded separately from numerical assessment.
    action = "PROCEED"

    if (assessment == "UNUSABLE") {
        action = "DO_NOT_ALIGN"
    } else if (assessment == "REVIEW") {
        if (align_review == "1") {
            action = "PROCEED_WITH_REVIEW"
        } else {
            action = "WITHHOLD_FOR_REVIEW"
        }
    }

    print $1, assessment, action,
          warnings[$1] + 0, failures[$1] + 0, reasons
}
' "$FLAGS_FILE" "$FASTQC_EXTRACT_DIR/paired_metrics.tsv" \
    > "$ASSESSMENT_FILE" || exit 1


# Include numerical issues and any FastQC module warnings/failures.
awk -F '\t' '
NR == 1 || $2 != "PASS" || $4 > 0 || $5 > 0
' "$ASSESSMENT_FILE" \
    > "$FASTQC_EXTRACT_DIR/samples_to_check.tsv" || exit 1

# Keep unusable inputs separate from samples awaiting review.
awk -F '\t' '
NR == 1 || $2 == "UNUSABLE"
' "$ASSESSMENT_FILE" \
    > "$FASTQC_EXTRACT_DIR/unusable_samples.tsv" || exit 1

# Write the sample name and mate filenames for permitted samples.
awk -F '\t' '
BEGIN {
    OFS = "\t"
    print "SampleName", "End1", "End2"
}

FILENAME == ARGV[1] {
    if (FNR > 1 &&
        ($3 == "PROCEED" || $3 == "PROCEED_WITH_REVIEW")) {
        permitted[$1] = 1
    }
    next
}

FNR == 1 { next }

$1 in permitted {
    print $1, $6, $7
}
' "$ASSESSMENT_FILE" "$METADATA_FILE" \
    > "$FASTQC_EXTRACT_DIR/alignment_samples.tsv" || exit 1

# Summarise the assessment.
awk -F '\t' '
NR == 1 { next }
NF {
    total++
    counts[$2]++

    if ($4 > 0 || $5 > 0) module_flags++

    if ($3 == "PROCEED" || $3 == "PROCEED_WITH_REVIEW") {
        permitted++
    }
    if ($3 == "WITHHOLD_FOR_REVIEW") withheld++
}

END {
    printf "\nSample quality assessment is complete.\n"
    printf "  Samples assessed                   : %d\n", total
    printf "  Numerical PASS                     : %d\n", counts["PASS"]
    printf "  Numerical FLAG                     : %d\n", counts["FLAG"]
    printf "  Numerical REVIEW                   : %d\n", counts["REVIEW"]
    printf "  Unusable inputs                    : %d\n", counts["UNUSABLE"]
    printf "  Samples with FastQC module flags   : %d\n", module_flags
    printf "  Permitted for alignment by policy  : %d\n", permitted
    printf "  Withheld pending review            : %d\n", withheld
}
' "$ASSESSMENT_FILE" || exit 1

printf '\nAssessment table: %s\n' "$ASSESSMENT_FILE"
printf 'FastQC flags: %s\n' "$FLAGS_FILE"
printf 'Samples to check: %s/samples_to_check.tsv\n' "$FASTQC_EXTRACT_DIR"
printf 'Unusable samples: %s/unusable_samples.tsv\n' "$FASTQC_EXTRACT_DIR"
printf 'Alignment sample list: %s/alignment_samples.tsv\n' "$FASTQC_EXTRACT_DIR"









