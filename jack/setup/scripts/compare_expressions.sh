#!/bin/bash

ALIGN_DIR="$1"
source ./config.txt

COUNT_DIR="$ALIGN_DIR/counts"
COMPARISONS_FILE="$WORK_DIR/prepare_inputs/comparisons.tsv"
EXPRESSION_DIR="$COUNT_DIR/expression"

mkdir -p "$EXPRESSION_DIR"


# Combine the sample count tables.
printf 'SampleName\tGeneID\tDescription\tReadAlignmentCount\n' \
    > "$EXPRESSION_DIR/combined_counts.tsv"

while read sample end1 end2; do
    file="$COUNT_DIR/$sample.counts.tsv"

    # if the count is missing for a file then alert as the mean would be affected
    if [[ ! -s "$file" ]]; then
        printf 'Missing or empty count table: %s\n' "$file" >&2
        exit 1
    fi

    awk -F '\t' -v sample="$sample" '
    BEGIN {
        OFS = "\t"
    }

    # Skip the file headers and combine all the counts into one tsv
    NR > 1 && NF {
        print sample, $1, $2, $3
    }
    ' "$file" >> "$EXPRESSION_DIR/combined_counts.tsv"

done < "$ALIGN_DIR/samples.tsv"


# Calculate the group means - i.e. the average number of counts per gene for each experimental group (combination of sample type, time, treatment)
awk -F '\t' '
BEGIN {
    OFS = "\t"

    print "GeneID", "Description", "SampleType", "Time",
          "Treatment", "SampleCount", "MeanCount"
}

# First file: store each samples experimental conditions.
FILENAME == ARGV[1] {
    if (FNR > 1 && NF) {
        sample_type[$1] = $2
        time[$1] = $4
        treatment[$1] = $5
    }
    next
}

# Second file: the combined count table.
FNR == 1 { next }
NF == 0 { next }

{
    sample = $1
    gene = $2

    # check the arrays stored from the first input file to get the sample group
    cell = sample_type[sample]
    hour = time[sample]
    induction = treatment[sample]

    #count each sample once to get the denominator - if it is not already in the counted array add it with the experimental conditions as the key and increment a counter
    if (!(sample in counted)) {
        sample_count[cell, hour, induction]++
        counted[sample] = 1
    }

    description[gene] = $3

    # Add the number in this sample to an array tagged with its gene and group - this is += so all the samples for the same gene and experimental group are added together. these will make up the numerator
    total[gene, cell, hour, induction] += $4
}

END {
    for (key in total) {
        # Recover the four parts of the combined array key using subsep to separate the key
        split(key, parts, SUBSEP)

        gene = parts[1]
        cell = parts[2]
        hour = parts[3]
        induction = parts[4]

        number = sample_count[cell, hour, induction]
        mean = total[key] / number

        print gene, description[gene], cell, hour,
              induction, number, mean
    }
}
' "$METADATA_FILE" "$EXPRESSION_DIR/combined_counts.tsv" \
    > "$EXPRESSION_DIR/group_means.tsv"


# get the list of comparisons from step 1
awk 'NR > 1 && NF { print }' "$COMPARISONS_FILE" \
    > "$EXPRESSION_DIR/comparisons_to_run.tsv" || exit 1

# prepare a file to store the fold changes
printf '' > "$EXPRESSION_DIR/fold_change_rows.tsv"

completed=0

while read id name type num_type num_time num_treatment den_type den_time den_treatment; do

    if [[ "$MODE" == "interactive" ]]; then
        printf '\nComparing %s\n' "$name" >&2
    fi

    # Pass the comparison identifiers and group definitions to AWK using -v for variable
    awk -F '\t' \
        -v comparison_id="$id" \
        -v comparison="$name" \
        -v comparison_type="$type" \
        -v num_type="$num_type" \
        -v num_time="$num_time" \
        -v num_treatment="$num_treatment" \
        -v den_type="$den_type" \
        -v den_time="$den_time" \
        -v den_treatment="$den_treatment" '
    BEGIN {
        OFS = "\t"
    }

    NR == 1 { next }

    {
        # Store the means for the numerator group
        if ($3 == num_type &&
            $4 == num_time &&
            $5 == num_treatment) {

            numerator[$1] = $7
            numerator_samples = $6
            description[$1] = $2
            genes[$1] = 1
        }

        # Store the means for the denominator group
        if ($3 == den_type &&
            $4 == den_time &&
            $5 == den_treatment) {

            denominator[$1] = $7
            denominator_samples = $6
            description[$1] = $2
            genes[$1] = 1
        }
    }

    END {
        # Calculate one result for every gene in this comparison group
        for (gene in genes) {
            num = numerator[gene]
            den = denominator[gene]

            # deal with zeros - flag if one or the other is zero and output NAs for the metrics
            if (num == 0 || den == 0) {
                ratio = "NA"
                direction = "NA"
                magnitude = "NA"

                if (num == 0 && den == 0) {
                    status = "BOTH_ZERO"
                } else if (den == 0) {
                    status = "DENOMINATOR_ZERO"
                    direction = "Higher"
                } else {
                    status = "NUMERATOR_ZERO"
                    direction = "Lower"
                    ratio = 0
                }

            } else {
                # Numerator divided by denominator.
                ratio = num / den
                status = "VALID"

                # scale the fold magnitude to account for decreases as well as increases
		# record the direction in a separate variable (and column) 
                if (ratio > 1) {
                    direction = "Higher"
                    magnitude = ratio
                } else if (ratio < 1) {
                    direction = "Lower"
                    magnitude = 1 / ratio
                } else {
                    direction = "Unchanged"
                    magnitude = 1
                }
            }

            # Write the result, including which comparison group
            print comparison_id, comparison, comparison_type,
                  gene, description[gene], num, den,
                  numerator_samples, denominator_samples,
                  ratio, direction, magnitude, status
        }
    }
    ' "$EXPRESSION_DIR/group_means.tsv" \
        >> "$EXPRESSION_DIR/fold_change_rows.tsv"

    completed=$((completed + 1))

done < "$EXPRESSION_DIR/comparisons_to_run.tsv"


# Write one header for the combined ranked table.
printf 'ComparisonID\tComparisonName\tComparisonType\tGeneID\tDescription\tNumeratorMean\tDenominatorMean\tNumeratorSamples\tDenominatorSamples\tFoldChange\tDirection\tFoldMagnitude\tStatus\n' \
    > "$EXPRESSION_DIR/fold_changes_ranked.tsv"

# Rank valid results by fold magnitude (column 12)
awk -F '\t' '$13 == "VALID"' \
    "$EXPRESSION_DIR/fold_change_rows.tsv" |
    sort -t $'\t' -k12,12nr -k4,4 -k1,1 \
    >> "$EXPRESSION_DIR/fold_changes_ranked.tsv"

# Add the zeros in as well
awk -F '\t' '$13 != "VALID"' \
    "$EXPRESSION_DIR/fold_change_rows.tsv" |
    sort -t $'\t' -k4,4 -k1,1 \
    >> "$EXPRESSION_DIR/fold_changes_ranked.tsv"


printf '\nExpression comparisons complete.\n'
printf '  Comparisons calculated : %d\n' "$completed"
printf 'Group means: %s/group_means.tsv\n' "$EXPRESSION_DIR"
printf 'Fold changes: %s/fold_changes_ranked.tsv\n' "$EXPRESSION_DIR"






















