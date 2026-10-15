#!/bin/bash

METADATA_FILE="$1"
WORK_DIR="$2"
FASTQ_DIR="$3"
PREPARE_DIR="$WORK_DIR/prepare_inputs"

mkdir -p "$PREPARE_DIR" || exit 1

# this will dynamically discover and count the groups in the current metadata. This script needs to be dynamic so that if another clone is added or another experimental condition then all the experimental condition combinations will be generated to compare all possible scenarios
awk -F '\t' 'NR > 1 && NF {print $2, $4, $5}' "$METADATA_FILE" |
sort |
uniq -c |
awk '
BEGIN {
    OFS = "\t"
    print "SampleType", "Time", "Treatment", "SampleCount"
}
{
    print $2, $3, $4, $1
}
' > "$PREPARE_DIR/group_inventory.tsv" || exit 1
#group inventory is a list of the groups and how many replicates are in each group so we can see whether there are enough for valid comparison

# Read the group inventory back in to output a message about what has been found.
awk -F '\t' '
NR == 1 { next }
NF {
    groups++
    samples += $4
    if (groups == 1 || $4 < minimum) minimum = $4
    if (groups == 1 || $4 > maximum) maximum = $4
}
END {
    printf "\nSamples have been grouped by cell line, time point and treatment.\n"
    printf "  %-20s : %d\n", "Experimental groups", groups
    printf "  %-20s : %d\n", "Biological samples", samples
    printf "  %-20s : %d to %d\n", "Samples per group", minimum, maximum
}
' "$PREPARE_DIR/group_inventory.tsv" || exit 1

printf '\nGroup inventory file: %s\n' "$PREPARE_DIR/group_inventory.tsv"

# Generate a list of experimental conditions dynamically (based on what is in the metadata file) so later we can loop the analysis and compare which experimental conditions saw the biggest fold differences in gene expression
awk -F '\t' '
BEGIN { OFS = "\t" }

#skip headers and blanks
NR == 1 { next }
NF == 0 { next }

{
    if ($3 == "Induced") induced[$1, $2] = 1
    if ($3 == "Uninduced") uninduced[$1, $2] = 1

    times[$2] = 1
    combinations[$1, $3] = 1
    exists[$1, $3, $2] = 1
    cell_lines[$1] = 1
    time_treatments[$2, $3] = 1
}

END {
    # Induced versus uninduced at matching cell line and time.
    for (key in induced) {
        if (key in uninduced) {
	    # use subsep to split the array key into its parts as stored above
            split(key, parts, SUBSEP)
            sample_type = parts[1]
            time = parts[2]

            description = "Induction_" sample_type "_" time "h_Induced_vs_Uninduced"
            print description, "Induction",
                  sample_type, time, "Induced",
                  sample_type, time, "Uninduced"
        }
    }

    # Later versus earlier at matching cell line and treatment.
    for (combination in combinations) {
        split(combination, parts, SUBSEP)
        sample_type = parts[1]
        treatment = parts[2]

        for (later in times) {
            for (earlier in times) {
                if ((later + 0) > (earlier + 0) &&
                    ((sample_type, treatment, later) in exists) &&
                    ((sample_type, treatment, earlier) in exists)) {

                    description = "Timing_" sample_type "_" treatment "_" later "h_vs_" earlier "h"
                    print description, "Timing",
                          sample_type, later, treatment,
                          sample_type, earlier, treatment
                }
            }
        }
    }

    # Compare each pair of cell lines once at matching time and treatment.
    for (combination in time_treatments) {
        split(combination, parts, SUBSEP)
        time = parts[1]
        treatment = parts[2]

        for (numerator in cell_lines) {
            for (denominator in cell_lines) {
                if (numerator < denominator &&
                    ((numerator, treatment, time) in exists) &&
                    ((denominator, treatment, time) in exists)) {

                    description = "CellLine_" numerator "_vs_" denominator "_" time "h_" treatment
                    print description, "CellLine",
                          numerator, time, treatment,
                          denominator, time, treatment
                }
            }
        }
    }
}
' "$PREPARE_DIR/group_inventory.tsv" |
sort -t $'\t' \
    -k2,2 -k3,3 -k5,5 -k4,4n \
    -k6,6 -k8,8 -k7,7n -k1,1 |
awk '
BEGIN {
    OFS = "\t"
    print "ComparisonID", "ComparisonName", "ComparisonType",
          "NumeratorSampleType", "NumeratorTime", "NumeratorTreatment",
          "DenominatorSampleType", "DenominatorTime", "DenominatorTreatment"
}
{
    # generate a unique ID for each comparison scenario
    print NR + 1000, $0
}
' > "$PREPARE_DIR/comparisons.tsv" || exit 1

# Output a message about which comparison groups were found.
awk -F '\t' '
NR == 1 { next }
NF {
    counts[$3]++
    total++
}
END {
    printf "\nEach planned comparison changes one factor while holding the other two constant.\n"
    printf "  %-52s : %d\n", "Induced versus uninduced at the same time", counts["Induction"]
    printf "  %-52s : %d\n", "Different times with the same treatment", counts["Timing"]
    printf "  %-52s : %d\n", "Different cell lines at the same time and treatment", counts["CellLine"]
    printf "  %-52s : %d\n", "Total comparisons", total
}
' "$PREPARE_DIR/comparisons.tsv" || exit 1

printf '\nComparison file: %s\n\n' "$PREPARE_DIR/comparisons.tsv"








