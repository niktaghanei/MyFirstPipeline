#!/bin/bash

METADATA_FILE="$1"
WORK_DIR="$2"
FASTQ_DIR="$3"
PREPARE_DIR="$WORK_DIR/prepare_inputs"

mkdir -p "$PREPARE_DIR" || exit 1

# Expected filenames come from the metadata.
awk -F '\t' '
#set output field separator to new line because the file names are in 2 columns but we want a list in 1 column to make looping easier later
BEGIN { OFS = "\n" }
NR > 1 && NF { print $6, $7 }
' "$METADATA_FILE" > "$PREPARE_DIR/expected_files.txt" || exit 1

# Actual filenames come from the FASTQ folder - they all end ($) with .fq.gz
# \ are needed to read the . as literal full stops
# exit as an error (1) if the files are not able to be read or written
ls "$FASTQ_DIR" |
awk '/\.fq\.gz$/' > "$PREPARE_DIR/actual_files.txt" || exit 1

awk '
#check we are in the first declared file - store the lines in an array as they contain the file names provided in the fastq folder
FILENAME == ARGV[1] {
    files[$0] = 1
    # Then skip to the next line - this will mean this section is only run on the lines of the first file and the next section is run only on the lines of the second file
    next
}
{
    total_count++

    #if the expected file is not in the list of files in the folder then increment a missing counter and store the file name in an array for later reporting
    if (!($0 in files)) {
        missing_count++
        missing_files[missing_count] = $0
    }
}
END {
    printf "\nChecked %d expected FASTQ filenames from the metadata against the FASTQ folder.\n", total_count

    if (total_count == 0) {
        print "No expected filenames were available to check."
        exit 1
    } else if (missing_count == 0) {
        print "All expected filenames were found."
    } else {
        printf "%d expected filenames were missing:\n", missing_count

        for (i = 1; i <= missing_count; i++) {
            print "  " missing_files[i]
        }
        exit 1
    }
}
' "$PREPARE_DIR/actual_files.txt" "$PREPARE_DIR/expected_files.txt" || exit 1

printf '\n'
