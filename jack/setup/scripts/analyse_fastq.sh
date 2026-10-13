#!/bin/bash

METADATA_FILE="$1"
FASTQC_EXTRACT_DIR="$2"

awk -F '\t' '
BEGIN {
    OFS = "\t"

    print "SampleName",
          "R1TotalReads", "R2TotalReads",
          "R1ReadLength", "R2ReadLength",
          "R1GCPercent", "R2GCPercent",
          "R1MinimumMedianQuality", "R2MinimumMedianQuality",
          "R1MinimumLowerQuartileQuality", "R2MinimumLowerQuartileQuality",
          "R1ModeSequenceQuality", "R2ModeSequenceQuality",
          "R1MaximumNPercentage", "R2MaximumNPercentage"
}

# First file: store metrics using filename and metric name as the key.
FILENAME == ARGV[1] {
    if (FNR > 1 && NF) {
        available[$1] = 1
        metrics[$1, "TotalReads"] = $2
        metrics[$1, "ReadLength"] = $3
        metrics[$1, "GCPercent"] = $4
        metrics[$1, "MinimumMedianQuality"] = $5
        metrics[$1, "MinimumLowerQuartileQuality"] = $6
        metrics[$1, "ModeSequenceQuality"] = $7
        metrics[$1, "MaximumNPercentage"] = $8
    }
    next
}

# Second file: read sample names and mate filenames.
FNR == 1 { next }
NF == 0 { next }

{
    print $1,
          metrics[$6, "TotalReads"], metrics[$7, "TotalReads"],
          metrics[$6, "ReadLength"], metrics[$7, "ReadLength"],
          metrics[$6, "GCPercent"], metrics[$7, "GCPercent"],
          metrics[$6, "MinimumMedianQuality"], metrics[$7, "MinimumMedianQuality"],
          metrics[$6, "MinimumLowerQuartileQuality"], metrics[$7, "MinimumLowerQuartileQuality"],
          metrics[$6, "ModeSequenceQuality"], metrics[$7, "ModeSequenceQuality"],
          metrics[$6, "MaximumNPercentage"], metrics[$7, "MaximumNPercentage"]
}

' "$FASTQC_EXTRACT_DIR/combined_metrics.txt" "$METADATA_FILE" \
    > "$FASTQC_EXTRACT_DIR/paired_metrics.tsv"

printf 'Paired metrics saved at: %s\n' "$FASTQC_EXTRACT_DIR/paired_metrics.tsv"
