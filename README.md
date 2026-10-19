# MyFirstPipeline

MyFirstPipeline is a Bash and AWK workflow for analysing paired-end RNA-seq data from a Trypanosoma congolense IL3000 RNAi experiment. It checks sequencing quality, aligns reads to the reference genome, counts alignments overlapping annotated genes, and calculates group mean counts and fold changes.

The files to run the pipeline are contained within an encrypted tar archive:

https://github.com/niktaghanei/MyFirstPipeline/blob/main/Unsupervised_Learners.MyFirstPipeline.tar.gz.cpt 

The encrypted archive contains the pipeline code and the project history file, all_the_things_I_did. Sequencing data and reference files are not included and must be accessible on the server.


## User Guide

### Requirements

Run on Linux with Git, Bash, AWK, ccrypt, tar, unzip, FastQC, Bowtie2, SAMtools and BEDtools available.

### 1 Download the submission

git clone https://github.com/niktaghanei/MyFirstPipeline.git

cd MyFirstPipeline

### 2 Decrypt and extract the pipeline

ccrypt -c Unsupervised_Learners.MyFirstPipeline.tar.gz.cpt \

    > Unsupervised_Learners.MyFirstPipeline.tar.gz

mkdir -p pipeline

tar -xzf Unsupervised_Learners.MyFirstPipeline.tar.gz  -C pipeline

cd pipeline

### 3 Configure the pipeline

Create a local config for the pipeline run to use, from the example:

cp examples/config.example.txt config.txt

vi config.txt

Update these settings:

- METADATA_FILE = Full path to the sample metadata file

- FASTQ_DIR = Full path to the paired FASTQ folder

- GENOME_FASTA = Full path to the reference genome FASTA

- GENE_BED = Full path to the gene annotation BED

- WORK_DIR = Folder for generated results

- THREADS = Number (which will depend on computing resources)

- MODE = "interactive" or "quiet"

Use absolute paths for the inputs. The metadata must follow the tab-delimited format shown in examples/metadata_template.tsv.

Interactive mode displays progress and pauses for input. Quiet mode runs without prompts and saves summaries to the report.

Leave the QC thresholds at their supplied values unless intentionally changing the assessment policies.

To save and exit vi, press Esc, type :wq, then press Enter.

### 4 Run the pipeline

Run from the extracted folder containing pipeline.sh and config.txt:

bash pipeline.sh

### 5 Review the results

The pipeline prints the report and results locations. Outputs include:

- Sample groups and planned comparisons.

- FastQC reports, metrics and sample assessments.

- Sorted and indexed alignments.

- Gene counts with descriptions.

- Group mean counts and ranked fold changes.




