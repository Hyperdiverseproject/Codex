Snakefile README

Overview
This Snakefile is designed for processing sequencing data and analyzing fragments. It performs the following steps:
Trimming and merging reads using fastp.
Demultiplexing reads based on barcodes file.
Preparing a tag database using BLAST.
Generating fasta files of the reads, running BLAST, and segregating fragments.
Assembling fragments using Trinity.
Calculating read depth using minimap.
Performing BLAST on NT (nucleotide) database.
Removing tags and filtering sequences.
Merging fragments and generating the final output.

Configuration
The configuration is stored in a YAML file named config.yaml.
Input data includes read pairs file, barcodes file, and the NCBI NT database (optional).
Usage

To run the Snakefile, use the following command:
snakemake -s Snakefile
Ensure that Snakemake and the required tools (fastp, Trinity, BLAST, minimap2) are installed and accessible in your environment.

Directory Structure
The Snakefile assumes a specific directory structure:

Input reads are organized in the specified format (Library_*_R1.fastq.gz, Library_*_R2.fastq.gz); required
Barcodes are provided in a text file (barcodes.txt), one line per well with the well name, the fragment 1 forward tag, then the fragment 1 reverse tag then the fragment 2 forward tag, then the fragment 2 reverse tag separated by a tabulation; required.
Output is organized into different steps (First_step, Second_step, ..., Final) within specified output directories.

Output
The final sequences are generated in the Final directory as {sample}_final.sequences.fa.
Various intermediate files are created in step-specific output directories.

Dependencies
fastp
Trinity
BLAST
minimap2
Ensure that these tools are installed and available in your PATH.

Additional Notes
Adjust the configuration parameters in config.yaml based on your dataset and analysis requirements.
Review the Snakefile rules and shell commands for specific details on each step.
