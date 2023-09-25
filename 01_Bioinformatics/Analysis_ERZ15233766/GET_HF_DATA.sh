#!/bin/bash
module load sra-tools/3.0.3 gnu-parallel/20161122
# Set the path to the directory where you want to store the downloaded files
DOWNLOAD_DIR="."
# Set the SRR accession numbers as an array
SRR=$(cat ${1})

# Read the SRR accession numbers from a file into an array
mapfile -t $SRR

# Define function to download and convert a single accession number
download_srr() {
    echo "starts downloading "$S1""
    prefetch "$1"
    fastq-dump --gzip --split-files --outdir "$DOWNLOAD_DIR" "$1"
}

# Export the function so it can be used by GNU Parallel
export -f download_srr

# Download and convert all SRR accession numbers in parallel using GNU Parallel
parallel -j $SLURM_CPUS_PER_TASK download_srr ::: "${SRR[@]}"