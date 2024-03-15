#!/bin/bash

# Function to simulate a sample 
create_sample() {
  sample_name=$1  # Name of the sample (e.g., "1_ctrl") 
  base_dir="$2/$1"

  # Create the base directory
  mkdir -p "$base_dir" 

  # Create fake FASTQ files in the structure
  for lane in L001 L002; do
    for read in R1 R2 I1 I2; do
        for index in 001; do
            file="${sample_name}_S${sample_name:0:1}_${lane}_${read}_${index}.fastq.gz"
            touch "$base_dir/$file" 
        done
    done
  done
}

if [ -z "$1" ] && [ -z "$2" ]; then
    echo -e "\nError: Please provide a base directory path and sample name \n"
    exit 1
fi

# Get sample name from second argument (assuming first argument is basedir)
sample_name=$2

create_sample "$sample_name" "$1"

