"""
This script is developped to generate sample sheet for curioseeker. 

Shaojun Xie on Aug 9, 2024 
"""

import sys
import os 

args = sys.argv

if len(args) != 3:
    print("Usage: script_name input_csv analysis_folder")
    sys.exit(1)

analysis_folder = args[2]

with open(args[1], "r") as fcsv:
    header = next(fcsv)
    header = header.rstrip()
    print(f"{header},fastq_1,fastq2")
    for line in fcsv:
        ele = line.rstrip().split(",")
        sample = ele[0]
        with open(f"{sample}.csv", "w+") as fout:
            fq1 = os.path.join(analysis_folder, f"fastq4curio/{sample}_R1.fastq.gz")
            fq2 = os.path.join(analysis_folder, f"fastq4curio/{sample}_R2.fastq.gz")
            fout.write(f"{header},fastq_1,fastq2\n")
            fout.write(f"{line.rstrip()},{fq1},{fq2}\n")
        print(f"{header},{fq1},{fq2}")
