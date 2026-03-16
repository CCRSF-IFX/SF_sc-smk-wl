
samples = fastq_samples
include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "prep_fastq_folder_ln.smk"
include: "fastqc4QC.smk"
include: "multiqc.smk"
samples = multi_samples
