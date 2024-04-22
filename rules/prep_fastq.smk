rule prep_fastq:
    input: unpack(getFirstFastqFile)
    output: R1 = "QC/Sample_{sample}/{sample}_R1.fastq.gz", R2 = "QC/Sample_{sample}/{sample}_R2.fastq.gz"
    shell: "ln -s {input.R1} {output.R1} && ln -s {input.R2} {output.R2}"
