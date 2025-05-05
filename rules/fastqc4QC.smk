
rule fastqc4QC:
    params: 
        prefix = "{sample}", 
        prefix2 = filterFastq4nopipe, 
        fastqc_outdir = "QC/Sample_{sample}/fastqc_outdir"
    output: 
        "QC/Sample_{sample}/fastqc_outdir/fastqc.log"
    container:
        program.fastqc
    threads: 8 # could be overwritten by the config file
    shell:
        """
        # find all the fastq files in the sample directory and perform fastqc on them
        fastq_files=$(find {params.prefix} -name "*.fastq.gz")
        fastqc --threads {threads} --outdir={params.fastqc_outdir} $fastq_files >{output} 2>&1
        """