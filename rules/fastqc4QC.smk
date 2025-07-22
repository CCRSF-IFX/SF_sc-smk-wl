import os
print(record_fastqpath)
rule fastqc4QC:
    # input: 
    #     log_prep_fq = rules.prep_fastq_folder.output.log
    params: 
        fastqc_outdir = os.path.join(analysis, "QC/Sample_{sample}/fastqc_outdir"),
        fastq_path = lambda wildcards: record_fastqpath[wildcards.sample],
    output: 
        os.path.join(analysis, "QC/Sample_{sample}/fastqc_outdir/fastqc.log")
    container:
        program.fastqc
    threads: 8 # could be overwritten by the config file
    shell:
        """
        mkdir -p {params.fastqc_outdir}
        echo "Running FastQC in the folder of: "
        pwd
        fastq_files=$(find {params.fastq_path} -name "*.fastq.gz" -or -name "*.fq.gz")
        fastqc --threads {threads} \
               --outdir={params.fastqc_outdir} \
               $fastq_files >{output} 2>&1
        """
