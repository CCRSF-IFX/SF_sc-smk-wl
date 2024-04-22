rule fastqscreen:
    input: R2 = rules.prep_fastq.output.R2
    output: two = "QC/Sample_{sample}/{sample}_R2_screen.png", two_txt = "QC/Sample_{sample}/{sample}_R2_screen.txt", two_tagged = "QC/Sample_{sample}/{sample}_R2.tagged.fastq.gz", two_tagged_filtered = "QC/Sample_{sample}/{sample}_R2.tagged_filter.fastq.gz", two_html = "QC/Sample_{sample}/{sample}_R2_screen.html"
    log: logname = "QC/Sample_{sample}/{sample}_fastq_screen.err"
    params: prefix = "QC/Sample_{sample}/"
    container: program.fastq_screen
    shell: "fastq_screen --outdir {params.prefix} --threads {clusterConfig[fastqscreen][threads]} --subset 5000000 --nohits --conf {program.conf} --aligner bowtie2 {input.R2} 2>{log.logname}"
