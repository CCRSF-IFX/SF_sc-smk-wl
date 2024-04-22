rule kraken:
    input: R2 = rules.prep_fastq.output.R2
    output: result = "QC/Sample_{sample}/{sample}.kraken", krona = "QC/Sample_{sample}/{sample}.kraken.krona", report = "QC/Sample_{sample}/{sample}.kraken.report.txt"
    log: err = "QC/Sample_{sample}/kraken.err"
    params: prefix = "QC/Sample_{sample}/{sample}.kraken"
    shell: "{program.kraken2} --gzip-compressed --threads {clusterConfig[kraken][threads]} --db {program.kraken2db} --output {params.prefix} --report {output.report} {input.R2} 2> {log.err} && cut -f2,3 {output.result} > {output.krona}"
