if forcecells:
    numcell = config.numcells
    numcells = config.numcells.split(',')
    dict2 = dict(zip(samples,numcells))

def count_force_numbers(wildcards):
    if forcecells:
        return('--force-cells=%s' % dict2[wildcards.sample])
    else:
        return('')

current_cellranger = program.atac_cellranger

rule atac_count:
    output: 
        "{sample}/outs/web_summary.html", 
        "{sample}/outs/filtered_peak_bc_matrix.h5", 
        "{sample}/outs/singlecell.csv", "{sample}/outs/fragments.tsv.gz"
    log: 
        err = "run_{sample}_10x_cellranger_atac.err", 
        log ="run_{sample}_10x_cellranger_atac.log"
    params: 
        prefix = "{sample}", 
        prefix2 = filterFastq, 
        cells_flag = count_force_numbers
    container: 
        program.cellranger_atac
    shell:
        """
rm -r {params.prefix}; cellranger-atac count --id={params.prefix} --fastqs={params.prefix2} --reference={reference.atac_reference} {params.cells_flag} 2>{log.err} 1>{log.log}
"""

rule aggregateCSV:
    input: 
        expand("{sample}/outs/web_summary.html", sample=samples)
    output: 
        "AggregatedDatasets.csv"
    shell: 
        """
python {analysis}/workflow/scripts/atac/generateAggregateCSV.py {analysis}
"""

rule atac_aggregate:
    input: 
        csv="AggregatedDatasets.csv"
    output: 
        touch("aggregate.complete")
    log: 
        err="run_10x_aggregate.err", 
        log="run_10x_aggregate.log"
    container: program.atac_cellranger
    shell: 
        """
cellranger-atac aggr --id=AggregatedDatasets --csv={input.csv} --reference={reference.atac_reference} --normalize=depth 2>{log.err} 1>{log.log}
"""

rule prep_fastq:
    input: unpack(getFirstFastqFile)
    output: 
        R1 = "QC/Sample_{sample}/{sample}_R1.fastq.gz", 
        R3 = "QC/Sample_{sample}/{sample}_R3.fastq.gz"
    shell: "ln -s {input.R1} {output.R1} && ln -s {input.R3} {output.R3}"

rule fastqscreen:
    input: 
        R1 = rules.prep_fastq.output.R1, 
        R3 = rules.prep_fastq.output.R3
    output: 
        one = "QC/Sample_{sample}/{sample}_R1_screen.png", 
        two = "QC/Sample_{sample}/{sample}_R3_screen.png", 
        one_txt = "QC/Sample_{sample}/{sample}_R1_screen.txt", 
        two_txt = "QC/Sample_{sample}/{sample}_R3_screen.txt", 
        one_html = "QC/Sample_{sample}/{sample}_R1_screen.html", 
        two_html = "QC/Sample_{sample}/{sample}_R3_screen.html"
    log: 
        logname = "QC/Sample_{sample}/{sample}_fastq_screen.err"
    container: 
        program.fastq_screen
    params: 
        prefix = "QC/Sample_{sample}/"
    shell: 
        """
fastq_screen --outdir {params.prefix} --threads 4 --subset 5000000 --nohits --conf {program.conf} --aligner bowtie2 {input.R1} {input.R3} 2>{log.logname}
"""

rule kraken:
    input: 
        R1 = rules.prep_fastq.output.R1, 
        R3 = rules.prep_fastq.output.R3
    output: 
        result = "QC/Sample_{sample}/{sample}.kraken", 
        krona = "QC/Sample_{sample}/{sample}.kraken.krona", 
        report = "QC/Sample_{sample}/{sample}.kraken.report.txt"
    log: 
        err = "QC/Sample_{sample}/kraken.err"
    params: 
        prefix = "QC/Sample_{sample}/{sample}.kraken"
    container: program.kraken2
    shell: 
        """
kraken2 --gzip-compressed --threads 4 --db {program.kraken2db} --output {params.prefix} --report {output.report} {input.R1} {input.R3} 2> {log.err} && cut -f2,3 {output.result} > {output.krona}
"""

rule multiqc_qc:
    input: 
        expand("QC/Sample_{sample}/{sample}_R3_screen.png", sample=samples), 
        expand("QC/Sample_{sample}/{sample}.kraken.report.txt", sample=samples)
    output: 
        "QC/" + project_name + "_multiqc.html"
    container: program.multiqc
    shell: 
        """
multiqc -f -c {program.multiqc_conf} -n {output} ./QC
"""

rule summaryFiles:
    input: 
        expand("{sample}/outs/web_summary.html", sample=samples)
    output: 
        "finalreport/metric_summary.xlsx", 
        expand("finalreport/summaries/{sample}_web_summary.html", sample=samples)
    shell: 
        """
python {analysis}/workflow/scripts/atac/generateSummaryFiles.py
"""
