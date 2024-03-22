
numcell = getattr(config, "numcells", False)
include_introns = getattr(config, "include_introns", True)

def count_expect_force():
    params_cell_number = dict()
    cells_flag = ""
    if forcecells:
        cells_flag = '--force-cells'
    else:
        ## that --expect is used
        if numcell:
            cells_flag = '--expect-cells'
    ## cell number to be auto-estimated
    if cells_flag == "":
        params_cell_num = dict(zip(samples, [ " " for sample in samples ]))
    else:
        numcells = numcell.split(',') 
        params_cell_num = dict(zip(samples, [ f"{cells_flag}={num}" for num in numcells ]))
    return params_cell_num


params_cell_number = count_expect_force()
sflog.info(params_cell_number)

current_cellranger = program.cellranger
 
rule count:
    output: "{sample}/outs/web_summary.html"
    log: err = "run_{sample}_10x_cellranger_count.err", log ="run_{sample}_10x_cellranger_count.log"
    params: batch = "-l nodes=1:ppn=16,mem=96gb", prefix = "{sample}", prefix2 = filterFastq, cells_flag=lambda wildcards: params_cell_number[wildcards.sample], include_introns = str(include_introns).lower()
    shell: "rm -r {params.prefix}; module load bcl2fastq2; {cmd_cellranger} count --include-introns {params.include_introns} --id={params.prefix} --sample={params.prefix} --fastqs={params.prefix2} {params.cells_flag} --transcriptome={reference.transcriptome} 2>{log.err} 1>{log.log}"

rule aggregateCSV:
    input: expand("{sample}/outs/web_summary.html", sample=samples)
    output: "AggregatedDatasets.csv"
    params: batch = "-l nodes=1:ppn=1"
    shell: "python generateAggregateCSV.py {analysis}"

rule aggregate:
    input: csv="AggregatedDatasets.csv"
    output: touch("aggregate.complete")
    log: err="run_10x_aggregate.err", log="run_10x_aggregate.log"
    params: batch = "-l nodes=1:ppn=16,mem=96gb"
    shell: "cellranger aggr --id=AggregatedDatasets --csv={input.csv} --normalize=mapped 2>{log.err} 1>{log.log}"

rule prep_fastq:
    input: unpack(getFirstFastqFile)
    output: R1 = "QC/Sample_{sample}/{sample}_R1.fastq.gz", R2 = "QC/Sample_{sample}/{sample}_R2.fastq.gz"
    shell: "ln -s {input.R1} {output.R1} && ln -s {input.R2} {output.R2}"

rule fastqscreen:
    input: R2 = rules.prep_fastq.output.R2
    output: two = "QC/Sample_{sample}/{sample}_R2_screen.png", two_txt = "QC/Sample_{sample}/{sample}_R2_screen.txt", two_tagged = "QC/Sample_{sample}/{sample}_R2.tagged.fastq.gz", two_tagged_filtered = "QC/Sample_{sample}/{sample}_R2.tagged_filter.fastq.gz", two_html = "QC/Sample_{sample}/{sample}_R2_screen.html"
    log: logname = "QC/Sample_{sample}/{sample}_fastq_screen.err"
    params: prefix = "QC/Sample_{sample}/"
    conda: "../envs/fastqscreen.yaml" 
    shell: "{program.fastq_screen} --outdir {params.prefix} --threads {clusterConfig[fastqscreen][threads]} --subset 5000000 --nohits --conf {program.conf} --aligner bowtie2 {input.R2} 2>{log.logname}"

rule kraken:
    input: R2 = rules.prep_fastq.output.R2
    output: result = "QC/Sample_{sample}/{sample}.kraken", krona = "QC/Sample_{sample}/{sample}.kraken.krona", report = "QC/Sample_{sample}/{sample}.kraken.report.txt"
    log: err = "QC/Sample_{sample}/kraken.err"
    params: prefix = "QC/Sample_{sample}/{sample}.kraken"
    shell: "{program.kraken2} --gzip-compressed --threads {clusterConfig[kraken][threads]} --db {program.kraken2db} --output {params.prefix} --report {output.report} {input.R2} 2> {log.err} && cut -f2,3 {output.result} > {output.krona}"

rule multiqc:
    input: expand("QC/Sample_{sample}/{sample}_R2_screen.png", sample=samples), expand("QC/Sample_{sample}/{sample}.kraken.report.txt", sample=samples)
    output: "QC/" + project_name + "_multiqc.html"
    params: batch = "-l nodes=1:ppn=4,mem=16g"
    shell: "{program.multiqc} -f -c {program.multiqc_conf} -n {output} ./QC"

rule copyScripts:
    output: directory("scripts")
    params: batch = "-l nodes=1:ppn=1"
    shell: "mkdir scripts; cp -r {program.rscripts} scripts; cp -r {program.pythonscripts} scripts"

rule summaryFiles:
    input: expand("{sample}/outs/web_summary.html", sample=samples)
    output: "finalreport/metric_summary.xlsx", expand("finalreport/summaries/{sample}_web_summary.html", sample=samples)
    params: batch = "-l nodes=1:ppn=1"
    shell: "{program.python2_7} {program.pythonscripts}/generateSummaryFiles.py"
