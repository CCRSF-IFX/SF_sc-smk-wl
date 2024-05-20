
numcell = getattr(config, "numcells", False)
include_introns = getattr(config, "include_introns", True)

def count_expect_force():
    params_cell_number = dict()
    cells_flag = ""
    if forcecells:
        cells_flag = '--force-cells'
    else:
        ## forcecells == False & numcell has values, meaning 
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
print(params_cell_number)

current_cellranger = program.cellranger

rule fastqc4report:
    output: "fastqc4report/{sample}/{sample}_fastqc.html"
    log: err = "run_{sample}_fastqc.err", log ="run_{sample}_fastqc.log"
    params: batch = "-l nodes=1:ppn=16,mem=96gb", prefix = "fastqc4report/{sample}", prefix2 = filterFastq4nopipe 
    shell: "rm -r {params.prefix}; mkdir -p {params.prefix};zcat {params.prefix2}*R1_001.fastq.gz |{program.fastqc} -o {params.prefix} --noextract -k 5 -t 8 -f fastq stdin:{wildcards.sample} 2>{log.err} 1>{log.log}"

rule multiqc4report:
    input:
        html = expand(rules.fastqc4report.output, sample = samples)
    output: 
         stat = "multiqc4report_data/multiqc_general_stats.txt",
         html = "multiqc4report.html"
    log: err = "run_multiqc4report.err", log ="run_multiqc4report.log"
    params: batch = "-l nodes=1:ppn=16,mem=96gb", 
    shell: """
{program.multiqc} -f  -n {output.html} ./fastqc4report/  2>{log.err} 1>{log.log} 
"""

rule prep_fastq:
    input: unpack(getFirstFastqFile)
    output: R1 = "QC/Sample_{sample}/{sample}_R1.fastq.gz", R2 = "QC/Sample_{sample}/{sample}_R2.fastq.gz"
    shell: "ln -s {input.R1} {output.R1} && ln -s {input.R2} {output.R2}"

rule fastqscreen:
    input: R2 = rules.prep_fastq.output.R2
    output: two = "QC/Sample_{sample}/{sample}_R2_screen.png", two_txt = "QC/Sample_{sample}/{sample}_R2_screen.txt", two_tagged = "QC/Sample_{sample}/{sample}_R2.tagged.fastq.gz", two_tagged_filtered = "QC/Sample_{sample}/{sample}_R2.tagged_filter.fastq.gz", two_html = "QC/Sample_{sample}/{sample}_R2_screen.html"
    log: logname = "QC/Sample_{sample}/{sample}_fastq_screen.err"
    params: prefix = "QC/Sample_{sample}/"
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
    input: rules.multiqc4report.output.stat
    output: "finalreport/metric_summary.xlsx",
    params: batch = "-l nodes=1:ppn=1"
    shell: "{program.python2_7} {program.pythonscripts}/generateSummaryFiles4nopipe.py"
