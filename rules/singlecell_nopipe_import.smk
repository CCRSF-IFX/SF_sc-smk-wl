
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
#print(params_cell_number)

current_cellranger = program.cellranger

rule fastqc4report:
    output: 
        "fastqc4report/{sample}/{sample}_fastqc.html"
    log: 
        err = "run_{sample}_fastqc.err", 
        log ="run_{sample}_fastqc.log"
    params: 
        prefix = "fastqc4report/{sample}", 
        prefix2 = filterFastq4nopipe 
    container: program.fastqc
    shell: 
        """
rm -r {params.prefix}; mkdir -p {params.prefix};zcat {params.prefix2}*R1_001.fastq.gz |fastqc -o {params.prefix} --noextract -k 5 -t 8 -f fastq stdin:{wildcards.sample} 2>{log.err} 1>{log.log}
"""

rule multiqc4report:
    input:
        html = expand(rules.fastqc4report.output, sample = samples)
    output: 
         stat = "multiqc4report_data/multiqc_general_stats.txt",
         html = "multiqc4report.html"
    log: err = "run_multiqc4report.err", log ="run_multiqc4report.log"
    container:
        program.multiqc
    shell:
        """
multiqc -f  -n {output.html} ./fastqc4report/  2>{log.err} 1>{log.log} 
"""

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "multiqc.smk"

rule summaryFiles:
    input: 
        rules.multiqc4report.output.stat
    output: 
        "finalreport/metric_summary.xlsx",
    shell: 
        """
python {analysis}/workflow/scripts/nopipe/generateSummaryFiles4nopipe.py
"""
