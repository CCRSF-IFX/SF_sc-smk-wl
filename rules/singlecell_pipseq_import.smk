import re

def is_pipseeker_version_greater_than_3(pipseeker_uri):
    # Extract version using regex
    match = re.search(r"fluent-pipseeker:(\d+)\.(\d+)\.(\d+)", pipseeker_uri)
    
    if match:
        major_version = int(match.group(1))
        return major_version >= 3
    else:
        raise ValueError("Invalid Pipseeker URI format")

def get_summary_script4pipseq_data():
    if is_pipseeker_version_greater_than_3(program.pipseeker):
        return "workflow/scripts/rna/python_scripts/generateSummaryFiles4pipseq_v3.py" 
    else:
        return "workflow/scripts/rna/python_scripts/generateSummaryFiles4pipseq.py"

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

current_cellranger = program.pipseeker
 
rule count:
    output: "{sample}/barcodes/barcode_whitelist.txt"
    log: err = "run_{sample}_fluent_pipseq_count.err", log ="run_{sample}_fluent_pipseq_count.log",
    params: 
        prefix = "{sample}", 
        prefix2 = filterFastq4pipseeker, 
        cells_flag = lambda wildcards: params_cell_number[wildcards.sample],
    container:
        program.pipseeker
    shell:
        """
rm -r {params.prefix}; pipseeker full --skip-version-check --chemistry {chemistry}  --fastq {params.prefix2}. --output-path {params.prefix} --star-index-path {reference.pipseq_reference}  {params.cells_flag} 2>{log.err} 1>{log.log}
"""

rule aggregateCSV:
    input: expand(rules.count.output, sample=samples) 
    output: "AggregatedDatasets.csv"
    params: batch = "-l nodes=1:ppn=1"
    shell: "python workflow/scripts/rna/python_scripts/generateAggregateCSV.py {analysis}"

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "multiqc.smk"

rule summaryFiles:
    input: 
        expand(rules.count.output, sample=samples),
    params:
        summary_script = get_summary_script4pipseq_data()
    output: 
        "finalreport/metric_summary.xlsx",
        expand("finalreport/summaries/{sample}_report.html", sample=samples)
    shell: "python {params.summary_script}"
