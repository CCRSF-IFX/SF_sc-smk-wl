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

current_cellranger = program.cellranger
 
rule count:
    output: "{sample}/outs/web_summary.html"
    log: err = "run_{sample}_10x_cellranger_count.err", log ="run_{sample}_10x_cellranger_count.log"
    params: 
        prefix = "{sample}", 
        prefix2 = filterFastq, 
        cells_flag=lambda wildcards: params_cell_number[wildcards.sample], 
        include_introns = str(include_introns).lower(),
        reference = get_reference_transcriptome, 
        option_flag = get_optional_flag_from_lib_csv 
    container: program.cellranger
    shell: 
        """
        rm -r {params.prefix}; 
        cellranger count {flag4cellranger_create_bam} --include-introns {params.include_introns} \
              --id={params.prefix} --sample={params.prefix} \
              --fastqs={params.prefix2} {params.cells_flag} \
              --transcriptome={params.reference} \
              --feature-ref={config.features} \
              {params.option_flag} 2>{log.err} 1>{log.log}"

rule aggregateCSV:
    input: expand("{sample}/outs/web_summary.html", sample=samples)
    output: "AggregatedDatasets.csv"
    shell: "python workflow/scripts/rna/python_scripts/generateAggregateCSV.py {analysis}"

rule aggregate:
    input: csv="AggregatedDatasets.csv"
    output: touch("aggregate.complete")
    log: err="run_10x_aggregate.err", log="run_10x_aggregate.log"
    container: program.cellranger 
    shell: "cellranger aggr --id=AggregatedDatasets --csv={input.csv} --normalize=mapped 2>{log.err} 1>{log.log}"

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "fastqc4QC.smk"
include: "multiqc.smk"
#include: "picard_rnaseqmetrics.smk"

rule copyScripts:
    output: directory("scripts")
    params: batch = "-l nodes=1:ppn=1"
    shell: "mkdir scripts; cp -r {program.rscripts} scripts; cp -r {program.pythonscripts} scripts"

rule summaryFiles:
    input: expand("{sample}/outs/web_summary.html", sample=samples)
    output: "finalreport/metric_summary.xlsx", expand("finalreport/summaries/{sample}_web_summary.html", sample=samples)
    params: batch = "-l nodes=1:ppn=1"
    shell: "python workflow/scripts/rna/python_scripts/generateSummaryFiles.py"

