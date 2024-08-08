numcell = getattr(config, "numcells", False)
def count_expect_force():
    params_cell_number = dict()
    cells_flag = ""
    if forcecells:
        cells_flag = '--force-cells'
        numcells = numcell.split(',') 
        params_cell_num = dict(zip(samples, [ f"{cells_flag} {num}" for num in numcells ]))
    else:
        params_cell_num = dict(zip(samples, [ f"" for sample in samples ]))
    return params_cell_num

params_cell_number = count_expect_force()
print(params_cell_number)

# this variable is used in `bin/currentsnake/single_cell/Snakefile_singlecell_rules`
# so I keep the old name. 
current_cellranger = program.pipseeker

rule count:
    output: "{sample}/barcodes/barcode_whitelist.txt"
    log: err = "run_{sample}_fluent_pipseq_count.err", log ="run_{sample}_fluent_pipseq_count.log",
    params: batch = "-l nodes=1:ppn=16,mem=96gb", prefix = "{sample}", prefix2 = filterFastq4pipseeker, cells_flag = lambda wildcards: params_cell_number[wildcards.sample],
    container: program.pipseeker
    shell: 
        """
rm -r {params.prefix}; pipseeker full --skip-version-check  --fastq {params.prefix2}. --output-path {params.prefix} --star-index-path {reference.pipseq_reference}  {params.cells_flag} 2>{log.err} 1>{log.log}
"""

rule aggregateCSV:
    input: expand(rules.count.output, sample=samples)
    output: "AggregatedDatasets.csv"
    params: batch = "-l nodes=1:ppn=1"
    shell: "{program.python2_7} {program.pythonscripts}/generateAggregateCSV.py {analysis}"

rule prep_fastq:
    input: unpack(getFirstFastqFile)
    output: R1 = "QC/Sample_{sample}/{sample}_R1.fastq.gz", R2 = "QC/Sample_{sample}/{sample}_R2.fastq.gz"
    shell: "ln -s {input.R1} {output.R1} && ln -s {input.R2} {output.R2}"

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "multiqc.smk"

print("XXX")
rule summaryFiles:
    input: expand(rules.count.output, sample=samples)
    output: "finalreport/metric_summary.xlsx"
    shell: "python workflow/scripts/rna/python_scripts/generateSummaryFiles4pipseq.py"

