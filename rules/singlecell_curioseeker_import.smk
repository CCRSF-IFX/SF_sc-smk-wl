import os

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

current_cellranger = "curioseeker:3.0.0"

rule generate_sample_csv:
    input:
        csv = config.libraries
    output:
        meta4curio = "metadata_for_curiodata.csv",
        csv_list = expand("{sample}.csv", sample = samples)
    shell:
        """
python workflow/scripts/curio/generate_sample_sheet.py {input.csv} {analysis} > {output.meta4curio} 
touch {output.csv_list}
"""


"""
`filterFastq4pipseeker` is used here to return the path to sample folder
with the symbolic links to the fastq files. 
First concatenate the R1 and R2 files, seperately, and then run curioseeker
"""
rule count:
    input: 
        meta4curio =  rules.generate_sample_csv.output.meta4curio
    output: 
        "{sample}/OUTPUT/{sample}/{sample}_Metrics.csv"
    log: 
        err = os.path.join(analysis, "run_{sample}_curioseeker.err"), 
        log = os.path.join(analysis, "run_{sample}_curioseeker.log"),
    params: 
        prefix = "{sample}", 
        prefix2 = filterFastq4nopipe,
        sample_sheet = os.path.join(analysis, "{sample}.csv"),
    shell: 
        """
rm -rf {analysis}/{params.prefix} {analysis}/{params.prefix}_work
mkdir -p fastq4curio/ {params.prefix}
# exit 1 is a shell command that terminates the script with an error status (1) 
cat {params.prefix2}*_R1_001.fastq.gz > fastq4curio/{params.prefix}_R1.fastq.gz || exit 1
cat {params.prefix2}*_R2_001.fastq.gz > fastq4curio/{params.prefix}_R2.fastq.gz || exit 1 
cd {params.prefix}
module load nextflow/23.08.0
nextflow run /mnt/ccrsf-ifx/Software/tools/curioseeker/curioseeker-v3.0.0/main.nf \
                --input {params.sample_sheet} \
                --outdir {analysis}/{params.prefix} \
                -work-dir {analysis}/{params.prefix}_work \
                --igenomes_base /mnt/ccrsf-ifx/RefGenomes/SingleCell_REF/Curioseeker \
                -resume \
                -profile slurm \
                -config /mnt/ccrsf-ifx/Software/tools/curioseeker/curioseeker-v3.0.0/curioseeker_slurm.config \
                2>{log.err} 1>{log.log}
"""

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "fastqc4QC.smk"
include: "multiqc.smk"

rule copyScripts:
    output: directory("scripts")
    params: batch = "-l nodes=1:ppn=1"
    shell: "mkdir scripts; cp -r {program.rscripts} scripts; cp -r {program.pythonscripts} scripts"

rule summaryFiles:
    input: expand(rules.count.output, sample = samples) 
    output: "finalreport/metric_summary.xlsx", 
    container: program.global_container
    shell: 
        """
python workflow/scripts/curio/generateSummaryFiles.py
"""
