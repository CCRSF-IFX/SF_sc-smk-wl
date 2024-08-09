
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
rule generate_sample_csv:
    input:
        csv = config.metadata
    output:
        meta4curio = "metadata_for_curiodata.csv",
        csv_list = expand("{sample}.csv", sample = samples)
    shell:
        """
python workflow/scripts/curio/generate_sample_sheet.py > {output.meta4curio} 
touch {output.csv_list}
"""


"""
`filterFastq4pipseeker` is used here to return the path to sample folder
with the symbolic links to the fastq files. 
First concatenate the R1 and R2 files, seperately, and then run curioseeker
"""
rule count:
    input: 
        sample_sheet = "{sample}.csv",
        meta4curio =  rules.generate_sample_csv.output.meta4curio
    output: 
        "{sample}/OUTPUT/{sample}/{sample}_Metrics.csv"
    log: 
        err = "run_{sample}_curioseeker.err", 
        log ="run_{sample}_curioseeker.log"
    params: prefix = "{sample}", prefix2 = filterFastq4pipseeker
    shell: 
        """
rm -r {params.prefix};
mkdir -p fastq4curio/ 
cat {params.prefix2}*_R1_001.fastq.gz > fastq4curio/{params.prefix}_R1.fastq.gz &
cat {params.prefix2}*_R2_001.fastq.gz > fastq4curio/{params.prefix}_R2.fastq.gz &
wait
module load nextflow singularity
nextflow run /mnt/ccrsf-ifx/Software/tools/curioseeker/curioseeker-v3.0.0/main.nf \
                --input {input.sample_sheet} \
                --outdir {params.prefix} \
                -work-dir {params.prefix}_work \
                --igenomes_base /mnt/ccrsf-ifx/Software/tools/curioseeker/References \
                -resume \
                -profile slurm \
                -config /mnt/ccrsf-ifx/Software/tools/curioseeker/curioseeker-v3.0.0/curioseeker_slurm.config
"""

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "multiqc.smk"

rule copyScripts:
    output: directory("scripts")
    params: batch = "-l nodes=1:ppn=1"
    shell: "mkdir scripts; cp -r {program.rscripts} scripts; cp -r {program.pythonscripts} scripts"

rule summaryFiles:
    input: expand(rules.count.output, sample = samples) 
    output: "finalreport/metric_summary.xlsx", 
    shell: "python workflow/scripts/curio/generateSummaryFiles.py"

