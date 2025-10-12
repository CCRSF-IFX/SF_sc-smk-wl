import re

current_cellranger = program.pipseeker

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "prep_fastq_folder_ln.smk"
include: "fastqc4QC.smk"
include: "multiqc.smk"
include: "prep_fastq_folder_concat.smk"

if getattr(config, "chemistry", False) is False:
    raise ValueError("chemistry is a required config parameter for parsebio")
if getattr(config, "kit", False) is False:
    raise ValueError("kit is a required config parameter for parsebio")
if getattr(reference, "transcriptome_parsebio", False) is False:
    raise ValueError("transcriptome_parsebio is a required config parameter for parsebio")
parsebio_sample_list = getattr(config, "parsebio_sample_list", "")
if parsebio_sample_list == "":
    raise ValueError("parsebio_sample_list is a required config parameter for parsebio")

rule split_pipe_all:
    input: 
        r1 = rules.make_fastq_concat.output.r1,
        r2 = rules.make_fastq_concat.output.r2,
    output: 
        summary_html = "{sample}/all-sample_analysis_summary.html"
    log: 
        err = "run_{sample}_parsebio_split_pipe_all.err", 
        log ="run_{sample}_parsebio_split_pipe_all.log",
    params:
        prefix = "{sample}",
        chemistry = config.chemistry,
        kit = config.kit,
    shell:
        """
rm -r {params.prefix}; 
{program.parsebio_split_pipe} \
    --mode all \
    --chemistry {params.chemistry} \
    --kit {params.kit} \
    --genome_dir {reference.transcriptome_parsebio} \
    --fq1 {input.r1} \
    --fq2 {input.r2} \
    --output_dir {params.prefix} \
    --samp_list {parsebio_sample_list} \
    2>{log.err} 1>{log.log}
"""

rule split_pipe_comb:
    input: 
        expand(rules.split_pipe_all.output, sample=samples)
    params:
        outdir = "split_pipe_comb",
        chemistry = config.chemistry,
        kit = config.kit,
        sublibs = expand(os.path.join(analysis, "{sample}"), sample=samples)
    output: 
        summary_html = "split_pipe_comb/all-sample_analysis_summary.html"
    log: 
        err = "run_parsebio_split_pipe_comb.err", 
        log ="run_parsebio_split_pipe_comb.log",

    shell:
        """
{program.parsebio_split_pipe} \
    --mode comb \
    --chemistry {config.chemistry} \
    --kit {config.kit} \
    --genome_dir {reference.transcriptome_parsebio} \
    --sublibraries {params.sublibs} \
    --output_dir {params.outdir} \
    --samp_list {parsebio_sample_list} \
    2>{log.err} 1>{log.log}
"""

rule summaryFiles:
    input: 
        rules.split_pipe_comb.output,
    params:
        summary_script = ""
    output: 
        "finalreport/metric_summary.xlsx",
        expand("finalreport/summaries/{sample}_report.html", sample=samples)
    shell: "python {params.summary_script}"
