import re

def get_tcr_parent_dir():
    import os
    import csv

    library_csv = getattr(config, "library", "libraries.csv")
    if not os.path.exists(library_csv):
        raise FileNotFoundError(f"Library CSV not found: {library_csv}")

    mapping = {}
    with open(library_csv, newline='') as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        # normalize header names for case-insensitive check
        normalized = [fn.lower().strip() for fn in fieldnames]
        if "tcr" not in normalized or "parent" not in normalized:
            raise ValueError("library CSV must contain at least two columns: 'tcr' and 'parent'")

        for row in reader:
            # case-insensitive access to columns
            row_ci = {k.lower().strip(): (v or "").strip() for k, v in row.items()}
            tcr_key = row_ci.get("tcr", "")
            parent_val = row_ci.get("parent", "")
            if not tcr_key:
                # skip empty tcr entries
                continue
            if tcr_key in mapping:
                raise ValueError(f"Duplicate tcr entry in {library_csv}: '{tcr_key}'")
            mapping[tcr_key] = parent_val

    return mapping

def get_tcr_parent_dir_from_dict(tem_dict, wildcards):
    tcr_key = wildcards if isinstance(wildcards, str) else wildcards.sample
    if tcr_key not in tem_dict:
        raise ValueError(f"TCR key '{tcr_key}' not found in the provided dictionary.")
    return tem_dict[tcr_key]

current_cellranger = program.parsebio_split_pipe
if "split-pipe" in current_cellranger or 'spipe' in current_cellranger:
    split_version = re.search(r"spipe_v(\d+\.\d+\.\d+)", current_cellranger)
    if split_version:
        version_info = split_version.group(1)
        print(f"Split-pipe version: {version_info}")
    else:
        version_info = "unknown"
        print(f"Could not extract version from: {current_cellranger}")
    current_cellranger = f"split-pipe:{version_info}"

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

# Check if TCR analysis is enabled
tcr_enabled = getattr(config, "tcr", False)

immune_genome = "none"
# If the reference indicates hg38 or hg39, force human immune genome
ref_name = getattr(config, "ref", "").lower()
print(ref_name)
if "hg38" in ref_name or "hg39" in ref_name:
    immune_genome = "human"
elif "mm10" in ref_name or "mm39" in ref_name:
    immune_genome = "mouse"

# Validate supported immune genomes when TCR analysis is enabled
if tcr_enabled:
    if immune_genome not in ("human", "mouse"):
        raise ValueError(f"Immune genome for TCR analysis must be 'human' or 'mouse'. Detected: '{immune_genome}'. Check config.parent_dir ('{ref_name}') to determine genome.")

    dict_tcr_parent = get_tcr_parent_dir()

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
            immune_genome = immune_genome,
            parent_dir = lambda wildcards: get_tcr_parent_dir_from_dict(dict_tcr_parent, wildcards),
        shell:
            """
rm -r {params.prefix}; 
{program.parsebio_split_pipe} \
    --mode all \
    --chemistry {params.chemistry} \
    --tcr_analysis \
    --immune_genome {params.immune_genome} \
    --parent_dir {params.parent_dir} \
    --output_dir {params.prefix} \
    --fq1 {input.r1} \
    --fq2 {input.r2} \
    2>{log.err} 1>{log.log}
"""
    rule split_pipe_comb:
        input: 
            expand(rules.split_pipe_all.output, sample=samples)
        params:
            outdir = "split_pipe_comb",
            immune_genome = immune_genome,
            parent_dir = get_tcr_parent_dir_from_dict(dict_tcr_parent, "split_pipe_comb"),
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
        --immune_genome {params.immune_genome} \
        --parent_dir {params.parent_dir} \
        --genome_dir {reference.transcriptome_parsebio} \
        --sublibraries {params.sublibs} \
        --output_dir {params.outdir} \
        2>{log.err} 1>{log.log}
    """

else:
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
        2>{log.err} 1>{log.log}
    """

rule summaryFiles:
    input: 
        rules.split_pipe_comb.output,
    params:
        summary_script = "workflow/scripts/parsebio/generateSummaryFiles.py"
    output: 
        touch("finalreport/metric_summary.xlsx"),
    shell: "python {params.summary_script}"
