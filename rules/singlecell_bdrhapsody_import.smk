import os
import shlex

current_cellranger = os.path.basename(program.bdrhapsody)

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "prep_fastq_folder_ln.smk"
include: "fastqc4QC.smk"
include: "multiqc.smk"
include: "prep_fastq_folder_concat.smk"


def get_sample_config_value(attr, sample, default=None):
    if not hasattr(config, attr):
        return default
    value = getattr(config, attr)
    if isinstance(value, dict):
        return value.get(sample, default)
    return value


def get_sample_path_list(attr, sample):
    value = get_sample_config_value(attr, sample)
    if value in (None, "", False):
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple)):
        return list(value)
    raise ValueError(f"Unsupported value for config.{attr}: {type(value)}")


def get_required_reference_archive(wildcards):
    archive = get_sample_config_value(
        "reference_archive",
        wildcards.sample,
        getattr(reference, "reference_archive_bdrhapsody", ""),
    )
    if archive in (None, ""):
        raise ValueError(
            "reference_archive is required for the BD Rhapsody pipeline. "
            "Set config.reference_archive or define reference.reference_archive_bdrhapsody."
        )
    return archive


def shell_join_args(arg_pairs):
    args = []
    for flag, value in arg_pairs:
        if value in (None, "", False):
            continue
        args.append(f"{flag} {shlex.quote(str(value))}")
    return " ".join(args)


def shell_join_repeated_flag(flag, values):
    return " ".join(f"{flag} {shlex.quote(str(value))}" for value in values)


def build_write_inputs_args(wildcards):
    sample = wildcards.sample
    scalar_pairs = [
        ("--file-field", f"Reference_Archive={get_required_reference_archive(wildcards)}"),
        ("--scalar", f"Run_Name={get_sample_config_value('run_name', sample, sample)}"),
        ("--scalar", f"Sample_Tags_Version={get_sample_config_value('sample_tags_version', sample)}"),
        ("--scalar", f"VDJ_Version={get_sample_config_value('vdj_version', sample)}"),
        ("--scalar", f"Cell_Calling_Data={get_sample_config_value('cell_calling_data', sample)}"),
        (
            "--scalar",
            "Cell_Calling_Bioproduct_Algorithm="
            f"{get_sample_config_value('cell_calling_bioproduct_algorithm', sample)}",
        ),
        (
            "--scalar",
            "Cell_Calling_ATAC_Algorithm="
            f"{get_sample_config_value('cell_calling_atac_algorithm', sample)}",
        ),
        ("--scalar", f"Exact_Cell_Count={get_sample_config_value('exact_cell_count', sample)}"),
        ("--scalar", f"Expected_Cell_Count={get_sample_config_value('expected_cell_count', sample)}"),
        ("--scalar", f"Maximum_Threads={get_sample_config_value('maximum_threads', sample)}"),
        ("--scalar", f"Custom_STAR_Params={get_sample_config_value('custom_star_params', sample)}"),
        ("--scalar", f"Custom_bwa_mem2_Params={get_sample_config_value('custom_bwa_mem2_params', sample)}"),
    ]
    args = shell_join_args(scalar_pairs)

    for attr, field_name in [
        ("targeted_reference", "Targeted_Reference"),
        ("abseq_reference", "AbSeq_Reference"),
        ("supplemental_reference", "Supplemental_Reference"),
        ("reads_atac", "Reads_ATAC"),
    ]:
        values = get_sample_path_list(attr, sample)
        if values:
            rendered = shell_join_repeated_flag(
                "--file-list-field",
                [f"{field_name}={value}" for value in values],
            )
            args = f"{args} {rendered}".strip()

    tag_names = get_sample_path_list("tag_names", sample)
    if tag_names:
        rendered = shell_join_repeated_flag(
            "--scalar-list",
            [f"Tag_Names={value}" for value in tag_names],
        )
        args = f"{args} {rendered}".strip()

    predefined_atac_peaks = get_sample_config_value("predefined_atac_peaks", sample)
    if predefined_atac_peaks not in (None, "", False):
        args = f"{args} --file-field {shlex.quote(f'Predefined_ATAC_Peaks={predefined_atac_peaks}')}".strip()

    bool_map = {
        "exclude_intronic_reads": "Exclude_Intronic_Reads",
        "generate_bam": "Generate_Bam",
        "long_reads": "Long_Reads",
    }
    for attr, field_name in bool_map.items():
        value = get_sample_config_value(attr, sample)
        if value in (True, False):
            args = f"{args} --scalar {shlex.quote(f'{field_name}={str(value).lower()}')}".strip()

    return args


rule write_bdrhapsody_inputs:
    input:
        r1=rules.make_fastq_concat.output.r1,
        r2=rules.make_fastq_concat.output.r2,
    output:
        "{sample}/pipeline_inputs.yml"
    params:
        extra_args=build_write_inputs_args,
        script=os.path.join(analysis, "workflow/scripts/bdrhapsody/write_inputs.py"),
    shell:
        """
python {params.script}     --output {output}     --read {input.r1}     --read {input.r2}     {params.extra_args}
"""


rule count:
    input:
        yaml=rules.write_bdrhapsody_inputs.output,
        r1=rules.make_fastq_concat.output.r1,
        r2=rules.make_fastq_concat.output.r2,
    output:
        web_summary="{sample}/outs/web_summary.html",
        metrics="{sample}/outs/metrics_summary.csv",
        matrix="{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        features="{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        barcodes="{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    log:
        err="run_{sample}_bdrhapsody.err",
        log="run_{sample}_bdrhapsody.log",
    params:
        rawdir=lambda wildcards: os.path.join(analysis, wildcards.sample, "bdrhapsody_raw"),
        script=os.path.join(analysis, "workflow/scripts/bdrhapsody/collect_outputs.py"),
        vendor=program.bdrhapsody,
    shell:
        r"""
set -euo pipefail
rm -rf {params.rawdir} {wildcards.sample}/outs
mkdir -p {params.rawdir}

{params.vendor} pipeline --no-parallel --outdir {params.rawdir} {input.yaml} > {log.log} 2> {log.err} || true

python {params.script}     --sample {wildcards.sample}     --rawdir {params.rawdir}     --outs {wildcards.sample}/outs     --vendor-log {log.log}     --vendor-err {log.err}
"""


rule summaryFiles:
    input:
        expand("{sample}/outs/metrics_summary.csv", sample=samples),
        expand("{sample}/outs/web_summary.html", sample=samples),
    output:
        "finalreport/metric_summary.xlsx",
        expand("finalreport/summaries/{sample}_web_summary.html", sample=samples),
    params:
        summary_script=os.path.join(analysis, "workflow/scripts/bdrhapsody/generateSummaryFiles.py"),
    shell:
        "python {params.summary_script}"
