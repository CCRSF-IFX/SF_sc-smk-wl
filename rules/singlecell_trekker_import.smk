"""Rules and configuration for Takara TrekkerU_CX primary analysis."""

import csv
import os
import subprocess
import sys
from pathlib import Path

import config
import program
import reference
from trekker_utils import (
    MATRIX_FILES,
    cellranger_fastqs,
    expected_report,
    load_libraries,
    validate_installation,
    write_vendor_samplesheet,
)


analysis = Path(config.analysis).resolve()
libraries_value = getattr(config, "libraries", "")
if not libraries_value:
    raise ValueError(
        "libraries is required for the Trekker pipeline; set it to the path of "
        "the completed libraries.csv in config.py"
    )
libraries = Path(libraries_value).expanduser().resolve()
output_root = Path(
    getattr(config, "trekker_output_dir", analysis / "trekker_out")
).resolve()
installation = Path(program.trekker_installation).resolve()
conda_environment = Path(program.trekker_conda_env).resolve()
cellranger_output_root = (analysis / "cellranger").resolve()
trekker_fastq_root = Path(
    getattr(config, "trekker_fastq_dir", analysis / "trekker_fastqs")
).resolve()
records = load_libraries(
    libraries,
    unaligned=config.unaligned,
    transcriptome=reference.transcriptome,
    cellranger_output_root=cellranger_output_root,
    trekker_fastq_root=trekker_fastq_root,
)
cellranger_records = records

# Takara specifies a host with at least 12 cores even though 8 positioning
# cores is the recommended samplesheet value.
trekker_threads = max(12, max(int(record["cores"]) for record in records.values()))
validate_installation(
    installation,
    {record["profile"] for record in records.values()},
    conda_environment=conda_environment,
)


def get_sample_record(wildcards):
    record = records.get(wildcards.sample)
    if record is None:
        raise ValueError(f"Unknown Trekker sample: {wildcards.sample}")
    return record


def required_inputs(wildcards):
    record = get_sample_record(wildcards)
    return [
        str(libraries),
        record["barcode_file"],
        record["fastq_1"],
        record["fastq_2"],
    ] + [str(Path(record["sc_outdir"]) / name) for name in MATRIX_FILES]


def cellranger_record(wildcards):
    record = cellranger_records.get(wildcards.sample)
    if record is None:
        raise ValueError(f"Sample does not request integrated Cell Ranger: {wildcards.sample}")
    return record


def cellranger_input_fastqs(wildcards):
    return cellranger_fastqs(cellranger_record(wildcards))


def report_path(record):
    return str(expected_report(output_root, record))


reports = [report_path(record) for record in records.values()]
metrics = [path.replace("_Trekker_Report.html", "_summary_metrics.csv") for path in reports]
completion_files = [str(Path(path).parent / ".trekker_complete") for path in reports]
cellranger_web_summaries = [
    str(cellranger_output_root / sample / "outs" / "web_summary.html")
    for sample in cellranger_records
]
cellranger_metrics = [
    str(cellranger_output_root / sample / "outs" / "metrics_summary.csv")
    for sample in cellranger_records
]
final_report = Path("finalreport")
final_metric_summary = str(final_report / "metric_summary.xlsx")
final_cellranger_summaries = [
    str(final_report / "summaries" / f"{sample}_web_summary.html")
    for sample in cellranger_records
]
final_trekker_summaries = [
    str(final_report / "summaries" / f"{sample}_Trekker_Report.html")
    for sample in records
]

# The common reporting/archive rules use ``samples`` for final analysis units,
# while the FASTQ manifest must retain every GEX and Trekker library prefix.
with libraries.open(newline="", encoding="utf-8-sig") as libraries_handle:
    fastq_samples = sorted(
        {row["Sample"].strip() for row in csv.DictReader(libraries_handle)}
    )

# Run the standard sequencing-QC suite for every source library prefix. The
# workflow later restores ``samples`` to the combined Trekker analysis names.
samples = fastq_samples
include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "prep_fastq_folder_ln.smk"
include: "fastqc4QC.smk"
include: "multiqc.smk"

samples = sorted(records)
aggregate = False
current_cellranger = installation.name


localrules: vendor_samplesheet


rule summaryFiles:
    input:
        libraries=str(libraries),
        cellranger_metrics=cellranger_metrics,
        cellranger_reports=cellranger_web_summaries,
        trekker_metrics=metrics,
        trekker_reports=reports,
    output:
        workbook=final_metric_summary,
        cellranger_reports=final_cellranger_summaries,
        trekker_reports=final_trekker_summaries,
    params:
        script="workflow/scripts/trekker/generateSummaryFiles.py",
    threads:
        1
    resources:
        mem_mb=4000,
        runtime_min=60,
    shell:
        """
        python {params.script:q} \
            --libraries {input.libraries:q} \
            --cellranger-root {cellranger_output_root:q} \
            --trekker-root {output_root:q} \
            --output-dir {final_report:q}
        """


rule cellranger_count:
    input:
        libraries=str(libraries),
        fastqs=cellranger_input_fastqs,
    output:
        web=str(cellranger_output_root / "{sample}" / "outs" / "web_summary.html"),
        barcodes=str(cellranger_output_root / "{sample}" / "outs" / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"),
        features=str(cellranger_output_root / "{sample}" / "outs" / "filtered_feature_bc_matrix" / "features.tsv.gz"),
        matrix=str(cellranger_output_root / "{sample}" / "outs" / "filtered_feature_bc_matrix" / "matrix.mtx.gz"),
    log:
        stdout="cellranger_logs/{sample}.log",
        stderr="cellranger_logs/{sample}.err",
    threads:
        32
    resources:
        mem_mb=773497,
        runtime_min=7200,
    params:
        cellranger_id=lambda wildcards: wildcards.sample,
        fastq_dir=lambda wildcards: ",".join(
            cellranger_record(wildcards)["cellranger_fastq_dirs"]
        ),
        fastq_sample=lambda wildcards: ",".join(
            cellranger_record(wildcards)["cellranger_fastq_samples"]
        ),
        transcriptome=lambda wildcards: cellranger_record(wildcards)["transcriptome"],
        output_dir=lambda wildcards: str(cellranger_output_root / wildcards.sample),
    singularity:
        program.cellranger
    shell:
        """
        mkdir -p {cellranger_output_root:q} cellranger_logs
        rm -rf -- {params.output_dir:q}
        cellranger count \
            --id={params.cellranger_id:q} \
            --output-dir={params.output_dir:q} \
            --fastqs={params.fastq_dir:q} \
            --sample={params.fastq_sample:q} \
            --transcriptome={params.transcriptome:q} \
            --chemistry=SC3Pv4 \
            --include-introns=true \
            --create-bam=true \
            --localcores={threads} \
            --localmem=740 \
            >{log.stdout:q} 2>{log.stderr:q}
        """


rule trekker_fastq_pair:
    input:
        libraries=str(libraries),
        r1=lambda wildcards: get_sample_record(wildcards)["trekker_r1_fastqs"],
        r2=lambda wildcards: get_sample_record(wildcards)["trekker_r2_fastqs"],
    output:
        r1=str(trekker_fastq_root / "{sample}_R1.fastq.gz"),
        r2=str(trekker_fastq_root / "{sample}_R2.fastq.gz"),
    threads:
        1
    resources:
        mem_mb=4000,
        runtime_min=240,
    shell:
        r"""
        set -euo pipefail
        mkdir -p {trekker_fastq_root:q}

        concat_or_link() {{
            destination="$1"
            shift
            if [ "$#" -eq 1 ]; then
                ln -sfn "$1" "$destination"
            else
                cat "$@" > "$destination"
            fi
        }}

        concat_or_link {output.r1:q} {input.r1:q}
        concat_or_link {output.r2:q} {input.r2:q}
        """


rule vendor_samplesheet:
    input:
        libraries=str(libraries),
        configuration="config.py",
    output:
        "trekker_samplesheets/{sample}.csv",
    run:
        record = records.get(wildcards.sample)
        if record is None:
            raise ValueError(f"Unknown Trekker sample: {wildcards.sample}")
        write_vendor_samplesheet(record, output[0])


rule trekker_ucx:
    input:
        samplesheet="trekker_samplesheets/{sample}.csv",
        required=required_inputs,
    output:
        report=str(output_root / "_{sample}" / "trekker_{sample}" / "output" / "{sample}_Trekker_Report.html"),
        metrics=str(output_root / "_{sample}" / "trekker_{sample}" / "output" / "{sample}_summary_metrics.csv"),
        complete=str(output_root / "_{sample}" / "trekker_{sample}" / "output" / ".trekker_complete"),
    log:
        "trekker_logs/{sample}.log",
    threads:
        trekker_threads
    resources:
        mem_mb=320000,
        runtime_min=5760,
    params:
        runtime=lambda wildcards: str(analysis / ".trekker_runtime" / wildcards.sample),
    run:
        get_sample_record(wildcards)
        command = [
            sys.executable,
            "workflow/scripts/trekker/trekker_utils.py",
            input.samplesheet,
            "--installation",
            str(installation),
            "--output-root",
            str(output_root),
            "--runtime-directory",
            params.runtime,
            "--conda-environment",
            str(conda_environment),
        ]
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(log[0], "w", encoding="utf-8") as log_handle:
            subprocess.run(
                command,
                check=True,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )
        Path(output.complete).touch()
