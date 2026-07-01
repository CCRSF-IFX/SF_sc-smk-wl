import csv
import os
import re
from pathlib import Path
from shutil import copyfile


def read_pixelator_samplesheet_rows(samplesheet=None):
    samplesheet = samplesheet or "samplesheet.pixiome.csv"
    if not samplesheet:
        sys.stderr.write("\npixelator_samplesheet is required for pixiome.\n\n")
        sys.exit(1)
    with open(samplesheet, newline="") as handle:
        reader = csv.DictReader(handle)
        required = ["sample", "sample_alias", "condition", "design", "fastq_1", "fastq_2"]
        fieldnames = reader.fieldnames or []
        missing = [column for column in required if column not in fieldnames]
        if "panel" not in fieldnames and "panel_file" not in fieldnames:
            missing.append("panel or panel_file")
        if missing:
            sys.stderr.write("\nMissing Pixelator samplesheet columns: " + ", ".join(missing) + "\n\n")
            sys.exit(1)
        return list(reader)


def read_pixelator_library_samples():
    libraries = getattr(config, "libraries", "libraries.csv")
    if not libraries:
        return []
    with open(libraries, newline="") as handle:
        reader = csv.DictReader(handle)
        if "sample" in (reader.fieldnames or []):
            return sorted({row["sample"].strip() for row in reader if row.get("sample", "").strip()})
        if "Name" in (reader.fieldnames or []):
            return sorted({row["Name"].strip() for row in reader if row.get("Name", "").strip()})
    return []


# matched_fastq_root and write_fastq_manifest are handled by singlecell_import.smk
# def matched_fastq_root(source_path):
#     source_path = os.path.abspath(source_path)
#     for index, fq_path in enumerate(fastqpath):
#         root = os.path.abspath(fq_path)
#         try:
#             if os.path.commonpath([source_path, root]) == root:
#                 return index, root
#         except ValueError:
#             continue
#     return 0, os.path.dirname(source_path)


# def write_fastq_manifest(manifest_path, samplesheet):
#     columns = [
#         "sample",
#         "read",
#         "lane",
#         "sample_number",
#         "run_name",
#         "run_label",
#         "fastq_root_index",
#         "fastq_root",
#         "source_layout",
#         "source_sample_folder",
#         "source_path",
#         "staged_path",
#     ]
#     rows = []
#     seen = set()
#     for row in read_pixelator_samplesheet_rows(samplesheet):
#         sample = row["sample"].strip()
#         for read, column in [("R1", "fastq_1"), ("R2", "fastq_2")]:
#             source_path = os.path.abspath(row[column].strip())
#             if source_path in seen:
#                 continue
#             seen.add(source_path)
#             index, root = matched_fastq_root(source_path)
#             run_name = run_name_for_fastq_path(index, root)
#             rows.append({
#                 "sample": sample,
#                 "read": read,
#                 "lane": "",
#                 "sample_number": "",
#                 "run_name": run_name,
#                 "run_label": run_name,
#                 "fastq_root_index": index + 1,
#                 "fastq_root": root,
#                 "source_layout": "pixelator_samplesheet",
#                 "source_sample_folder": "",
#                 "source_path": source_path,
#                 "staged_path": source_path,
#             })
#     os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
#     with open(manifest_path, "w") as manifest:
#         manifest.write("\t".join(columns) + "\n")
#         for row in sorted(rows, key=lambda item: (item["sample"], item["source_path"])):
#             manifest.write("\t".join(str(row[column]) for column in columns) + "\n")


if not hasattr(config, "samples"):
    if getattr(config, "libraries", ""):
        samples = read_pixelator_library_samples()
    else:
        samples = sorted({row["sample"].strip() for row in read_pixelator_samplesheet_rows()})


# rule fastq_manifest:
#     input:
#         "samplesheet.pixiome.csv"
#     output:
#         FASTQ_MANIFEST_PATH
#     run:
#         write_fastq_manifest(output[0], input[0])


def get_pixelator_pipeline_dir(wildcards=None):
    pipeline_dir = getattr(config, "nf_pixelator_path", getattr(program, "nf_pixelator_path", None))
    if pipeline_dir:
        return os.path.abspath(str(pipeline_dir))
    return None

current_cellranger = get_pixelator_pipeline_dir()


def get_pixelator_pipeline_version(wildcards=None):
    pipeline_dir = get_pixelator_pipeline_dir(wildcards)
    if pipeline_dir:
        match = re.search(r"v\d+\.\d+\.\d+", Path(pipeline_dir).name)
        if match:
            return match.group(0)
    return getattr(program, "nf_pixelator_version", None)


def get_pixelator_nextflow(wildcards=None):
    return getattr(config, "pixelator_nextflow", getattr(program, "nextflow", "nextflow"))


rule pixelator_generated_samplesheet:
    input:
        lambda wildcards: config.libraries
    output:
        "pixelator_samplesheet.csv"
    params:
        fastqs=",".join(os.path.abspath(path) for path in fastqpath)
    shell:
        "python workflow/scripts/pixiome/generate_samplesheet.py --libraries {input} --fastq-paths {params.fastqs} --output {output}"


rule pixelator_samplesheet:
    input:
        rules.pixelator_generated_samplesheet.output
    output:
        "samplesheet.pixiome.csv"
    shell:
        "python workflow/scripts/pixiome/prepare_samplesheet.py --input {input} --output {output}"


rule pixelator_params:
    input:
        rules.pixelator_samplesheet.output
    output:
        "params.pixiome.yaml"
    run:
        container = getattr(reference, "pixelator_container", None)
        with open(output[0], "w") as handle:
            handle.write(f'input: "{os.path.abspath(input[0])}"\n')
            handle.write(f'outdir: "{analysis}"\n')
            handle.write(f'technology: "{getattr(config, "pixelator_technology", "proxiome-v1")}"\n')
            if container:
                handle.write(f'pixelator_container: "{container}"\n')


rule pixelator_nextflow_config:
    output:
        "nextflow.pixiome.config"
    run:
        source_config = getattr(config, "nf_pixelator_config", "workflow/config/nf_pixelator.config")
        copyfile(source_config, output[0])


rule count:
    input:
        samplesheet=rules.pixelator_samplesheet.output,
        params_yaml=rules.pixelator_params.output,
        nf_config=rules.pixelator_nextflow_config.output,
    output:
        summary="pixelator/experiment-summary.html",
        layout=expand("pixelator/{sample}.layout.pxl", sample=samples),
        reports=expand("pixelator/analysis/{sample}.report.json", sample=samples),
    log:
        err="run_pixelator.err",
        log="run_pixelator.log",
    params:
        pipeline_dir=get_pixelator_pipeline_dir,
        nextflow=get_pixelator_nextflow,
        technology=lambda wildcards: getattr(config, "pixelator_technology", "proxiome-v1"),
        workdir=os.path.join(analysis, "work", "pixiome"),
    shell:
        r"""
set -euo pipefail
if [ -d "{params.pipeline_dir}" ] && [ ! -f "{params.pipeline_dir}/main.nf" ]; then
    echo "ERROR: pixelator pipeline main.nf not found in {params.pipeline_dir}" >&2
    exit 1
fi
if type module >/dev/null 2>&1; then
    module load nextflow || true
    module load singularity || true
fi
mkdir -p "{params.workdir}"
{params.nextflow} -c {input.nf_config} run "{params.pipeline_dir}" \
    -profile singularity \
    --technology "{params.technology}" \
    -params-file {input.params_yaml} \
    -work-dir "{params.workdir}" \
    -resume \
    > {log.log} 2> {log.err}
"""


rule summaryFiles:
    input:
        summary=rules.count.output.summary,
        reports=rules.count.output.reports,
        samplesheet=rules.pixelator_samplesheet.output,
    output:
        "finalreport/metric_summary.xlsx",
        "finalreport/summaries/experiment-summary.html",
    shell:
        "python workflow/scripts/pixiome/generateSummaryFiles.py"


# if external == False:
#     rule archive:
#         input:
#             metadata=report_result,
#             summaryFiles="finalreport/metric_summary.xlsx",
#             fastq_manifest=FASTQ_MANIFEST_PATH,
#             pixelator_summary="pixelator/experiment-summary.html",
#         output:
#             touch("archive_setup.complete")
#         params:
#             fastqs=",".join(os.path.abspath(path) for path in fastqpath),
#             runs=",".join(run_names),
#             meta2json=os.path.join(analysis, "scripts/SF_scDMEarchive/cli/meta2json_single_cell.py"),
#         log:
#             "archive.log"
#         shell:
#             "cd {one_up}; python {params.meta2json} --pipeline {config.pipeline} -m {input.metadata} -r {params.runs} -f {params.fastqs} --fastq-manifest {config.analysis}/{input.fastq_manifest} -c {config.analysis} > {log}"
