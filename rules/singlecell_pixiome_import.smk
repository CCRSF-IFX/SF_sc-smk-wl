import os
import re
from pathlib import Path
from shutil import copyfile

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
        container = getattr(config, "pixelator_container", getattr(program, "pixelator_container", None))
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
        source_config = getattr(config, "pixelator_config", "")
        if source_config:
            copyfile(source_config, output[0])
        else:
            Path(output[0]).write_text(
                "process {\n"
                "    executor = 'slurm'\n"
                "    queue = 'largemem'\n"
                "    clusterOptions = '--nodes=1 --exclusive --mem=0'\n"
                "}\n"
                "executor {\n"
                "    queueSize = 50\n"
                "    pollInterval = '30 sec'\n"
                "}\n"
                "report { overwrite = true }\n"
                "timeline { overwrite = true }\n"
                "trace { overwrite = true }\n"
                "dag { overwrite = true }\n"
            )


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
