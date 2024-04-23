import os
import pandas as pd
img_path = os.path.dirname(os.path.abspath(config.images))
metadata_df = pd.read_csv(config.images)
def check_spatial_files_exist(df, img_path):
    slide_files = set(df['Slide'].tolist())
    missing_slide_files = [file for file in slide_files if not os.path.exists(os.path.join(img_path, file))]
    if missing_slide_files:
        sflog.error("Missing Slide files:" + ",".join(missing_slide_files))
        sys.exit("Pipeline stopped: missing Slides files")
    else:
        sflog.info("All Slide files deteceted.")
    img_files = set(df['Image'].tolist())
    missing_img_files = [file for file in img_files if not os.path.exists(os.path.join(img_path, file))]
    if missing_img_files:
        sflog.error("Missing Image files:" + ",".join(missing_img_files))
        sys.exit("Pipeline stopped: missing Image files")
    else:
        sflog.info("All Image files deteceted.")

sflog.info("Check image file existense...")
check_spatial_files_exist(metadata_df, img_path)

def get_image_data(wildcards):
    print(wildcards)
    sample_info = metadata_df[metadata_df['Sample'] == wildcards[0]]
    print(sample_info)
    if len(sample_info) > 0:
        slide_name = os.path.splitext(sample_info['Slide'].iloc[0])[0]
        slide_file = os.path.join(img_path, sample_info['Slide'].iloc[0])
        area = sample_info['Area'].iloc[0]
        image = os.path.join(img_path, sample_info['Image'].iloc[0])
        #print(" ".join([slide_name, slide_file, area, image]))
        return [slide_name, slide_file, area, image]
    else:
        sflog.error(f"{wildcards} not complete in {config.images}")

flag_probe_set = ""
probe_set = getattr(reference, config.spatial_method, "") 
if probe_set is not "":
    flag_probe_set = f"--probe-set={probe_set}"

current_cellranger = program.spaceranger
 
rule count:
    output: "{sample}/outs/web_summary.html"
    log: err = "run_{sample}_10x_spaceranger_count.err", log ="run_{sample}_10x_spaceranger_count.log"
    params: batch = "-l nodes=1:ppn=16,mem=96gb", prefix = "{sample}", prefix2 = filterFastq, image_info  = get_image_data
    container: program.spaceranger
    shell: "rm -r {params.prefix}; spaceranger count {flag4spaceranger_create_bam}  --id={params.prefix} --sample={params.prefix} --fastqs={params.prefix2} --image={params.image_info[3]} --slide={params.image_info[0]}  --slidefile={params.image_info[1]} --area={params.image_info[2]}  {flag_probe_set}  --transcriptome={reference.transcriptome} --reorient-images=true  2>{log.err} 1>{log.log}"

rule aggregateCSV:
    input: expand("{sample}/outs/web_summary.html", sample=samples)
    output: "AggregatedDatasets.csv"
    params: batch = "-l nodes=1:ppn=1"
    shell: "python workflow/scripts/spatial/generateAggregateCSV.py {analysis}"

rule aggregate:
    input: csv="AggregatedDatasets.csv"
    output: touch("aggregate.complete")
    log: err="run_10x_aggregate.err", log="run_10x_aggregate.log"
    params: batch = "-l nodes=1:ppn=16,mem=96gb"
    container: program.spaceranger 
    shell: "spaceranger aggr --id=AggregatedDatasets --csv={input.csv} --normalize=mapped 2>{log.err} 1>{log.log}"

include: "prep_fastq.smk"
include: "fastqscreen.smk"
include: "kraken.smk"
include: "multiqc.smk"

#rule copyScripts:
#    output: directory("scripts")
#    params: batch = "-l nodes=1:ppn=1"
#    shell: "mkdir scripts; cp -r {program.rscripts} scripts; cp -r {program.pythonscripts} scripts"

rule summaryFiles:
    input: expand("{sample}/outs/web_summary.html", sample=samples)
    output: "finalreport/metric_summary.xlsx", expand("finalreport/summaries/{sample}_web_summary.html", sample=samples)
    params: batch = "-l nodes=1:ppn=1"
    shell: "python workflow/scripts/spatial/generateSummaryFiles.py"
