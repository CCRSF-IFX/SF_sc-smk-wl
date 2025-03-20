import os
import pandas as pd
img_path = os.path.dirname(os.path.abspath(config.images))
metadata_df = pd.read_csv(config.images)
    
def check_spatial_files_exist(df, img_path):
    """
    Validates the existence of required spatial files (Slide, Image, etc.) in the metadata.

    Args:
        df (pd.DataFrame): Metadata DataFrame containing file paths.
        img_path (str): Base directory for image files.

    Raises:
        SystemExit: If any required file is missing.
    """
    # Define the columns to check and their descriptions
    columns_to_check = [
        ("slide", "Slide files"),
        ("image", "Image files"),
        ("cytaimage", "CytAssist image files"),
        ("darkimage", "Dark background fluorescence image files"),
        ("colorizedimage", "Composite colored fluorescence image files"),
    ]

    # Loop through each column and check for missing files
    for column, description in columns_to_check:
        if column in df.columns:
            files = set(df[column].dropna().tolist())
            missing_files = [file for file in files if not os.path.exists(os.path.join(img_path, file))]
            if missing_files:
                sys.exit(f"Pipeline stopped because of missing {description}: " + ", ".join(missing_files))
        else:
            if column == "slide":  # Slide is mandatory
                sys.exit(f"Pipeline stopped: '{column}' column is missing in the metadata file.")
            else:
                print(f"Warning: '{column}' column is missing in the metadata file. Skipping {description} checks.")

    print("All required spatial files are present.")

check_spatial_files_exist(metadata_df, img_path)

def get_params_from_metatab(wildcards):
    """
    Retrieves parameters for Spaceranger from the metadata for a given sample.

    Args:
        wildcards (dict): Wildcards containing the sample name.

    Returns:
        list: A list of parameters to be passed to Spaceranger.
    """
    sr_params = []
    sample_info = metadata_df[metadata_df['sample'] == wildcards.sample]

    if len(sample_info) > 0:
        # Check if the slide is unknown
        unknown_slide = False
        if 'UnknownSlide' in sample_info.columns:
            unknown_slide = sample_info['UnknownSlide'].iloc[0]
            unknown_slide = str(unknown_slide).lower() == "true"

        if unknown_slide:
            sr_params.append("--unknown-slide")
        else:
            # Get slide name and file
            if 'slide' in sample_info.columns and not pd.isna(sample_info['slide'].iloc[0]):
                slide_name = os.path.splitext(sample_info['slide'].iloc[0])[0]
                slide_file = os.path.join(img_path, sample_info['slide'].iloc[0])
                if not os.path.exists(slide_file):
                    sys.exit(f"Slide file {slide_file} does not exist")
                sr_params.append(f"--slide={slide_name} --slidefile={slide_file}")
            else:
                sys.exit(f"Slide information is missing for sample {wildcards.sample}")

            # Get area
            area = sample_info['area'].iloc[0]
            if pd.isna(area) or area == "None":
                sys.exit(f"Area not defined for sample {wildcards.sample}")
            sr_params.append(f"--area={area}")
        with_loupe_alignment = False
        # Define optional columns and their corresponding flags
        optional_columns = [
            ("cytaimage", "--cytaimage"),
            ("image", "--image"),
            ("darkimage", "--darkimage"),
            ("colorizedimage", "--colorizedimage"),
            ("loupe-alignment", "--loupe-alignment"),
        ]

        # Check if at least one image-related column is provided
        image_columns_provided = any(
            column in sample_info.columns and not pd.isna(sample_info[column].iloc[0])
            for column, _ in optional_columns
        )
        if not image_columns_provided:
            sys.exit(f"No image-related columns (CytaImage, Image, DarkImage, ColorizedImage) are provided for sample {wildcards.sample}")

        # Loop through optional columns and add parameters if files exist
        for column, flag in optional_columns:
            if column in sample_info.columns and not pd.isna(sample_info[column].iloc[0]):
                file_path = os.path.join(img_path, sample_info[column].iloc[0])
                if not os.path.exists(file_path):
                    sys.exit(f"{column} file {file_path} does not exist")
                sr_params.append(f"{flag}={file_path}")
                if column == "loupe-alignment":
                    with_loupe_alignment = True
        if with_loupe_alignment:
            sr_params.append("--reorient-images=false")
        else:
            sr_params.append("--reorient-images=true")
        return " ".join(sr_params)
    else:
        sys.exit(f"Sample {wildcards.sample} not found in metadata file {config.images}")

flag_probe_set = ""
probe_set = getattr(reference, config.spatial_method, "") 
if probe_set != "":
    flag_probe_set = f"--probe-set={probe_set}"

current_cellranger = program.spaceranger
 
rule count:
    output: "{sample}/outs/web_summary.html"
    log: err = "run_{sample}_10x_spaceranger_count.err", log ="run_{sample}_10x_spaceranger_count.log"
    params: batch = "-l nodes=1:ppn=16,mem=96gb", prefix = "{sample}", prefix2 = filterFastq, params_from_metatab  = get_params_from_metatab
    container: program.spaceranger
    shell: "rm -r {params.prefix}; spaceranger count {flag4spaceranger_create_bam}  --id={params.prefix} --sample={params.prefix} --fastqs={params.prefix2}  {flag_probe_set}  --transcriptome={reference.transcriptome} {params.params_from_metatab} 2>{log.err} 1>{log.log}"

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
