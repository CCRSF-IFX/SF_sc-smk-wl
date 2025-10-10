copydir = "/mnt/ccrsf-ifx/Report_archive/report_archive_singlecell/"
active_scripts = "/mnt/ccrsf-ifx/Software/scripts/bin/"

# Last checking: Mar 5, 2024
# Avaliable versions of cellranger:
# 9.0.1;9.0.0;8.0.1;8.0.0;7.2.0;7.1.0; 
# https://hub.docker.com/r/ccrsfifx/cellranger/tags
cellranger      = "docker://ccrsfifx/cellranger:9.0.1"
# Available versions of cellranger-atac:
# 2.1.0
cellranger_atac = "docker://ccrsfifx/cellranger-atac:2.1.0"
# Available versions of cellranger-arc:
# 2.0.2
cellranger_arc  = "docker://ccrsfifx/cellranger-arc:2.0.2"
# Available versions of spaceranger:
# 4.0.1;3.1.3;3.0.0
spaceranger     = "docker://ccrsfifx/spaceranger:4.0.1"
Renv            = "docker://ccrsfifx/sc-smk-wl:r1.0.0"
Renv4rmd        = "docker://ccrsfifx/sc-rmd-report:1.0.0"
pipseeker       = "docker://public.ecr.aws/w3e1n2j6/fluent-pipseeker:3.3.0"

multiqc         = "docker://multiqc/multiqc:v1.12"
multiqc_conf    = "/mnt/ccrsf-ifx/Software/tools/MultiQC/config_multiQC.yaml"
fastq_screen    = "docker://quay.io/biocontainers/fastq-screen:0.15.2--pl5321hdfd78af_0"
conf = "workflow/config/fastq_screen_slurm.conf"
kraken2 = "docker://staphb/kraken2:2.1.3"
kraken2db = "/mnt/ccrsf-ifx/RefGenomes/Kraken2/Kraken2_Complete_0820"
global_container = "docker://ccrsfifx/sc-smk-wl:py1.0.0"
fastqc = "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
picard = "docker://quay.io/biocontainers/picard:3.4.0--hdfd78af_0"

tenx_cloud_token_path = "/mnt/ccrsf-ifx/RefGenomes/SingleCell_REF/10X/CLI/txg"
