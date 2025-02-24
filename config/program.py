copydir = "/mnt/ccrsf-ifx/Report_archive/report_archive_singlecell/"
active_scripts = "/mnt/ccrsf-ifx/Software/scripts/bin/"

cellranger      = "docker://ccrsfifx/cellranger:8.0.1"
cellranger_atac = "docker://ccrsfifx/cellranger-atac:2.1.0"
cellranger_arc  = "docker://ccrsfifx/cellranger-arc:2.0.2"
spaceranger     = "docker://ccrsfifx/spaceranger:3.0.0"
Renv            = "docker://ccrsfifx/sc-smk-wl:r1.0.0"
Renv4rmd        = "docker://ccrsfifx/sc-rmd-report:1.0.0"

multiqc         = "docker://multiqc/multiqc:v1.12"
multiqc_conf    = "/mnt/ccrsf-ifx/Software/tools/MultiQC/config_multiQC.yaml"
fastq_screen    = "docker://quay.io/biocontainers/fastq-screen:0.15.2--pl5321hdfd78af_0"
pipseeker       = "docker://public.ecr.aws/w3e1n2j6/fluent-pipseeker:3.3.0"
conf = "workflow/config/fastq_screen_slurm.conf"
kraken2 = "docker://staphb/kraken2:2.1.3"
kraken2db = "/mnt/ccrsf-ifx/RefGenomes/Kraken2/Kraken2_Complete_0820"
global_container = "docker://ccrsfifx/sc-smk-wl:py1.0.0"
fastqc = "docker://quay.io/biocontainers/fastqc:0.11.8--2"
