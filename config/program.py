copydir = "/mnt/ccrsf-ifx/Report_archive/report_archive_singlecell/"
active_scripts = "/mnt/ccrsf-ifx/Software/scripts/bin/"

cellranger     = "docker://ccrsfifx/cellranger:8.0.0"
cellranger_arc = "docker://ccrsfifx/cellranger-arc2.0.2"
spaceranger    = "docker://ccrsfifx/spaceranger:3.0.0"

#atac_cellranger = "/mnt/ccrsf-ifx/Software/tools/GemCode/cellranger-atac-2.1.0/cellranger-atac"
#atac_cellranger_old = "/mnt/ccrsf-ifx/Software/tools/GemCode/cellranger-atac-1.2.0/cellranger-atac"
#arc_cellranger = "/mnt/ccrsf-ifx/Software/tools/GemCode/cellranger-arc-2.0.0/cellranger-arc"

#general_pythonscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/scripts"
##10x RNA specific
#rscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/rna/sctransform/R_scripts_active"
#pythonscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/rna/sctransform/python_scripts"

## PIPseq specific
#pipseeker = "/mnt/ccrsf-ifx/Software/tools/PIPseeker/pipseeker-v2.1.4-linux/pipseeker"
#star_path = "/mnt/ccrsf-ifx/Software/tools/star/STAR-2.7.10a_alpha_220818/source/STAR"

##10x ATAC specific
#atac_rscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/atac/scripts/R_scripts/"
#atac_pythonscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/atac/scripts/python_scripts/"

##10x ATAC old (v1.2.0) specific
#atac_rscripts_old = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/atac_1.2.0/scripts/R_scripts/"
#atac_pythonscripts_old = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/atac_1.2.0/scripts/python_scripts/"

##10x Feature Barcode specific
#fb_pythonscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/featurebarcode/scripts/python_scripts/"

##10x VDJ specific
#vdj_rscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/vdj/scripts/R_scripts/"
#vdj_pythonscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/vdj/scripts/python_scripts/"

##10x Multiome specific
#arc_pythonscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/multiome/scripts/python_scripts/"

##10x Multi specific
#multi_pythonscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/multi/scripts/python_scripts/"

##10x Multi-Sample RNA specific
#multirna_rscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/multi_rna/scripts/R_scripts"
#multirna_pythonscripts = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/multi_rna/scripts/python_scripts"

##Installed packages
#seurat = urd = "export R_LIBS=/mnt/nasapps/development/R/r_libs/3.5.0;export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib;export PATH=/mnt/nasapps/development/python/miniconda3/4.5.12/bin:$PATH;/mnt/nasapps/development/R/3.5.0/bin/R"
#clusterpro = scran = enrichr = nozzle = "export R_LIBS=/mnt/nasapps/development/R/r_libs/3.5.0;/mnt/nasapps/development/R/3.5.0/bin/R"
#pvelocyto = "module load samtools; export PATH=/mnt/ccrsf-ifx/Software/tools/Anaconda/3.6/install/bin:$PATH; export PATH=/opt/nasapps/development/samtools/1.7/bin/:$PATH; export PYTHONPATH=/mnt/ccrsf-ifx/Software/tools/velocyto/lib/python3.6/site-packages/; /mnt/ccrsf-ifx/Software/tools/velocyto/bin/velocyto"
#svelocyto = "export PATH=/opt/nasapps/development/R/3.5.0/bin:/opt/nasapps/development/perl/5.22.1/bin:/mnt/nasapps/development/python/3.7.1/bin:/mnt/nasapps/modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/mnt/nasapps/production/perl/lib/perl5; export R_LIBS=/mnt/nasapps/development/R/r_libs/3.6.0; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/python/miniconda3/4.5.12/bin:$PATH; export LD_LIBRARY_PATH=/mnt/nasapps/development/jags/4.3.0/lib:$LD_LIBRARY_PATH; export PATH=/mnt/nasapps/development/python/3.7.7/bin:$PATH; unset DISPLAY; /mnt/nasapps/development/python/3.7.7/bin/python"
#singler = "export R_LIBS=/mnt/nasapps/development/R/r_libs/3.6.0; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/python/miniconda3/4.5.12/bin:$PATH; export LD_LIBRARY_PATH=/mnt/nasapps/development/jags/4.3.0/lib:$LD_LIBRARY_PATH; /mnt/nasapps/development/R/3.6.0/bin/R"
#nozzleEdit = "export PATH=/mnt/ccrsf-ifx/Software/tools/Anaconda/2.7/install/bin:$PATH;export PYTHONPATH=/mnt/ccrsf-ifx/Software/tools/Anaconda/2.7/install/lib/python2.7/site-packages;python"
#scRNABatchQC = "export R_LIBS=/mnt/nasapps/development/R/r_libs/4.1.1; export PATH=/mnt/nasapps/development/pandoc/2.8.0.1/bin:$PATH; /mnt/nasapps/development/R/4.1.1/bin/R"
#reactomeGSA = "export R_LIBS=/mnt/nasapps/development/R/r_libs/4.1.1; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/pandoc/2.8.0.1/bin:$PATH; export LD_LIBRARY_PATH=/mnt/nasapps/development/jags/4.3.0/lib:$LD_LIBRARY_PATH; /mnt/nasapps/development/R/4.1.1/bin/R"
#immunarch = "module load gcc; export R_LIBS=/mnt/nasapps/development/R/r_libs/4.0.2; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib:$LD_LIBRARY_PATH; export PATH=/mnt/nasapps/development/pandoc/2.8.0.1/bin:$PATH; export LD_LIBRARY_PATH=/mnt/nasapps/development/jags/4.3.0/lib:$LD_LIBRARY_PATH; /mnt/nasapps/development/R/4.0.2/bin/R"
#signac = "export R_LIBS=/mnt/nasapps/development/R/r_libs/3.6.0; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/python/miniconda3/4.5.12/bin:$PATH; /mnt/nasapps/development/R/3.6.0/bin/R"
#signac = "export R_LIBS=/mnt/nasapps/development/R/r_libs/4.0.2; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/pandoc/2.8.0.1/bin:$PATH; export LD_LIBRARY_PATH=/mnt/nasapps/development/jags/4.3.0/lib:$LD_LIBRARY_PATH; /mnt/nasapps/development/R/4.0.2/bin/R"
#python2_7 = "export PATH=/mnt/ccrsf-ifx/Software/tools/Anaconda/2.7/install/bin:$PATH; export PYTHONPATH=/mnt/ccrsf-ifx/Software/tools/Anaconda/2.7/install/lib/python2.7/site-packages; python"
#python3 = "export PATH=/mnt/ccrsf-ifx/Software/tools/Anaconda/3.7/bin:$PATH; python"
#fastmnn = "export R_LIBS=/mnt/nasapps/development/R/r_libs/3.6.0; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/python/miniconda3/4.5.12/bin:$PATH; /mnt/nasapps/development/R/3.6.0/bin/R"
#scanorama = mast = edger = limmatrend = "export R_LIBS=/mnt/nasapps/development/R/r_libs/3.5.0; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/python/miniconda3/4.5.12/bin:$PATH; /mnt/nasapps/development/R/3.5.0/bin/R"
#infercnv = "export R_LIBS=/mnt/nasapps/development/R/r_libs/3.6.0; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/python/miniconda3/4.5.12/bin:$PATH; export LD_LIBRARY_PATH=/mnt/nasapps/development/jags/4.3.0/lib; /mnt/nasapps/development/R/3.6.0/bin/R"
#paga = "export PATH=/opt/nasapps/development/R/3.5.0/bin:/opt/nasapps/development/perl/5.22.1/bin:/mnt/nasapps/development/python/3.7.1/bin:/mnt/nasapps/modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/mnt/nasapps/production/perl/lib/perl5; export R_LIBS=/mnt/nasapps/development/R/r_libs/3.6.0; export LD_LIBRARY_PATH=/mnt/nasapps/development/hdf5/1.10.5/lib; export PATH=/mnt/nasapps/development/python/miniconda3/4.5.12/bin:$PATH; export LD_LIBRARY_PATH=/mnt/nasapps/development/jags/4.3.0/lib:$LD_LIBRARY_PATH; export PATH=/mnt/nasapps/development/python/3.6.5/bin:$PATH; unset DISPLAY; /mnt/nasapps/development/python/3.6.5/bin/python"
#genrich = "/mnt/ccrsf-ifx/Software/tools/genrich/Genrich-0.6/Genrich"
#multiqc = "/mnt/ccrsf-ifx/Software/tools/Anaconda/3.7/bin/multiqc"
multiqc = "docker://multiqc/multiqc:v1.12"
multiqc_conf = "/mnt/ccrsf-ifx/Software/tools/MultiQC/config_multiQC.yaml"
#fastq_screen = "/mnt/ccrsf-ifx/Software/tools/fastq_screen/0.11.4/bin/fastq_screen"
fastq_screen = "docker://quay.io/biocontainers/fastq-screen:0.15.2--pl5321hdfd78af_0"
conf = "workflow/config/fastq_screen_slurm.conf"
#kraken2 = "/mnt/ccrsf-ifx/Software/tools/kraken2/2.0.7/bin/kraken2"
kraken2 = "docker://staphb/kraken2:2.1.3"
kraken2db = "/mnt/ccrsf-ifx/RefGenomes/Kraken2/Kraken2_Complete_0820"
global_container = "docker://ccrsfifx/sc-smk-wl:py1.0.0"
fastqc = "docker://quay.io/biocontainers/fastqc:0.11.8--2"


###############################################################
#x  workflow_img4pipseq = "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/pipseq/img/SingleCell_RNA_PIPseq.png"
