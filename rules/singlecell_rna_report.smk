import os

rule report_rmd:
    input:
        os.path.join(analysis, "{sample}/singleR/seur_10x_cluster_singler.rds")
    params: 
        dir = os.path.join(analysis, '{sample}/'), 
        out = analysis, 
        smp = "{sample}"
    output: 
        os.path.join(analysis, "{sample}_sc_report.html")
    container: program.Renv4rmd
    shell: 
        """
Rscript -e 'rmarkdown::render("{analysis}/workflow/scripts/rna/single_cell_report.Rmd", params = list(project = "{params.smp}", dir = "{params.dir}", sample = "{params.smp}"), output_file = "{output}", output_dir = "{params.out}")'
"""
