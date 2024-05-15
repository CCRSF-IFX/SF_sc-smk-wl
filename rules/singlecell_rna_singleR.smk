import os
def load_gene_list(config_ref):
    if config_ref == "mm10":
        return os.path.join(analysis, "workflow/data/marker_genelist/mouse_gene_list.csv")
    elif config_ref == "hg38":
        return os.path.join(analysis, "workflow/data/marker_genelist/human_gene_list.csv")
    else:
        print("Error: Invalid config.ref value. Allowed values are 'mm10' or 'hg38'.")
        exit(1)

gene_list_file = load_gene_list(config.ref)

rule singleR_annotation:
    """
    SingleR annotation
    """
    input:
        seur = rules.seurat_proc.output.seur 
    params:
        outdir = os.path.join(analysis, "{sample}/singleR/"),
    log:
        os.path.join(analysis, "{sample}/singleR/singleR.log")
    output:
        seur = os.path.join(analysis, "{sample}/singleR/seur_10x_cluster_singler.rds")
    container: program.Renv
    shell:
        """
Rscript {analysis}/workflow/scripts/rna/sc_singleR.prod.R --genome={config.ref} --markerList={gene_list_file} --rds={input.seur} --outdir={params.outdir} > {log} 2>&1
"""
