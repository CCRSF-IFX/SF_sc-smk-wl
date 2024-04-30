rule seurat_proc:
    """
    """
    input:
        h5 = rules.count.output 
    params:
        fil_mtx = os.path.join(analysis, "{sample}/outs/filtered_feature_bc_matrix/"),
        outdir = os.path.join(analysis, "{sample}/seurat/"),
    log:
        os.path.join(analysis, "{sample}/seurat/seurat.log")
    output:
        seur = os.path.join(analysis, "{sample}/seurat/seur_10x_cluster_object.rds")
    container: program.Rseurat
    shell:
        """
Rscript {analysis}/workflow/scripts/rna/sc_seurat.prod.R --genome={config.ref} --data.dir={params.fil_mtx}  --outdir={params.outdir} > {log} 2>&1
"""
