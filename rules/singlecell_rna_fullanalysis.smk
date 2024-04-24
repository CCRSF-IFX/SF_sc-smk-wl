rule seurat_proc:
    """
    """
    input:
        h5 = rules.count.output 
    params:
        fil_mtx = os.path.join(analysis, "{sample}/outs/filtered_feature_bc_matrix/"),
        outdir = os.path.join(analysis, "{sample}/outs/seurat/"),
    output:
        log = os.path.join(analysis, "{sample}/outs/seurat/seurat.log")
    container: program.Rseurat 
    shell:
        """
Rscript {analysis}/scripts/rna/sc_seurat.prod.R --genome={config.ref} --data.dir={params.fil_mtx}  --outdir={params.outdir} > {output.log} 2>&1
"""
