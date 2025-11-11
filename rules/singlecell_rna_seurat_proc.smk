if config.pipeline == "parsebio":
    rule seurat_proc:
        """
        """
        input:
            split_pipe_comb_output = rules.split_pipe_comb.output,
        params:
            fil_mtx = lambda wildcards: parsebio_samples[wildcards.parsebio_sample][0],
            outdir = lambda wildcards: parsebio_samples[wildcards.parsebio_sample][1],
        log:
            os.path.join(analysis, "split_pipe_comb/{parsebio_sample}/seurat/seurat.log")
        output:
            seur = os.path.join(analysis, "split_pipe_comb/{parsebio_sample}/seurat/seur_10x_cluster_object.rds")
        container: program.Renv
        shell:
            """
    Rscript {analysis}/workflow/scripts/rna/sc_seurat.prod.R --genome={config.ref} --data.dir={params.fil_mtx}  --outdir={params.outdir} > {log} 2>&1
    """
else:
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
        container: program.Renv
        shell:
            """
    Rscript {analysis}/workflow/scripts/rna/sc_seurat.prod.R --genome={config.ref} --data.dir={params.fil_mtx}  --outdir={params.outdir} > {log} 2>&1
    """
