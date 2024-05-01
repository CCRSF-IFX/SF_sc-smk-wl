aggregate = False

current_cellranger = program.cellranger

# Force the analysis to be carried out for a particular chain type
chain_type = getattr(config, "chain", "auto")
chain_params = f"--chain {chain_type}" if chain_type != "auto" else ""

rule count:
    output: "{sample}/outs/web_summary.html"
    log: err = "run_{sample}_10x_cellranger_vdj.err", log ="run_{sample}_10x_cellranger_vdj.log"
    params: batch = "-l nodes=1:ppn=16,mem=96gb", prefix = "{sample}", prefix2 = filterFastq,
    container: program.cellranger
    shell: "rm -r {params.prefix}; cellranger vdj --id={params.prefix} --sample={params.prefix} {chain_params} --reference={reference.vdj_reference} --fastqs={params.prefix2} 2>{log.err} 1>{log.log}"

rule summaryFiles:
    input: 
        expand("{sample}/outs/web_summary.html", sample=samples)
    output: 
        "finalreport/metric_summary.xlsx", 
        expand("finalreport/summaries/{sample}_web_summary.html", sample=samples)
    shell: 
        """
python {analysis}/workflow/scripts/vdj/generateSummaryFiles.py 
"""
