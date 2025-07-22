def picard_rnametrics_output():
    """
    Function to conditionally include output files based on the pipeline.
    This function checks the pipeline type and returns the appropriate output files.
    If the pipeline is "rna" or "pipseq", it returns the rnametrics file.
    """
    if config.pipeline == "rna" or config.pipeline == "pipseq":
        return expand("QC/Sample_{sample}/{sample}_rnametrics.txt", sample=samples)
    else:
        return []

rule multiqc:
    input: 
        expand("QC/Sample_{sample}/{sample}_R2_screen.png", sample=samples), 
        expand("QC/Sample_{sample}/{sample}.kraken.report.txt", sample=samples),
        expand(rules.fastqc4QC.output, sample=samples)
        #picard_rnametrics_output()
    output: 
        "QC/" + project_name + "_multiqc.html"
    container: 
        program.multiqc
    shell: 
        """
multiqc -f -c {program.multiqc_conf} -n {output} ./QC
"""
