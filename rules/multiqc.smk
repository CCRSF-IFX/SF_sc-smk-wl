rule multiqc:
    input: expand("QC/Sample_{sample}/{sample}_R2_screen.png", sample=samples), expand("QC/Sample_{sample}/{sample}.kraken.report.txt", sample=samples)
    output: "QC/" + project_name + "_multiqc.html"
    params: batch = "-l nodes=1:ppn=4,mem=16g"
    container: program.multiqc
    shell: "multiqc -f -c {program.multiqc_conf} -n {output} ./QC"

