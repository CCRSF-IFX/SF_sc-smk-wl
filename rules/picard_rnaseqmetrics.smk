import os

rule picard_rnametrics:
    input:
        bam = rules.count.output,
        refflat = "{sample}/path/to/refFlat.txt", 
        rRNA_interval = "{sample}/path/to/rRNA_interval_list.txt", 
    params:
        bam = os.path.join(analysis, "{sample}/outs/possorted_genome_bam.bam"), 
        bai = os.path.join(analysis, "{sample}/outs/possorted_genome_bam.bam.bai"),
        bam_link = os.path.join(analysis, "QC/Sample_{sample}/{sample}.bam"),
        bai_link = os.path.join(analysis, "QC/Sample_{sample}/{sample}.bam.bai"),
    output:
        rna_metrics = "QC/Sample_{sample}/{sample}_rnametrics.txt"
    container: 
        program.picard
    log:
        "QC/Sample_{sample}/{sample}_rnametrics.log"
    shell:
        """
        ln -s {input.bam} {params.bam_link}
        ln -s {params.bai} {params.bai_link}
        picard CollectRnaSeqMetrics \
            INPUT={input.bam_link} \
            OUTPUT={output.rna_metrics} \
            STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
            REF_FLAT={input.refflat} \
            RIBOSOMAL_INTERVALS={input.genome} \
            VALIDATION_STRINGENCY=SILENT \
            RIBOSOMAL_INTERVALS={input.rRNA_interval} \
            STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND
            > {log} 2>&1
        rm {params.bam_link}
        rm {params.bai_link}
        """