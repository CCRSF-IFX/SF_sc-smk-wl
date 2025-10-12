
rule prep_fastq_folder_ln:
    """
    Rule to prepare the fastq folder for each sample.
    It creates a symlink to the fastq files.
    The symbolic links will be used by fastQC rule, 
    some of the rules that needs to merge the fastq 
    files from different runs
    """
    output:
        touch("fastq/{sample}/.prepared")
    run:
        sample = wildcards.sample
        prep_fastq_folder_ln(sample)