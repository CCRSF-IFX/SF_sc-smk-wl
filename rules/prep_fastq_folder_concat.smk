
# Put this once in your Snakefile if you want all rules to use bash:
# shell.executable("/bin/bash")

rule make_fastq_concat:
    """
    Concatenate FASTQ files for each sample.
    """
    input: 
        rules.prep_fastq_folder_ln.output,
    params:
        prefix = "{sample}",
        prefix2 = lambda wildcards: record_fastqpath[wildcards.sample],
    output:
        r1 = "fastq_concat/{sample}_R1.fastq.gz",
        r2 = "fastq_concat/{sample}_R2.fastq.gz",
    shell: r"""
        set -euo pipefail
        shopt -s nullglob
        export LC_ALL=C
        mkdir -p fastq_concat

        concat_or_link () {{
          local pattern="$1"
          local dest="$2"
          local -a files
          files=($pattern)

          case "${{#files[@]}}" in
            0)
              echo "ERROR: No files match: $pattern" >&2
              exit 1
              ;;
            1)
              echo "Linking single file:"
              printf "  %s -> %s\n" "${{files[0]}}" "$dest"
              ln -sf "${{files[0]}}" "$dest"
              ;;
            *)
              echo "Concatenating ${{#files[@]}} files in this order:"
              printf "  %s\n" "${{files[@]}}"
              cat "${{files[@]}}" > "$dest"
              ;;
          esac
        }}

        # R1
        concat_or_link "{params.prefix2}*_R1_001.fastq.gz" "fastq_concat/{params.prefix}_R1.fastq.gz"

        # R2
        concat_or_link "{params.prefix2}*_R2_001.fastq.gz" "fastq_concat/{params.prefix}_R2.fastq.gz"
    """
