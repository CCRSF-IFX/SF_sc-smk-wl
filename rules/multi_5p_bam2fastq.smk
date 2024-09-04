import subprocess
import re
indir = config["indir"]
outdir = config["outdir"]
LIB, SAMPLE = glob_wildcards(os.path.join(indir,"{lib}/outs/per_sample_outs/{sample}/count/sample_alignments.bam"))
LIB = list(set(LIB))
print(LIB)
print(SAMPLE)

# Conversion function in Python
def convert_hto_to_a03(wildcards):
    # Extract the numeric part from the HTO code
    numeric_part = wildcards.sample.split('_')[-1]
    # Convert it to two digits by adding a leading zero if necessary
    two_digit_number = numeric_part.zfill(2)
    # Combine with the prefix "A03"
    result = "A03" + two_digit_number
    return result

def prep_gex_folder(wildcards):
    bam = os.path.join(indir, f"{wildcards.lib}/outs/per_sample_outs/{wildcards.sample}/count/sample_alignments.bam")
    #LIB_ID = subprocess.check_output(f'samtools view -H {bam}' + r'|grep -oP '"library_id":\K\d+(?=,"library_type":"Gene Expression")', shell=True).decode()
    samtools_output = subprocess.check_output(f'samtools view -H {bam}', shell=True).decode()
    #print(samtools_output)
    pattern = r'@CO\s+library_info:\{"library_id":(\d+),"library_type":"Gene Expression".*\}'
    match =re.search(pattern, samtools_output)
    # If a match is found, extract the library_id
    LIB_ID = int(match.group(1))
    fq_path = f'{outdir}/bamtofastq/{wildcards.lib}_{wildcards.sample}/{wildcards.lib}_{LIB_ID}_1_*/*fastq.gz'
    #print(fq_path)
    fq_list = subprocess.check_output(f'readlink -f {fq_path}', shell=True).decode().split("\n") 
    cmd4mkdir = {}
    cmd4ln_s  = []
    a03_id = convert_hto_to_a03(wildcards) 
    print(a03_id)
    # 22NK2FLT3/outs/fastq_path/22NK2FLT3
    for fq in fq_list:
        flowcell_id = os.path.dirname(fq).split("_")[-1]
        if len(flowcell_id) == 9: 
            basename = os.path.basename(fq) 
            basename = re.sub("bamtofastq_S", f'{wildcards.lib}_{a03_id}_S', basename)
            tem_cmd = f"mkdir -p {outdir}/raw_fq/{flowcell_id}/outs/fastq_path/{flowcell_id}/{wildcards.lib}_{a03_id}/"
            cmd4mkdir[tem_cmd] = 1
            tem_cmd_lns = f"ln -s {fq} {outdir}/raw_fq/{flowcell_id}/outs/fastq_path/{flowcell_id}/{wildcards.lib}_{a03_id}/{basename}"
            cmd4ln_s.append(tem_cmd_lns)
    return "\n".join(cmd4mkdir.keys()) + "\n" + "\n".join(cmd4ln_s)

rule all:
    input:
        expand(os.path.join(outdir, "bamtofastq/{lib}-{sample}.log"), lib=LIB,  sample = SAMPLE),
        expand(os.path.join(outdir, "raw_fq/{lib}-{sample}.log"), lib=LIB,  sample = SAMPLE),

rule bamtofastq:
    input:
        bam = os.path.join(indir,"{lib}/outs/per_sample_outs/{sample}/count/sample_alignments.bam") 
    params:
        outdir = os.path.join(outdir, "bamtofastq/{lib}_{sample}/")
    output:
        log = os.path.join(outdir, "bamtofastq/{lib}-{sample}.log")
    shell:
        """
/mnt/ccrsf-ifx/Software/tools/GemCode/cellranger-8.0.1/lib/bin/bamtofastq {input.bam} {params.outdir} > {output.log}
"""

rule prep_raw_fq_folder:
    input:
        bam = os.path.join(indir,"{lib}/outs/per_sample_outs/{sample}/count/sample_alignments.bam"),
        log = os.path.join(outdir, "bamtofastq/{lib}-{sample}.log")
    params:
        outdir = os.path.join(outdir, "bamtofastq/{lib}_{sample}/"),
        #a03_code = convert_hto_to_a03
        cmd = prep_gex_folder
    output:
        log = os.path.join(outdir, "raw_fq/{lib}-{sample}.log")
    shell:
        """
module load samtools
#LIB_ID=$(samtools view -H {input.bam} |grep -oP '"library_id":\\K\\d+(?=,"library_type":"Gene Expression")') 
#readlink -f {outdir}/bamtofastq/{wildcards.lib}_{wildcards.sample}/{wildcards.lib}_${{LIB_ID}}_1_*/*fastq.gz
{params.cmd}
touch {output.log}
"""

#print(wildcards)
