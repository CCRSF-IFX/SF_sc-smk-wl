import os 
## RNA
fastqpath4rna=os.path.abspath(config["rna"]["fastqpath"])
projectname4rna = config["rna"]["projectname"]
ref = config["rna"]["ref"]
genome = config["rna"]["genome"]

rule test_sc_rna_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
    params:
        dir4test = os.path.join(outdir_abspath,config["rna"]["projectname"]+"_rna_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_scRNA_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f "{params.dir4test}")
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | python {input.run_snakemake4sc} rna -f {fastqpath4rna} -r {ref} -g {genome} > {output.log} 2>&1
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} rerun 
"""

