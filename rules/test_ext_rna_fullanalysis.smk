import os 
## RNA
fastqpath4rna_full = os.path.abspath(config["rna_fullanalysis"]["fastqpath"])
projectname4rna_full = config["rna_fullanalysis"]["projectname"]
ref4rna_full = config["rna_fullanalysis"]["ref"]
genome4rna_full = config["rna_fullanalysis"]["genome"]

rule test_sc_rna_fullanalysis:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
    params:
        dir4test = os.path.join(outdir_abspath,config["rna_fullanalysis"]["projectname"]+"_rna_fullanalysis/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_scRNA_fullanalysis.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f "{params.dir4test}")
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | python {input.run_snakemake4sc} rna -f {fastqpath4rna_full} -r {ref4rna_full} -g {genome4rna_full} --fullanalysis > {output.log} 2>&1
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} rerun
"""
