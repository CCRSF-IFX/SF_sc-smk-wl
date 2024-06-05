import os 

run_info4rna = config["rna"]["archive_run"]
fastqpath4rna=os.path.abspath(config["rna"]["fastqpath"])
projectname4rna = config["rna"]["projectname"]
ref = config["rna"]["ref"]

## no pipe use 'rna' data as an example
rule test_sc_nopipe_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        metadata = os.path.abspath(config["rna"]["metadata"])
    params:
        dir4test = os.path.join(outdir_abspath,config["rna"]["projectname"]+"_nopipe_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_nopipe_default.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4rna} nopipe {ref}  -p {projectname4rna} -a {run_info4rna} {param_test_email} > {output.log} 2>&1
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""
