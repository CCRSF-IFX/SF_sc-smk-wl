import os 

## atac
run_info4atac = config["atac"]["archive_run"]
fastqpath4atac = os.path.abspath(config["atac"]["fastqpath"])
projectname4atac = config["atac"]["projectname"]
ref4atac = config["atac"]["ref"]
metadata = os.path.abspath(config["atac"]["metadata"])

rule test_sc_atac_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        metadata = metadata
    params:
        dir4test = os.path.join(outdir_abspath, config["atac"]["projectname"] + "_atac_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_atac_default.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4atac} atac {ref4atac}  -p {projectname4atac} -a {run_info4atac} {param_test_email} > {output.log} 2>&1
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun
"""

