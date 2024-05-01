import os 
## vdj
run_info4vdj = config["vdj"]["archive_run"]
fastqpath4vdj= os.path.abspath(config["vdj"]["fastqpath"])
projectname4vdj = config["vdj"]["projectname"]
ref4vdj = config["vdj"]["ref"]


rule test_sc_vdj_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
    params:
        dir4test = os.path.join(outdir_abspath, config["vdj"]["projectname"] + "_vdj_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_vdj_default.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} {fastqpath4vdj} vdj {ref4vdj} -p {projectname4vdj} -a {run_info4vdj} {param_test_email} > {output.log} 2>&1
"""

rule test_sc_vdj_chain:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
    params:
        dir4test = os.path.join(outdir_abspath, config["vdj"]["projectname"] + "_vdj_chain/Analysis/"),
        chain = config["vdj"]["chain"]
    output:
        log = os.path.join(outdir_abspath, "test_sc_vdj_chain.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} {fastqpath4vdj} vdj {ref4vdj} --chain {params.chain} -p {projectname4vdj} -a {run_info4vdj} {param_test_email} > {output.log} 2>&1
"""
