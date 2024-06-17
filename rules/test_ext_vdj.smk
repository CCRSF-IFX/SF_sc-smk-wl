import os 
## vdj
fastqpath4vdj= os.path.abspath(config["vdj"]["fastqpath"])
ref4vdj = config["vdj"]["ref"]
genome4vdj = config["vdj"]["genome"]
chain_type = config["vdj"]["chain"]

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
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} vdj --fastq {fastqpath4vdj} --reference {ref4vdj} --genome {genome4vdj} --chain {chain_type} > {output.log} 2>&1
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
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} {fastqpath4vdj} vdj {ref4vdj} --chain {params.chain}  > {output.log} 2>&1
"""
