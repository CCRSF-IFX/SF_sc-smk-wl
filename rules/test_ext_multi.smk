import os

run_info4multi = config["multi"]["archive_run"]
fastqpath4multi=config["multi"]["fastqpath"]
genome4multi = config["multi"]["genome"]
ref4multi = config["multi"]["ref"]

rule test_sc_multi_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = os.path.abspath(config["multi"]["lib_csv"]),
    params:
        dir4test = os.path.join(outdir_abspath, config["multi"]["projectname"] + "_multi_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_multi_default.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} multi -f {fastqpath4multi} -r {ref4multi} -g {genome4multi} > {output.log} 2>&1
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

