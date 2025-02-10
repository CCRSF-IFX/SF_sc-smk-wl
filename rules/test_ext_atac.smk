import os 

## atac
fastqpath4atac = os.path.abspath(config["atac"]["fastqpath"])
projectname4atac = config["atac"]["projectname"]
ref4atac = config["atac"]["ref"]
genome_atac = config["atac"]["genome"]

rule test_sc_atac_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
    params:
        dir4test = os.path.join(outdir_abspath, config["atac"]["projectname"] + "_atac_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_atac_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f {params.dir4test}/../)
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc}  atac  -r {ref4atac} -f {fastqpath4atac} -g {genome_atac}  > {output.log} 2>&1
printf "chemistry='ARC-v1'\\n" >> config.py 
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} rerun
"""

