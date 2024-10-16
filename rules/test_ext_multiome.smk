import os 

## multiome
fastqpath4multiome = [os.path.abspath(tempath) for tempath in config["multiome"]["fastqpath"].split(",")]
projectname4multiome = config["multiome"]["projectname"]
ref4multiome = config["multiome"]["ref"]
genome4multiome = config["multiome"]["genome"]
lib_csv = os.path.abspath(config["multiome"]["lib_csv"]), 

rule test_sc_multiome_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = lib_csv, 
    params:
        dir4test = os.path.join(outdir_abspath, config["multiome"]["projectname"] + "_multiome_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_multiome_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f "{params.dir4test}")
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} multiome -f {fastqpath4multiome} -r {ref4multiome}  -g {genome4multiome} --library_config {input.lib_csv} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} rerun
"""

