import os 

## multiome
run_info4multiome = config["multiome"]["archive_run"]
fastqpath4multiome = ",".join([os.path.abspath(tempath) for tempath in config["multiome"]["fastqpath"].split(",")])
projectname4multiome = config["multiome"]["projectname"]
ref4multiome = config["multiome"]["ref"]
lib_csv = os.path.abspath(config["multiome"]["lib_csv"]), 
metadata = os.path.abspath(config["atac"]["metadata"])

rule test_sc_multiome_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = lib_csv, 
        metadata = metadata
    params:
        dir4test = os.path.join(outdir_abspath, config["multiome"]["projectname"] + "_multiome_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_multiome_default.log")
    shell:
        """
mkdir -p {params.dir4test}
rm -rf $(readlink -f {params.dir4test}/../)
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4multiome} multiome {ref4multiome}  -p {projectname4multiome} -a {run_info4multiome} {param_test_email} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun
"""

rule test_sc_multiome_exclude_exons:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = config["multiome"]["lib_csv"]
    params:
        dir4test = os.path.join(outdir_abspath, config["multiome"]["projectname"] + "_multiome_exclude_exons/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_multiome_exclude_exons.log")
    shell:
        """
mkdir -p {params.dir4test}
rm -rf $(readlink -f {params.dir4test}/../)
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4multiome} multiome {ref4multiome} --exclude-introns  -p {projectname4multiome} -a {run_info4multiome} {param_test_email} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

