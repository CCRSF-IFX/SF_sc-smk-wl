import os

run_info4fixedrna_singleplex = config["fixedrna_singleplex"]["archive_run"]
fastqpath4fixedrna_singleplex = ",".join([os.path.abspath(tempath) for tempath in config["fixedrna_singleplex"]["fastqpath"].split(",")])
projectname4fixedrna_singleplex = config["fixedrna_singleplex"]["projectname"]
ref4fixedrna_singleplex = config["fixedrna_singleplex"]["ref"]
probe_set4singlepelx = config["fixedrna_singleplex"]["probe_set"]

rule test_sc_fixedrna_singleplex_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = os.path.abspath(config["fixedrna_singleplex"]["lib_csv"]),
        metadata = os.path.abspath(config["fixedrna_singleplex"]["metadata"])
    params:
        dir4test = os.path.join(outdir_abspath, config["fixedrna_singleplex"]["projectname"] + "_fixedrna_singleplex_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_fixedrna_singleplex_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f {params.dir4test}/../)
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4fixedrna_singleplex} fixedrna {ref4fixedrna_singleplex} --probe_set {probe_set4singlepelx} --singleplex -p {projectname4fixedrna_singleplex} -a {run_info4fixedrna_singleplex} {param_test_email} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""
