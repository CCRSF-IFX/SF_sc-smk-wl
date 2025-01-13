import os

fastqpath4fixedrna_singleplex = [os.path.abspath(tempath) for tempath in config["fixedrna_singleplex"]["fastqpath"].split(",")]
projectname4fixedrna_singleplex = config["fixedrna_singleplex"]["projectname"]
genome4fixedrna_singleplex = config["fixedrna_singleplex"]["genome"]
probe_set4singlepelx = config["fixedrna_singleplex"]["probe_set"]
lib_csv4singleplex = config["fixedrna_singleplex"]["lib_csv"]

rule test_sc_fixedrna_singleplex_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = os.path.abspath(config["fixedrna_singleplex"]["lib_csv"]),
    params:
        dir4test = os.path.join(outdir_abspath, config["fixedrna_singleplex"]["projectname"] + "_fixedrna_singleplex_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_fixedrna_singleplex_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f {params.dir4test}/)
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} fixedrna -f {fastqpath4fixedrna_singleplex} -g {genome4fixedrna_singleplex} --probe_set {probe_set4singlepelx} --library_config {lib_csv4singleplex} --singleplex > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} rerun 
"""
