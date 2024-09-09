import os

run_info4multi = config["multi"]["archive_run"]
fastqpath4multi = ",".join([os.path.abspath(tempath) for tempath in config["multi"]["fastqpath"].split(",")])
projectname4multi = config["multi"]["projectname"]
ref4multi = config["multi"]["ref"]

rule test_sc_multi_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = os.path.abspath(config["multi"]["lib_csv"]),
        metadata = os.path.abspath(config["multi"]["metadata"])
    params:
        dir4test = os.path.join(outdir_abspath, config["multi"]["projectname"] + "_multi_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_multi_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f {params.dir4test}/../)
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4multi} multi {ref4multi}  -p {projectname4multi} -a {run_info4multi} {param_test_email} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

rule test_sc_multi_exclude_introns:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = os.path.abspath(config["multi"]["lib_csv"])
    params:
        dir4test = os.path.join(outdir_abspath, config["multi"]["projectname"] + "_multi_exclude_introns/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_multi_exclude_introns.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f {params.dir4test}/../)
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4multi} multi {ref4multi} --exclude-introns -p {projectname4multi} -a {run_info4multi} {param_test_email} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

rule test_sc_multi_force_cell:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = os.path.abspath(config["multi"]["lib_csv"])
    params:
        dir4test = os.path.join(outdir_abspath, config["multi"]["projectname"] + "_multi_forcecell/Analysis/"),
        cell_number = config["multi"]["cell_number"]
    output:
        log = os.path.join(outdir_abspath, "test_sc_multi_forcecell.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4multi} multi {ref4multi} --force {params.cell_number} -p {projectname4multi} -a {run_info4multi} {param_test_email} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

rule test_sc_multi_expect_cell:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = os.path.abspath(config["multi"]["lib_csv"]),
        metadata = os.path.abspath(config["multi"]["metadata"]),
    params:
        dir4test = os.path.join(outdir_abspath, config["multi"]["projectname"] + "_multi_expectcell/Analysis/"),
        cell_number = config["multi"]["cell_number"]
    output:
        log = os.path.join(outdir_abspath, "test_sc_multi_expectcell.log")
    shell:
        """
mkdir -p {params.dir4test}
cp {input.metadata} {params.dir4test}/../ 
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4multi} multi {ref4multi} --expect {params.cell_number} -p {projectname4multi} -a {run_info4multi} {param_test_email} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

if "crispr" in config: 
    ## crispr
    run_info4crispr   = config["crispr"]["archive_run"]
    fastqpath4crispr  = config["crispr"]["fastqpath"]
    projectname4crispr = config["crispr"]["projectname"]
    ref4crispr = config["crispr"]["ref"]
    rule test_sc_crispr_default:
        input:
            run_snakemake4sc = config["path2run_snakemake_sc"],
            lib_csv = config["crispr"]["lib_csv"]
        params:
            dir4test = os.path.join(outdir_abspath, config["crispr"]["projectname"] + "_crispr_default/Analysis/"),
        output:
            log = os.path.join(outdir_abspath, "test_sc_crispr_default.log")
        shell:
            """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4crispr} multi {ref4crispr}  -p {projectname4crispr} -a {run_info4crispr} {param_test_email} > {output.log} 2>&1
sed -i 's/""/"libraries.csv"/g' {params.dir4test}/config.py
cp {input.lib_csv} {params.dir4test}/libraries.csv
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun
"""
