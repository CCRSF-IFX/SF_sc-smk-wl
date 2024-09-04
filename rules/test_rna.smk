import os 
## RNA
run_info4rna = config["rna"]["archive_run"]
fastqpath4rna=os.path.abspath(config["rna"]["fastqpath"])
projectname4rna = config["rna"]["projectname"]
ref = config["rna"]["ref"]

rule test_sc_rna_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        metadata = os.path.abspath(config["rna"]["metadata"])
    params:
        dir4test = os.path.join(outdir_abspath,config["rna"]["projectname"]+"_rna_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_scRNA_default.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4rna} rna {ref}  -p {projectname4rna} -a {run_info4rna} {param_test_email} > {output.log} 2>&1
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

rule test_sc_rna_exclude_introns:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        metadata = os.path.abspath(config["rna"]["metadata"])
    params:
        dir4test = os.path.join(outdir_abspath, config["rna"]["projectname"] + "_rna_exclude_introns/Analysis/")
    output:
        log = os.path.join(outdir_abspath, "test_scRNA_exclude_introns.log")
    shell:
        """
mkdir -p {params.dir4test}
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} {fastqpath4rna} rna {ref} --exclude-introns -p {projectname4rna} -a {run_info4rna} {param_test_email} > {output.log} 2>&1 

"""

rule test_sc_rna_force_cell:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        metadata = os.path.abspath(config["rna"]["metadata"])
    params:
        dir4test = os.path.join(outdir_abspath, config["rna"]["projectname"] + "_rna_forcecell/Analysis/"),
        cell_number = config["rna"]["cell_number"]
    output:
        log = os.path.join(outdir_abspath, "test_scRNA_force_cell.log")
    shell:
        """
mkdir -p {params.dir4test}
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} {fastqpath4rna} rna {ref} --force {params.cell_number} -p {projectname4rna} -a {run_info4rna}  {param_test_email} > {output.log} 2>&1
"""

rule test_sc_rna_expect_cell:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"]
    params:
        dir4test = os.path.join(outdir_abspath, config["rna"]["projectname"] + "_rna_expectcell/Analysis/"),
        cell_number = config["rna"]["cell_number"]
    output:
        log = os.path.join(outdir_abspath, "test_scRNA_expect_cell.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} {fastqpath4rna} rna {ref} --expect {params.cell_number} -p {projectname4rna} -a {run_info4rna}  {param_test_email} > {output.log} 2>&1
"""
