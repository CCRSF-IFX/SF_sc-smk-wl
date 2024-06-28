import os

## spatial
print(config)
print(config["spatial"])
run_info4spatial = config["spatial_polyA"]["archive_run"]
fastqpath4spatial=os.path.abspath(config["spatial_polyA"]["fastqpath"])
projectname4spatial = config["spatial_polyA"]["projectname"]
ref4spatial = config["spatial_polyA"]["ref"]
spatial_method = config["spatial_polyA"]["spatial_method"]
imagemeta = os.path.abspath(config["spatial_polyA"]["imagemeta"])

rule test_sc_spatial_polyA_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        metadata = os.path.abspath(config["spatial_polyA"]["metadata"])
    params:
        dir4test = os.path.abspath(os.path.join(outdir_abspath,config["spatial_polyA"]["projectname"]+"_spatial_default/Analysis/")),
    output:
        log = os.path.join(outdir_abspath, "test_spatial_polyA_default.log")
    shell:
        """
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4spatial} spatial {ref4spatial} -p {projectname4spatial} --spatial-method {spatial_method} --images {imagemeta}  -a {run_info4spatial} {param_test_email} > {output.log} 2>&1
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

