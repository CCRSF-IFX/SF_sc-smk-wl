import os

## spatial
#print(config)
#print(config["spatial_polyA"])
run_info4spatial_polyA = config["spatial_polyA"]["archive_run"]
fastqpath4spatial_polyA = os.path.abspath(config["spatial_polyA"]["fastqpath"])
projectname4spatial_polyA = config["spatial_polyA"]["projectname"]
ref4spatial_polyA = config["spatial_polyA"]["ref"]
spatial_method_polyA = config["spatial_polyA"]["spatial_method"]
imagemeta_polyA = os.path.abspath(config["spatial_polyA"]["imagemeta"])

rule test_sc_spatial_polyA_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        metadata = os.path.abspath(config["spatial_polyA"]["metadata"])
    params:
        dir4test = os.path.abspath(os.path.join(outdir_abspath, config["spatial_polyA"]["projectname"]+"_spatial_default/Analysis/")),
    output:
        log = os.path.join(outdir_abspath, "test_spatial_polyA_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f {params.dir4test}/../)
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4spatial_polyA} spatial {ref4spatial_polyA} -p {projectname4spatial_polyA} --spatial-method {spatial_method_polyA} --images {imagemeta_polyA}  -a {run_info4spatial_polyA} {param_test_email} > {output.log} 2>&1
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

