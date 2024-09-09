import os

## spatial
#print(config)
print(config["spatial"])
run_info4spatial = config["spatial"]["archive_run"]
fastqpath4spatial=os.path.abspath(config["spatial"]["fastqpath"])
projectname4spatial = config["spatial"]["projectname"]
ref4spatial = config["spatial"]["ref"]
spatial_method = config["spatial"]["spatial_method"]
imagemeta = config["spatial"]["imagemeta"]


rule test_sc_spatial_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        metadata = os.path.abspath(config["spatial"]["metadata"])
    params:
        dir4test = os.path.abspath(os.path.join(outdir_abspath,config["spatial"]["projectname"]+"_spatial_default/Analysis/")),
    output:
        log = os.path.join(outdir_abspath, "test_scRNA_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f {params.dir4test}/../)
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | {input.run_snakemake4sc} {fastqpath4spatial} spatial {ref4spatial} -p {projectname4spatial} --spatial-method {spatial_method} --images {imagemeta}  -a {run_info4spatial} {param_test_email} > {output.log} 2>&1
cp {input.metadata} {params.dir4test}/../
touch archive_setup.complete
cd {params.dir4test} && echo "{submit_job}" | {input.run_snakemake4sc} --rerun 
"""

