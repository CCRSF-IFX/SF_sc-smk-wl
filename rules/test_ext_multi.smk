import os

fastqpath4multi= config["multi"]["fastqpath"].split(",")
genome4multi = config["multi"]["genome"]
ref4multi = config["multi"]["ref"]

rule test_sc_multi_default:
    input:
        run_snakemake4sc = config["path2run_snakemake_sc"],
        lib_csv = os.path.abspath(config["multi"]["lib_csv"]),
    params:
        dir4test = os.path.join(outdir_abspath, config["multi"]["projectname"] + "_multi_default/Analysis/"),
    output:
        log = os.path.join(outdir_abspath, "test_sc_multi_default.log")
    shell:
        """
if [ -d "{params.dir4test}" ]; then
    rm -rf $(readlink -f "{params.dir4test}")
fi
mkdir -p {params.dir4test}
cd {params.dir4test} && echo "n" | python {input.run_snakemake4sc} multi -f {fastqpath4multi} -r {ref4multi} -g {genome4multi} -l {input.lib_csv} > {output.log} 2>&1
cd {params.dir4test} && echo "{submit_job}" | python {input.run_snakemake4sc} rerun 
"""



"""
/mnt/ccrsf-static/Analysis/xies4/github_repos/pipeline_dev_test
usage: scMaestro multi [-h] -f [FASTQPATH [FASTQPATH ...]] -r REFERENCE FOLDER
                       -g GENOME [--chain {auto,TR,IG}] [--fullanalysis]
                       [--vdj_ref VDJ_REF] -l LIBRARY_CONFIG [--cmo] [--count]

optional arguments:
  -h, --help            show this help message and exit
  -f [FASTQPATH [FASTQPATH ...]], --fastqs [FASTQPATH [FASTQPATH ...]]
                        Path(s) to fastq files, multiple paths can be provided
                        together. eg. "-f path1 path2"
  -r REFERENCE FOLDER, --reference REFERENCE FOLDER
                        Reference folder for alignment. VDJ reference is
                        needed if subcommand "vdj" is used.
  -g GENOME, --genome GENOME
                        Genome build, e.g. "hg38", "mm10"
  --chain {auto,TR,IG}  Force the analysis to be carried out for a particular
                        chain type.
  --fullanalysis        Run full analysis pipeline
  --vdj_ref VDJ_REF     Reference folder for VDJ data
  -l LIBRARY_CONFIG, --library_config LIBRARY_CONFIG
                        CSV file with the library configration.
  --cmo                 CMO information will be used for multi analysis
  --count               Run cellranger count for projects with HTO libraries
                        with more 10 individuals mixed.
"""
