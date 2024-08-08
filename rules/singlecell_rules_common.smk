import os
def get_workflow_img(wildcards):
    img_path = ""
    if config.pipeline == "pipseq": 
        img_abs_path = os.path.join(os.getcwd(), "workflow/images/SingleCell_RNA_PIPseq.png")
        img_path = f" -w {img_abs_path}"
    return img_path

testing = getattr(config, "test_email", False)
args4wreport_test = "--test" if testing else ""
yields = getattr(config, "yields", False)
args4wreport_yields = f"--yields {yields}" if yields else ""

software = "Cellranger"
if config.pipeline == "pipseq":
    software = "PIPseeker"

def output_list_web_summary(wildcards):
    """Get the list of web summary from cellranger count or cellranger multi"""
    if config.pipeline == "pipseq":
        return expand("{sample}/barcodes/barcode_whitelist.txt", sample=samples)
    elif config.pipeline == "multi":
        return expand("{sample}/outs/per_sample_outs/{sample}/web_summary.html", sample=samples)
    else:
        return expand("{sample}/outs/web_summary.html", sample=samples)

if config.pipeline == "pipseq" or config.pipeline == "nopipe":
    aggregate = False

active_script_folder = ""
if hasattr(program, "active_scripts"):
    active_script_folder = program.active_scripts

if external == False:
    if testing:
        rule copy:
            input: xreport_result, wreport_result#, expand("finalreport/{sample}.html", sample=samples)
            output: touch('copy.complete')
            log: copy_result
            shell: "mkdir -p {program.copydir}/pipeline_testing/{run_name}/{project_name}; cp -r finalreport {program.copydir}/pipeline_testing/{run_name}/{project_name} >> {log}; cd {one_up}; cp -v *.docx {program.copydir}/pipeline_testing/{run_name}/{project_name} >> {log}; cp -v *.xlsx {program.copydir}/pipeline_testing/{run_name}/{project_name} >> {log}"
    else:
        rule copy:
            input: xreport_result, wreport_result#, expand("finalreport/{sample}.html", sample=samples)
            output: touch('copy.complete')
            log: copy_result
            shell: "mkdir -p {program.copydir}/{run_name}/{project_name}; cp -r finalreport {program.copydir}/{run_name}/{project_name} >> {log}; cd {one_up}; cp -v *.docx {program.copydir}/{run_name}/{project_name} >> {log}; cp -v *.xlsx {program.copydir}/{run_name}/{project_name} >> {log}"
    
    if aggregate:
        sflog.info(aggregate)
        rule archive:
            input: metadata = report_result, aggr_log = "run_10x_aggregate.log"
            output: touch('archive_setup.complete')
            params: fastqs = ",".join([os.path.dirname(name.rstrip('/')) for name in flowcells.values()]), runs = ','.join([j for i in flowcells for j in run_names if i in j])
            log: "archive.log"
            shell: "cd {one_up}; python {analysis}/workflow/scripts/meta2json_single_cell_v0.1.py -m {input.metadata} -r {params.runs} -f {params.fastqs} -c {config.analysis} -a {config.analysis}/AggregatedDatasets > {log}"
    else:   ## aggregate is false 
        rule archive:
            input:  metadata = report_result, summaryFiles = "finalreport/metric_summary.xlsx"
            output: touch('archive_setup.complete')
            params: fastqs = ",".join([os.path.dirname(name.rstrip('/')) for name in flowcells.values()]), runs = ','.join([j for i in flowcells for j in run_names if i in j])
            log: "archive.log"
            shell: "cd {one_up}; python {analysis}/workflow/scripts/meta2json_single_cell_v0.1.py -m {input.metadata} -r {params.runs} -f {params.fastqs} -c {config.analysis} > {log}"

    rule report:
        output: 
            report_result
        shell: 
            """
cd {one_up}; perl {active_script_folder}/run_GenerateAllReports.pl -flowcell {flowcell} -p {csas} -s 2 -u {unaligned} -f {run_name} -t 0 -m 1 -e 2 -r 1 -x verajc:test1test2 -software {software} ; rm -f {xreport_result}
"""

    rule xreport:
        input: report_result, metric="finalreport/metric_summary.xlsx"
        output: xreport_result
        shell: "cp {input.metric} {xreport_result}"

    rule wreport:
        input: metadata = report_result, excel = "finalreport/metric_summary.xlsx"
        output: wreport_result
        params: runs = ','.join(run_names), pipeline = config.pipeline, workflow_img = get_workflow_img
        shell: "cd {one_up}; python {analysis}/workflow/scripts/SF_scWordReport/run_wordreport_sc.py -e {analysis}/{input.excel} -m {input.metadata} -c {current_cellranger} -p {params.pipeline} -r {params.runs} {params.workflow_img} {args4wreport_test} {args4wreport_yields}"

    if testing: 
        onsuccess:
            success = "Yes"
            al = analysis
            shell(r'''python {analysis}/workflow/scripts/sendmail4test.py {success} {al} {testing} > snakemake.stats''')
    
        onerror:
            success = "No"
            al = analysis
            shell(r'''python {analysis}/workflow/scripts/sendmail4test.py {success} {al} {testing}''')
    else:
        onsuccess:
            success = "Yes"
            al = analysis
            shell(r'''python {analysis}/workflow/scripts/sendmail.py {success} {al} > snakemake.stats''')
    
        onerror:
            success = "No"
            al = analysis
            shell(r'''python {analysis}/workflow/scripts/sendmail.py {success} {al}''') 
