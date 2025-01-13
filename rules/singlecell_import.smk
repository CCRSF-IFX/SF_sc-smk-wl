import json
import packaging.version
from xml.dom import minidom
import shutil
from pathlib import Path
import pandas as pd
import xml.etree.ElementTree as ET


include: "runParametersImport" 
<<<<<<< HEAD

## Path().absolute() return PosixPath object.
## PosixPath object cause import issue. So str()
## is used here to convert PosixPath to regular path
sys.path.insert(0, str(Path().absolute()))
import config
import reference
import program

if 'sflog' not in globals():
    import logging as sflog
    debug_bool = getattr(config, "debug", False)  
    if debug_bool == True:
        sflog.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=sflog.DEBUG)
    else: 
        sflog.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=sflog.INFO)
else:
    pass
    sflog.debug("sflog has been imported")
=======
>>>>>>> main

sflog.debug(str(Path().absolute()))

# by default, config is an empty dictionary.
# if --configfile is provided, meaning that external users
# are using the workflow
#print("Is config a dictionary: " + str(isinstance(config, dict)))
#print("Is config empty: " + str(not config))

<<<<<<< HEAD
=======
## Path().absolute() return PosixPath object.
## PosixPath object cause import issue. So str()
## is used here to convert PosixPath to regular path
sys.path.insert(0, str(Path().absolute())) 
#print(str(Path().absolute()))
import config
import reference
import program
>>>>>>> main

def get_bool4internal(): 
    program_filep = os.path.join(config.analysis, "program.py")
    try:
        with open(program_filep) as progf:
            for line in progf:
                if "copydir = " in line or "active_scripts =" in line:
                   return False 
        return True
    except FileNotFoundError:
        #print(f"{program_filep} not found")
        return False

external = get_bool4internal()
<<<<<<< HEAD
sflog.debug(f"Workflow is used by external user: {external}.")
=======
#print(f"Workflow is used by external user: {external}.")
>>>>>>> main

#container: program.global_container

localrules: count

csas = re.search("CS[0-9]{6}", config.analysis).group(0) if re.search("CS[0-9]{6}", config.analysis) else os.path.basename(config.analysis.strip('/'))
unaligned = config.unaligned[0]

#Get the run names
if hasattr(config, "runs"):
    run_names = config.runs.split(',')
else:
    run_names = list()
    for unaligned in config.unaligned:
      if unaligned.split('outs')[0].split('/')[-3] == 'Unaligned':
        run_names.append(unaligned.split('outs')[0].split('/')[-4])
      else:
        run_names.append(unaligned.split('outs')[0].split('/')[-3])

#Try to get the most recent run, used in copy rule
run_names = list(set(run_names))
run_names.sort()
run_name = run_names[-1]

fastqpath = config.unaligned
analysis = config.analysis
one_up = '/'.join(config.analysis.rstrip('/').split('/')[:-1])

#with open(config.analysis + '/cluster.json') as file:
#    clusterConfig = json.load(file)

forcecells = getattr(config, "forcecells", False)

#Get project name
if hasattr(config, "projectname"):
    project_name = config.projectname
else:
    project_name = os.path.basename(config.analysis.strip('/'))

#Get sample names
if hasattr(config, "samples"):
    samples = config.samples
else:
    sample = list(set([os.path.basename(file).split('.')[0] for file in list(itertools.chain.from_iterable([glob.glob(i + '/*') for i in fastqpath]))]))
    samps = []
    for item in sample:
        if len(re.findall(r"(\S*)_S\d+_L0\d{2}_[RI]", item)) > 0:
          samps.append(re.findall(r"(\S*)_S\d+_L0\d{2}_[RI]", item)[0])
        else:
          samps.append(item)
    ## Add if else in the list comprehension for the sample name with 'Sample_' in the middle.
    ## For example, Tube_1___Sample_3__GEX_library in CS033737 
    samples = list(set(s.replace('Sample_', '') if s.startswith('Sample_') else s for s in set(samps)))
    samples = sorted(samples)

# Move fastq files to subfolders with the names of the sample.
from subprocess import Popen, PIPE
reformat_paths = [i for i in fastqpath if len(glob.glob(os.path.join(i, "*", ""))) == 0]
for i in reformat_paths:
    for sample in [j for j in samples if len(glob.glob(os.path.join(i, j + "*_S[0-9]*_L0[0-9]*_*fastq.gz"))) > 0]:
        process = Popen("mkdir %s/%s; mv %s/%s_S[0-9]*_L0[0-9]*_*fastq.gz %s/%s" % (i, sample, i, sample, i, sample), shell=True, stdout=PIPE, stderr=PIPE)
        process.wait()

#For deliverFastq
def flowcellPath(wildcards):
    return flowcells[wildcards.flowcell]

#Function used to get the fastq paths for each sample
def filterFastq(wildcards):
    paths = [i + "/%s" % wildcards.sample for i in fastqpath]
    paths += [i + "/Sample_%s" % wildcards.sample for i in fastqpath]
    return(','.join([j for j in paths if os.path.exists(j)] + [os.path.dirname(j) for j in paths if not os.path.exists(j) if len(glob.glob(j + "*")) > 0]))

def getFirstFastqFile(wildcards):
    path = filterFastq(wildcards).split(',')[0]
    suffix = '_001.fastq.gz'
    return({i: glob.glob(f'{path}/*_{i}{suffix}')[0] for i in ['R1', 'R2', 'R3'] if len(glob.glob(f'{path}/*_{i}{suffix}')) > 0})

def filterFastq4pipseeker(wildcards):
    """
    Prepare the folders for pipseeker 
    """
    path_fq_new = "fastq/Sample_%s/" % wildcards.sample
    process = Popen("mkdir -p %s" % (path_fq_new), shell=True, stdout=PIPE, stderr=PIPE)
    process.wait()
    cnt_fq_file = 0
    for index, fq_path in enumerate(fastqpath):
        path_sample = os.path.join(fq_path, "Sample_%s" % wildcards.sample) 
        for fastq_file in glob.glob(os.path.join(path_sample,  "*fastq.gz")):
            cnt_fq_file = cnt_fq_file + 1 
            basename_fastq = os.path.basename(fastq_file)
            basename_fastq_new = "%s_%s" % (run_names[index], basename_fastq) 
            cmd = "ln -s %s %s%s" % (fastq_file, path_fq_new, basename_fastq_new) 
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            process.wait()
    if cnt_fq_file == 0: 
        sys.stderr.write("\nNo fastq detected. Please check it out! \n\n")
        sys.exit()
    return(path_fq_new)

def filterFastq4nopipe(wildcards):
    """
    Prepare the folders for nopipe
    """
    path_fq_new = "fastq/%s/" % wildcards.sample
    process = Popen("mkdir -p %s" % (path_fq_new), shell=True, stdout=PIPE, stderr=PIPE)
    process.wait()
    cnt_fq_file = 0
    for index, fq_path in enumerate(fastqpath):
        path_sample = os.path.join(fq_path, "%s" % wildcards.sample)
        for fastq_file in glob.glob(os.path.join(path_sample,  "*fastq.gz")):
            cnt_fq_file = cnt_fq_file + 1
            basename_fastq = os.path.basename(fastq_file)
            basename_fastq_new = "%s_%s" % (run_names[index], basename_fastq)
            cmd = "ln -s %s %s%s" % (fastq_file, path_fq_new, basename_fastq_new)
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            process.wait()
    if cnt_fq_file == 0:
        sys.stderr.write("\nNo fastq detected. Please check it out! \n\n")
        sys.exit()
    return(path_fq_new)

#Setting aggregate flag, this gets turned off for certain pipelines
aggregate = True

#include: program.runParametersImport
flowcell = run_names[-1][-9:]
flowcells = dict()
for run_name in run_names:
    try:
        (RTAVersion, flowcellRunParameters, workFlowType, flowcellMode, chemistry, chemistryVersion, xmlRunParametersPath, xmlRunInfoPath) = runParametersXmlPath(run_name)
    except IOError as error:
        xmlRunParametersPath = "Unknown"
        xmlRunInfoPath = "Unknown"
        sys.stdout.write("\n\n\nExecution failed: " + str(error) +"\n")
        sys.stdout.write("No RunParameters.xml and RunInfo.xml will be archied and flowcell ID is geneated from the analysis folder\n\n\n")
    except:
        sys.stdout.write("Unexpected error:" + str(sys.exc_info()[0]) +"\n")
    if flowcellRunParameters == "Unknown":
        flowcellRunParameters = run_name[-9:]
    for path in fastqpath:
        if flowcellRunParameters in path:
            flowcells[flowcellRunParameters] = path
    flowcell = flowcellRunParameters if flowcell in flowcellRunParameters else flowcell

#Create file names
#flowcell = os.path.basename(config.unaligned[0].strip('/'))
<<<<<<< HEAD
#flowcells = {os.path.basename(i.strip('/')): i for i in config.unaligned}
=======


#flowcells = {os.path.basename(i.strip('/')): i for i in config.unaligned}

def get_flowcell_name_from_reports(path):
    """
    Tries to retrieve the flowcell name from the Reports/html directory for a given path.

    Args:
        path (str): The directory path to check for flowcell in Reports/html.
    
    Returns:
        str: The flowcell name if found, otherwise None.
    """
    ## bcl2fastq output folder
    reports_html_path = os.path.join(path, "../Reports/html")
    if os.path.exists(reports_html_path):
        try:
            # Find a folder in html directory that matches the flowcell naming pattern
            flowcell_name = next(
                d for d in os.listdir(reports_html_path) 
                    if os.path.isdir(os.path.join(reports_html_path, d))
            )
            return flowcell_name
        except StopIteration:
            return None  # No valid flowcell name found
    ## bcl converter 
    runinf_xml = os.path.join(path, "../Reports/RunInfo.xml") 
    if os.path.exists(runinf_xml):
        try:
            # Find a folder in html directory that matches the flowcell naming pattern
            tree = ET.parse(runinf_xml)
            root = tree.getroot()
            # Find the Flowcell element
            flowcell = root.find(".//Flowcell")
            return flowcell.text if flowcell is not None else None
        except ET.ParseError as e:
            print(f"Error parsing XML: {e}")
            return None  # No valid flowcell name found
    return None

# Combine with your existing logic for `config.unaligned`
flowcells = {}
for path in config.unaligned:
    base_name = os.path.basename(path)
    flowcell_name = base_name  # Default to the base name

    # Check if the flowcell name needs to be extracted from Reports/html
    alternative_name = get_flowcell_name_from_reports(path)
    if alternative_name:
        flowcell_name = alternative_name
    
    flowcells[flowcell_name] = path

>>>>>>> main
#cfile = one_up + "/" + project_name+"_"+'_'.join(flowcells)+".count.tar"
report_result = one_up + "/" + project_name + "_" + flowcell + "_Metadata.txt"
#print(one_up)
wreport_result = one_up + "/" + project_name + "_" + flowcell + ".docx"
xreport_result = one_up + "/" + project_name + "_" + flowcell + ".xlsx"
copy_result = one_up + "/" + project_name + "_" + flowcell + "_copy.txt"

#print(flowcells)

rule_all_append = []
if hasattr(config, 'archive'):
    if config.archive:
        archive = True
        rule_all_append += ["archive_setup.complete"]

def is_version_greater_than(version_string, target_version):
    """Checks if a version string is greater than a specific version.

    Args:
        version_string: The version string to check.
        target_version: The version to compare

    Returns:
        True if the version is greater than 8.0.0, False otherwise.
    """
    try:
        version = packaging.version.parse(version_string.split(":")[-1])
        return version >= packaging.version.parse(target_version)
    except packaging.version.InvalidVersion:
        return False  # Handle cases where the version string is invalid

flag4cellranger_create_bam = ""
if is_version_greater_than(program.cellranger, "8.0.0"):
    flag4cellranger_create_bam = " --create-bam=true "

flag4spaceranger_create_bam = ""
if is_version_greater_than(program.spaceranger, "3.0.0"):
    flag4spaceranger_create_bam = " --create-bam=true "

def get_reference_transcriptome(wildcards):
    ref_path = reference.transcriptome
    if hasattr(config, 'libraries'):
        lib_df = pd.read_csv(config.libraries)
        if "transcriptome" in lib_df.columns:
            tem_ref = lib_df.loc[lib_df['Sample'] == wildcards.sample, 'transcriptome']
            if not tem_ref.empty:
                ref_path = tem_ref.iloc[0]
            else:
                sys.stderr.write(f"\n'No information found in 'transcriptome' column of 'libraries.csv' for sample '{wildcards.sample}'\n\n")
                sys.exit()
            if not hasattr(config, 'aggregate') or config.aggregate == True:
                sys.stderr.write("\n'transcriptome' column in 'libraries.csv' detected, indicating more than one reference exist. Please set `aggregate` as `False`\n\n")
                sys.exit()
    return ref_path

def generate_option_flag(wildcards, lib_df, column_name, optional_flag):
    if column_name in lib_df.columns:
        tem_vals = lib_df.loc[lib_df['Sample'] == wildcards.sample, column_name]
        if not tem_vals.empty:
            tem_val = tem_vals.iloc[0]
            optional_flag = f'{optional_flag} --{column_name}={tem_val}'
        else:
            sys.stderr.write(f"\n'No information found in '{column_name}' column of 'libraries.csv' for sample '{wildcards.sample}'\n\n")
            sys.exit()
    return optional_flag

def get_optional_flag_from_lib_csv(wildcards):
    optional_flag = ''
    flags = ['force-cells', 'expect-cells', 'chemistry']
    if hasattr(config, 'libraries'):
        lib_df = pd.read_csv(config.libraries)
        for flag in flags:
            optional_flag = generate_option_flag(wildcards, lib_df, flag, optional_flag)
    return optional_flag
