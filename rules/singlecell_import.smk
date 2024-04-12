import json
import packaging.version
from xml.dom import minidom
include: "runParametersImport" 
if 'sflog' not in globals():
    import logging as sflog
    sflog.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=sflog.INFO)
else:
    pass
    sflog.debug("sflog has been imported")

# by default, config is an empty dictionary.
# if --configfile is provided, meaning that external users
# are using the workflow
sflog.debug("Is config a dictionary: " + str(isinstance(config, dict)))
sflog.debug("Is config empty: " + str(not config))
if isinstance(config, dict) and not config:
    import shutil
    from pathlib import Path
    ## Path().absolute() return PosixPath object.
    ## PosixPath object cause import issue. So str()
    ## is used here to convert PosixPath to regular path
    sys.path.insert(0, str(Path().absolute())) 
    sflog.debug(str(Path().absolute()))
    import config
    import reference
    import program
    external = False
else:
    sflog.debug("--configfile is used")
    external = True
    pass

cmd_cellranger = program.cellranger if external == False else "cellranger"

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

with open(config.analysis + '/cluster.json') as file:
    clusterConfig = json.load(file)

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
    sflog.info("Samples to be analyzed: " + " ".join(samples))

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
#flowcells = {os.path.basename(i.strip('/')): i for i in config.unaligned}
#cfile = one_up + "/" + project_name+"_"+'_'.join(flowcells)+".count.tar"
report_result = one_up + "/" + project_name + "_" + flowcell + "_Metadata.txt"
wreport_result = one_up + "/" + project_name + "_" + flowcell + ".docx"
xreport_result = one_up + "/" + project_name + "_" + flowcell + ".xlsx"
copy_result = one_up + "/" + project_name + "_" + flowcell + "_copy.txt"

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

flag_create_bam = ""
if is_version_greater_than(program.cellranger, "8.0.0"):
    flag_create_bam = " --create-bam=true "
