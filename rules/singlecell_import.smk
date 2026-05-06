import json
import packaging.version
from xml.dom import minidom
import shutil
from pathlib import Path
import pandas as pd
import xml.etree.ElementTree as ET
import glob
import itertools
import os
import re
import sys
from collections import Counter, defaultdict
import uuid

include: "runParametersImport" 

def create_pipeline_id_file():
    """
    Creates a hidden file with a unique pipeline ID (e.g.: 550e8400-e29b-41d4-a716-446655440000)
    The ID is generated using UUID for uniqueness.
    Returns the generated pipeline ID.
    """
    # Generate unique pipeline ID using UUID
    pipeline_id = str(uuid.uuid4())
    
    # Create hidden file with the pipeline ID
    hidden_file_path = os.path.join(analysis, ".pipeline_id")
    
    try:
        with open(hidden_file_path, 'w') as f:
            f.write(pipeline_id)
        print(f"Pipeline ID created: {pipeline_id}")
        return pipeline_id
    except Exception as e:
        sys.stderr.write(f"Error creating pipeline ID file: {e}\n")
        return None

def get_pipeline_id():
    """
    Retrieves the pipeline ID from the hidden file if it exists,
    otherwise creates a new one.
    """
    hidden_file_path = os.path.join(analysis, ".pipeline_id")
    
    if os.path.exists(hidden_file_path):
        try:
            with open(hidden_file_path, 'r') as f:
                pipeline_id = f.read().strip()
            return pipeline_id
        except Exception as e:
            sys.stderr.write(f"Error reading pipeline ID file: {e}\n")
            return create_pipeline_id_file()
    else:
        return create_pipeline_id_file()

# by default, config is an empty dictionary.
# if --configfile is provided, meaning that external users
# are using the workflow
#print("Is config a dictionary: " + str(isinstance(config, dict)))
#print("Is config empty: " + str(not config))

## Path().absolute() return PosixPath object.
## PosixPath object cause import issue. So str()
## is used here to convert PosixPath to regular path
sys.path.insert(0, str(Path().absolute())) 
#print(str(Path().absolute()))
import config
import reference
import program

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
#print(f"Workflow is used by external user: {external}.")

if config.pipeline == "curioseeker":
    localrules: count

csas = re.search("CS[0-9]{6}", config.analysis).group(0) if re.search("CS[0-9]{6}", config.analysis) else os.path.basename(config.analysis.strip('/'))
unaligned = config.unaligned[0]

def infer_run_name_from_unaligned_path(unaligned_path):
    parts = Path(unaligned_path.rstrip("/")).parts
    if "Unaligned" in parts:
        unaligned_index = len(parts) - 1 - parts[::-1].index("Unaligned")
        if unaligned_index > 0:
            return parts[unaligned_index - 1]
    if "outs" in parts:
        outs_index = len(parts) - 1 - parts[::-1].index("outs")
        if outs_index > 0:
            return parts[outs_index - 1]
    if len(parts) >= 3:
        return parts[-3]
    return os.path.basename(unaligned_path.rstrip("/"))

def make_unique_run_labels(run_names_for_paths):
    totals = Counter(run_names_for_paths)
    seen = defaultdict(int)
    labels = []
    for index, run_name_for_path in enumerate(run_names_for_paths):
        base_label = run_name_for_path or f"fastq_root_{index + 1}"
        seen[base_label] += 1
        if totals[base_label] > 1:
            labels.append(f"{base_label}_root{seen[base_label]}")
        else:
            labels.append(base_label)
    return labels

#Get the run names
if hasattr(config, "runs"):
    configured_run_names = [run.strip() for run in config.runs.split(',') if run.strip()]
else:
    configured_run_names = [infer_run_name_from_unaligned_path(unaligned) for unaligned in config.unaligned]

run_names_orig = [
    configured_run_names[index] if index < len(configured_run_names) else infer_run_name_from_unaligned_path(unaligned)
    for index, unaligned in enumerate(config.unaligned)
]
run_labels_by_fastqpath = make_unique_run_labels(run_names_orig)

#Try to get the most recent run, used in copy rule
run_names = list(dict.fromkeys(run_names_orig))  # Removes duplicates while preserving order
run_names.sort()
run_name = run_names[-1]

fastqpath = config.unaligned
analysis = config.analysis
one_up = '/'.join(config.analysis.rstrip('/').split('/')[:-1])
FASTQ_MANIFEST_PATH = os.path.join("fastq", "fastq_manifest.tsv")

# Initialize pipeline unique ID
pipeline_id = get_pipeline_id()

#with open(config.analysis + '/cluster.json') as file:
#    clusterConfig = json.load(file)

forcecells = getattr(config, "forcecells", False)

#Get project name
if hasattr(config, "projectname"):
    project_name = config.projectname
else:
    project_name = os.path.basename(config.analysis.strip('/'))

def sample_name_from_fastq_entry(entry):
    basename = os.path.basename(entry).split('.')[0]
    match = re.match(r"(.+)_S\d+(?:_L\d{3})?_[RI]\d?_001$", basename)
    return match.group(1) if match else basename

def parse_fastq_filename(fastq_file, fallback_sample=None):
    basename = os.path.basename(fastq_file)
    match = re.match(r"(.+)_S(\d+)(?:_L(\d{3}))?_([RI]\d?)_001\.fastq\.gz$", basename)
    if match:
        sample, sample_number, lane, read = match.groups()
    else:
        sample = fallback_sample if fallback_sample else basename.replace(".fastq.gz", "")
        sample_number = ""
        lane = ""
        read = ""
    if sample.startswith("Sample_"):
        sample = sample.replace("Sample_", "", 1)
    return {
        "sample": sample,
        "sample_number": sample_number or "",
        "lane": lane or "",
        "read": read or "",
    }

def detect_sample_folder(fq_path, sample):
    for sample_folder in [f"Sample_{sample}", sample]:
        if os.path.exists(os.path.join(fq_path, sample_folder)):
            return sample_folder
    return None

def direct_fastq_files_for_sample(fq_path, sample):
    fastq_files = set()
    for sample_prefix in [sample, f"Sample_{sample}"]:
        fastq_files.update(glob.glob(os.path.join(fq_path, f"{sample_prefix}_S[0-9]*_L[0-9][0-9][0-9]_*fastq.gz")))
        fastq_files.update(glob.glob(os.path.join(fq_path, f"{sample_prefix}_S[0-9]*_*fastq.gz")))
    return sorted(fastq_files)

def fastq_files_for_sample(fq_path, sample, sample_folder=None):
    if sample_folder:
        return sorted(glob.glob(os.path.join(fq_path, sample_folder, "*fastq.gz")))
    return direct_fastq_files_for_sample(fq_path, sample)

def run_name_for_fastq_path(index, fq_path):
    if index < len(run_names_orig) and run_names_orig[index]:
        return run_names_orig[index]
    return infer_run_name_from_unaligned_path(fq_path)

def run_label_for_fastq_path(index, fq_path):
    if index < len(run_labels_by_fastqpath) and run_labels_by_fastqpath[index]:
        return run_labels_by_fastqpath[index]
    run_name = run_name_for_fastq_path(index, fq_path)
    return run_name or f"fastq_root_{index + 1}"

def build_fastq_manifest_rows():
    rows = []
    seen = set()
    for sample in samples:
        for index, fq_path in enumerate(fastqpath):
            run_name = run_name_for_fastq_path(index, fq_path)
            run_label = run_label_for_fastq_path(index, fq_path)
            detected_sample_folder = detect_sample_folder(fq_path, sample)
            source_layout = "sample_folder" if detected_sample_folder else "loose"
            for fastq_file in fastq_files_for_sample(fq_path, sample, detected_sample_folder):
                source_path = os.path.abspath(fastq_file)
                if source_path in seen:
                    continue
                seen.add(source_path)
                parsed = parse_fastq_filename(fastq_file, fallback_sample=sample)
                staged_name = f"{run_label}_{os.path.basename(fastq_file)}" if run_label else os.path.basename(fastq_file)
                rows.append({
                    "sample": sample,
                    "read": parsed["read"],
                    "lane": parsed["lane"],
                    "sample_number": parsed["sample_number"],
                    "run_name": run_name,
                    "run_label": run_label,
                    "fastq_root_index": index + 1,
                    "fastq_root": os.path.abspath(fq_path),
                    "source_layout": source_layout,
                    "source_sample_folder": detected_sample_folder or "",
                    "source_path": source_path,
                    "staged_path": os.path.abspath(os.path.join(analysis, "fastq", sample, staged_name)),
                })
    return sorted(rows, key=lambda row: (row["sample"], row["run_name"], row["source_path"]))

def write_fastq_manifest(manifest_path):
    columns = [
        "sample",
        "read",
        "lane",
        "sample_number",
        "run_name",
        "run_label",
        "fastq_root_index",
        "fastq_root",
        "source_layout",
        "source_sample_folder",
        "source_path",
        "staged_path",
    ]
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
    with open(manifest_path, "w") as manifest:
        manifest.write("\t".join(columns) + "\n")
        for row in build_fastq_manifest_rows():
            manifest.write("\t".join(str(row[column]) for column in columns) + "\n")

#Get sample names
if hasattr(config, "samples"):
    samples = config.samples
else:
    sample = list(set([os.path.basename(file).split('.')[0] for file in list(itertools.chain.from_iterable([glob.glob(i + '/*') for i in fastqpath]))]))
    samps = []
    for item in sample:
        samps.append(sample_name_from_fastq_entry(item))
    ## Add if else in the list comprehension for the sample name with 'Sample_' in the middle.
    ## For example, Tube_1___Sample_3__GEX_library in CS033737 
    samples = list(set(s.replace('Sample_', '') if s.startswith('Sample_') else s for s in set(samps)))
    samples = sorted(samples)

from subprocess import Popen, PIPE

#For deliverFastq
def flowcellPath(wildcards):
    return flowcells[wildcards.flowcell]

#Function used to get the fastq paths for each sample
def filterFastq(wildcards):
    paths = [i + "/%s" % wildcards.sample for i in fastqpath]
    paths += [i + "/Sample_%s" % wildcards.sample for i in fastqpath]
    return(','.join([j for j in paths if os.path.exists(j)] + [os.path.dirname(j) for j in paths if not os.path.exists(j) if len(glob.glob(j + "*")) > 0]))

def getFirstFastqFile(wildcards):
    suffix = '_001.fastq.gz'
    fastq_files = []
    for fq_path in fastqpath:
        detected_sample_folder = detect_sample_folder(fq_path, wildcards.sample)
        fastq_files.extend(fastq_files_for_sample(fq_path, wildcards.sample, detected_sample_folder))
    first_fastq_files = {}
    for read in ['R1', 'R2', 'R3']:
        read_fastq_files = sorted([fastq_file for fastq_file in fastq_files if fastq_file.endswith(f'_{read}{suffix}')])
        if read_fastq_files:
            first_fastq_files[read] = read_fastq_files[0]
    return first_fastq_files

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
        run_label = run_label_for_fastq_path(index, fq_path)
        for fastq_file in glob.glob(os.path.join(path_sample,  "*fastq.gz")):
            cnt_fq_file = cnt_fq_file + 1 
            basename_fastq = os.path.basename(fastq_file)
            basename_fastq_new = "%s_%s" % (run_label, basename_fastq) 
            symlink_path = os.path.join(path_fq_new, basename_fastq_new)
            if os.path.islink(symlink_path):
                os.unlink(symlink_path)  # Unlink the existing symbolic link
            cmd = "ln -s %s %s" % (fastq_file, symlink_path) 
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            process.wait()
    if cnt_fq_file == 0: 
        sys.stderr.write("\nNo fastq detected. Please check it out! \n\n")
        sys.exit()
    return(path_fq_new)

def filterFastq4nopipe(wildcards):
    """
    Prepare the folders for pipseeker or nopipe.
    Automatically detects whether "Sample_" prefix is used in folder structure.
    """
    has_fastqs = any(
        detect_sample_folder(fq_path, wildcards.sample) or direct_fastq_files_for_sample(fq_path, wildcards.sample)
        for fq_path in fastqpath
    )
    if not has_fastqs:
        sys.stderr.write(f"\nError: No FASTQ folder found for sample {wildcards.sample}. Check the directory structure.\n\n")
        sys.exit(1)

    # Stage symlinks under the workflow sample id and keep source-folder detection separate.
    path_fq_new = f"fastq/{wildcards.sample}/"

    # Create directory if it doesn't exist
    os.makedirs(path_fq_new, exist_ok=True)
    # Remove existing files in the directory before creating symlinks
    # This is to ensure that the directory is clean before adding new symlinks
    if os.path.exists(path_fq_new):
        for file in os.listdir(path_fq_new):
            file_path = os.path.join(path_fq_new, file)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)  # Remove file or symbolic link
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)  # Remove directory

    cnt_fq_file = 0
    for index, fq_path in enumerate(fastqpath):
        detected_sample_folder = detect_sample_folder(fq_path, wildcards.sample)
        run_label = run_label_for_fastq_path(index, fq_path)
        for fastq_file in fastq_files_for_sample(fq_path, wildcards.sample, detected_sample_folder):
            cnt_fq_file += 1
            basename_fastq = os.path.basename(fastq_file)
            basename_fastq_new = f"{run_label}_{basename_fastq}"
            symlink_path = os.path.join(path_fq_new, basename_fastq_new)
            os.symlink(fastq_file, symlink_path)

    if cnt_fq_file == 0:
        sys.stderr.write("\nNo FASTQ files detected. Please check it out!\n\n")
        sys.exit(1)

    return path_fq_new

record_fastqpath = {}
record_fastqfiles = defaultdict(set)
def prep_fastq_folder_ln(sample, get_dict_only=False):
    """
    Prepare the folders for fastq symlinks and record the paths and files used.
    Updates the global record_fastqpath and record_fastqfiles dictionaries only if get_dict_only is False.
    """
    global record_fastqpath, record_fastqfiles
    has_fastqs = any(
        detect_sample_folder(fq_path, sample) or direct_fastq_files_for_sample(fq_path, sample)
        for fq_path in fastqpath
    )
    if not has_fastqs:
        sys.stderr.write(f"\nError: No FASTQ folder found for sample {sample}. Check the directory structure.\n\n")
        sys.exit(1)

    # Stage symlinks under the workflow sample id and keep source-folder detection separate.
    path_fq_new = os.path.join(analysis, f"fastq/{sample}/")
    if not get_dict_only:
        if os.path.exists(path_fq_new):
            for file in os.listdir(path_fq_new):
                file_path = os.path.join(path_fq_new, file)
                if os.path.islink(file_path) or os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
        else:
            os.makedirs(path_fq_new, exist_ok=True)

    cnt_fq_file = 0
    temp_fastqfiles = set()
    for index, fq_path in enumerate(fastqpath):
        detected_sample_folder = detect_sample_folder(fq_path, sample)
        run_label = run_label_for_fastq_path(index, fq_path)
        for fastq_file in fastq_files_for_sample(fq_path, sample, detected_sample_folder):
            cnt_fq_file += 1
            basename_fastq = os.path.basename(fastq_file)
            basename_fastq_new = f"{run_label}_{basename_fastq}"
            symlink_path = os.path.join(path_fq_new, basename_fastq_new)
            if not get_dict_only:
                os.symlink(fastq_file, symlink_path)
            temp_fastqfiles.add(symlink_path)

    if cnt_fq_file == 0:
        sys.stderr.write("\nNo FASTQ files detected. Please check it out!\n\n")
        sys.exit(1)

    if get_dict_only:
        record_fastqpath[sample] = path_fq_new
        record_fastqfiles[sample] = temp_fastqfiles

    return path_fq_new

for sample in samples:
    prep_fastq_folder_ln(sample, get_dict_only=True)  # Prepare the fastq folder for each sample

rule fastq_manifest:
    output:
        FASTQ_MANIFEST_PATH
    run:
        write_fastq_manifest(output[0])

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

#cfile = one_up + "/" + project_name+"_"+'_'.join(flowcells)+".count.tar"
report_result = one_up + "/" + project_name + "_" + flowcell + "_Metadata.txt"
wreport_result = one_up + "/" + project_name + "_" + flowcell + ".docx"
xreport_result = one_up + "/" + project_name + "_" + flowcell + ".xlsx"
copy_result = one_up + "/" + project_name + "_" + flowcell + "_copy.txt"

#print(flowcells)

rule_all_append = [FASTQ_MANIFEST_PATH]
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

def process_config_attr(config, attr, sample, flags):
    """
    Processes a given attribute in the config object and updates the flag list.
    Exits the program if any errors are encountered.

    :param config: Configuration object with an optional attribute.
    :param attr: The attribute name to process (e.g., 'cmo', 'reference').
    :param sample: Sample key used if the attribute is a dictionary.
    :param flag: List to append the processed flag.
    """
    if not hasattr(config, attr):
        return

    attr_value = getattr(config, attr)
    if isinstance(attr_value, str):
        flags.append(f"--{attr} {os.path.abspath(attr_value)}")
        if not os.path.isfile(attr_value):
            sys.exit(f"Error: file {attr_value} not found. Please fix the issue.")
    
    elif isinstance(attr_value, dict):
        if sample in attr_value:
            file_path = attr_value[sample]
            flags.append(f"--{attr} {os.path.abspath(file_path)}")
            if not os.path.isfile(file_path):
                sys.exit(f"Error: file {file_path} not found. Please fix the issue.")
        else:
            sys.exit(f"Error: dictionary detected for '{attr}' but key '{sample}' is missing. Please fix the issue.")
 
    else:
        sys.exit(f"Error: '{attr}' should be either string or a dictionary. Please fix the issue.")

    ## https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi#hashing
    ## https://tracker.nci.nih.gov/browse/CCRIFXDEV-33
    ## https://github.com/CCRSF-IFX/SF_sc-smk-wl/blob/fdf45ad657e3267860d75e621f203c8fb40de76b/rules/singlecell_import.smk#L464
    if attr == "cmo": 
        if not hasattr(config, "hashedabc"):
            config.hashedabc = False
        if config.hashedabc == True:
            flags.append("--hashedabc")
    ## In the future release, `hashedabc` should be removed. 
    if attr == "cmo": 
        if not hasattr(config, "hashing_with_abc"):
            config.hashing_with_abc = False
        if config.hashing_with_abc == True:
            flags.append("--hashing_with_abc")

def get_parsebio_sample_names(sample_list_path):
    """
    Parse a list of ParseBio sample names and construct their corresponding input/output directories.

    This function reads a tab-delimited text file containing sample identifiers (one per line)
    and returns a dictionary mapping each sample name to its associated directories used
    in downstream processing (e.g., DGE-filtered matrices and Seurat outputs).

    Lines that are blank or begin with '//' are ignored.

    Parameters
    ----------
    sample_list_path : str
        Path to a text file listing sample names. Each line should contain at least
        one tab-separated value representing the sample name.

    Returns
    -------
    dict
        A dictionary where:
            - Keys are sample names (str)
            - Values are lists of two directory paths (str):
                [0]: Path to the DGE-filtered matrix directory
                [1]: Path to the Seurat output directory
    """
    parsebio_samples = {}
    with open(sample_list_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('//'):
                continue
            parts = line.split('\t')
            if len(parts) >= 1:
                parsebio_samples[parts[0]] = [os.path.join(analysis, "split_pipe_comb", parts[0], "DGE_filtered/"),
                                               os.path.join(analysis, "split_pipe_comb", parts[0], "seurat/")]
    return parsebio_samples
