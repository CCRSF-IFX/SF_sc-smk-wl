#!/usr/bin/env python3
import glob
import os
import argparse 
import re
from collections import defaultdict
import psutil
import logging
logging.basicConfig(level = logging.INFO, format='%(asctime)s %(name)s  %(levelname)s: %(message)s', datefmt='%d %b %Y %H:%M:%S')
config = {}

suffix = '_001.fastq.gz'

def print_cmd():
    my_process = psutil.Process( os.getpid() )
    print("# Command line: " + ' '.join(my_process.cmdline()))

def remove_leading_number_trailing_str(sample_name):
    '''
    Remove leading number and trailing lib infor from sample name. E.g "12_SNSCC_2_TCR" = "SNSCC_2"
    ''' 
    sample_name = re.sub("^\d+_", "", sample_name)
    return "_".join(sample_name.split('_')[:-1])

def get_run_folder(path):
    oRunName = re.compile("(\d{6,})_([A-Z0-9]+)_(\d+)_([-A-Z0-9]+)")
    sRunName = oRunName.search(path)
    sRunName_str = path[sRunName.start():sRunName.end()]
    #print(sRunName_str)
    path_sub1, path_sub2 = path.split(sRunName_str+"/")
    #print(path_sub1)
    #print(path_sub2)
    return os.path.join(path_sub1, sRunName_str), path_sub2, sRunName_str

def get_path_list(path_list):
    paths = []
    with open(path_list, "r") as fpath:
        for line in fpath:
            paths.append(line.rstrip("\n"))
    return paths

def print_slurm_config():
    print("""
#!/bin/bash
#SBATCH --partition=norm
#SBATCH --ntasks=30
#SBATCH --mem=20g
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --no-requeue    
""")

def main(options):
    # get the fastq folders 
    paths = options.path.split(",")
    samples = options.sample.split(',')
    with open(options.output, 'w') as fout: 
        for path in paths:
            print("#The folders in " + path + ": ")
            ## list all sub-folders in the unaligned folder
            dir_list = [dI for dI in os.listdir(path) if os.path.isdir(os.path.join(path, dI))] 
            #logging.info(dir_list)
            flow_cell_id = os.path.basename(path)
            run_folder, subfolder, run_name = get_run_folder(path)
            #print([run_folder, subfolder, run_name])
            print(f"rsync -av {run_folder} {options.outdir} --exclude {subfolder} &")
            for dir_name in dir_list:
                sample_name = dir_name
                sample_name4gex_vdj = remove_leading_number_trailing_str(sample_name)
                fq_files = [glob.glob(f'{path}/{sample_name}/{sample_name}_*_{i}{suffix}')[0] for i in ['I1', 'I2', 'R1', 'R2', 'R3'] if len(glob.glob(f'{path}/{sample_name}/{sample_name}_*_{i}{suffix}')) > 0]
                for fqf in fq_files:
                    basename = os.path.basename(fqf)
                    fqf_folder = os.path.join(*[options.outdir, run_name, subfolder, sample_name])
                    fqf_new = os.path.join(fqf_folder, basename)
                    print(f'mkdir -p {fqf_folder}\nzcat {fqf} |head -1000000 |gzip - > {fqf_new} &')
            print("wait")
if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description = "Prepare libraries.csv for `mutlti`")
    parser.add_argument('-p', '--path', help = 'List of paths to the demultiplexed files', required = True)
    parser.add_argument('-s', '--sample', help = 'List of sampels to use', required = True)
    parser.add_argument('--outdir', help = 'Outdir', required = True) 
    parser.add_argument('-o', '--output', help = 'Output file for library info', default = "libraries.csv")

    options = parser.parse_args()
    #print(options.filter_str)
    print_slurm_config() 
    print_cmd()
    os.makedirs(options.outdir, exist_ok = True)
    main(options) 
