#!/bin/bash
#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4g
#SBATCH --mail-user=ccrsfifx@nih.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=96:00:00
#SBATCH --no-requeue
#SBATCH --job-name=submit.sh

module load miniconda
eval "$(conda shell.bash activate)"
conda activate /mnt/ccrsf-static/Analysis/xies4/condaenv4plsc
module load singularity 
snakemake --jobname 's.{jobid}.{rulename}' --profile workflow/profile/slurm/  -e cluster-generic >& snakemake.log
