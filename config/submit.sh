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

conda activate /mnt/ccrsf-ifx/Software/tools/conda_env4scwf
module load singularity/4.1.5
snakemake --jobname 's.{jobid}.{rulename}' --profile workflow/profile/slurm/  -e cluster-generic > snakemake.log 2>&1
