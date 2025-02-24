#!/bin/bash
#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4g
#SBATCH --mail-user=$USER@nih.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=96:00:00
#SBATCH --no-requeue
#SBATCH --job-name=scMastro


snakemake --jobname 's.{jobid}.{rulename}' --profile workflow/profile/slurm/  -e cluster-generic > snakemake.log 2>&1 
