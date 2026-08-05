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

#conda activate /mnt/ccrsf-ifx/Software/tools/conda_env4scwf
#module load singularity/4.1.5
SF_SC_WORKFLOW_VERSION="{sf_sc_workflow_version}"

module purge
module load /mnt/ccrsf-ifx/Software/modulefiles/sfifx.tcl
# Check if the module file exists before loading
if module avail sf-sc-workflow/"${SF_SC_WORKFLOW_VERSION}" 2>&1 | grep -q "${SF_SC_WORKFLOW_VERSION}"; then
    module load sf-sc-workflow/"${SF_SC_WORKFLOW_VERSION}"
else
    echo "ERROR: sf-sc-workflow/${SF_SC_WORKFLOW_VERSION} module not found."
    exit 1
fi

module list
SINGULARITY_WORK="/scratch/ccrsf_scratch/scratch/ccrsfifx/SINGULARITY_WORK"
export SINGULARITY_CACHEDIR="${SINGULARITY_WORK}/cache/${USER}"
mkdir -p "$SINGULARITY_CACHEDIR"
chmod 700 "$SINGULARITY_CACHEDIR" 
chmod g-s "$SINGULARITY_CACHEDIR" 
snakemake --jobname 's.{jobid}.{rulename}' --profile workflow/profile/slurm/  -e cluster-generic > snakemake.log 2>&1
