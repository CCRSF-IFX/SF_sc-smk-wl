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

snakemake --jobname 's.{jobid}.{rulename}' --latency-wait 600 -j 100 --rerun-incomplete --restart-times 0 --keep-going --stats snakemake.stats --printshellcmds --cluster-config cluster.json --cluster "sbatch --partition={cluster.partition} --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem}  --time={cluster.time} --no-requeue"  --cluster-status /mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/cluster_status.py >& snakemake.log
