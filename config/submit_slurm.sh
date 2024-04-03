#!/bin/bash
#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4g
#SBATCH --mail-user=tsai-wei.shen@nih.gov,biraj.shrestha@nih.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=96:00:00
#SBATCH --no-requeue
#SBATCH --job-name=submit.sh

#2019/12/27 shent2: --no-requeue is added
#2020/02/19 chenv3: added in --restart-times 4
#export PATH=/mnt/nasapps/development/python/3.7.1/bin:$PATH
#export PATH=/mnt/nasapps/development/python/3.11.0/bin:$PATH
#snakemake --jobname 's.{jobid}.{rulename}' -k --stats snakemake.stats --rerun-incomplete --restart-times 4 -j 300 --cluster 'qsub {params.batch}'  >&  snakemake.log
snakemake --jobname 's.{jobid}.{rulename}' --latency-wait 600 -j 100 --rerun-incomplete --restart-times 4 --keep-going --stats snakemake.stats --printshellcmds --cluster-config cluster.json --cluster "sbatch --partition={cluster.partition} --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem}  --time={cluster.time} --no-requeue"  --cluster-status /mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/cluster_status.py >& snakemake.log
