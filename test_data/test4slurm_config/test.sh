# A minimal demo to test the configuration for slurm job submission 
snakemake -s Snakefile --profile ../../profile/slurm/ -e cluster-generic 
rm -rf jobs/ .snakemake/ results/ slurm-* 
