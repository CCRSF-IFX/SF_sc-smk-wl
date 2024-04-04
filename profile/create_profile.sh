# create config directory that snakemake searches for profiles (or use something else)
profile_dir="./slurm/"
mkdir -p "$profile_dir"
# use cookiecutter to create the profile in the config directory
template="gh:Snakemake-Profiles/slurm"
cookiecutter --output-dir "$profile_dir" "$template"
