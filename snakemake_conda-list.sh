#!/bin/bash

# check for snakemake
command -v snakemake >/dev/null 2>&1 || { echo "snakemake is not in your PATH"; exit 1; }

# check for conda envs
if [[ ! -d .snakemake ]] || [[ -z "$(ls -A .snakemake/conda/)" ]]; then
    echo "No conda envs found!"
    echo "To create the envs, run: 'snakemake --use-conda --create-envs-only -F'"
    exit 1
fi

# list all conda envs
for X in $(snakemake --list-conda-envs -F | tail -n +3)
do
    conda list -p $X 2>/dev/null || echo "#--- conda env: $X ---#"
done

