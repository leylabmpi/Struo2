#!/bin/bash

# user input
if [ "$#" -lt 2 ]; then
    echo "Usage: snakemake_clean.sh config_file [preview|delete]"
    echo "Description: delete all snakemake-generated files"
    echo " preview => preview files to delete"
    echo " delete => delete files"
    echo "NOTE: this only deletes files that snakemake knows about, but that's all that's needed to fully restart the snakemake pipeline"
    exit
fi


FILES=`snakemake --summary --rerun-incomplete --configfile $1 | tail -n+2 | cut -f1`
if [ $1 == "preview" ]; then
    echo "#-- Files to delete --#"
    printf '%s\n' "${FILES[@]}"
elif [ $1 == "delete" ]; then
    echo "#-- Deleting the following files --#"
    printf '%s\n' "${FILES[@]}"
    rm -rf $FILES
else
    echo "$1 not recoginized"
fi 
