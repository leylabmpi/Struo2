#!/bin/bash

# user input
if [ "$#" -lt 2 ]; then
    echo "snakemake_sge.sh config.yaml jobs ..."
    echo " config.yaml : snakemake config"
    echo " jobs : number of parallel qsub jobs"
    echo " ... : additional arguments passed to snakemake"
    exit
fi

# check for snakemake
command -v snakemake >/dev/null 2>&1 || { echo "snakemake is not in your PATH"; exit 1; }

# set args
CONFIG=$1
JOBS=$2

# snakemake call
WORKDIR=`pwd`
snakemake -f \
	  --profile sge \
	  --use-conda \
	  --configfile $CONFIG \
	  --jobs $JOBS \
	  --local-cores $JOBS \
	  --printshellcmds \
	  --resources temp=$JOBS \
	  --directory $WORKDIR \
	  "${@:3}"


