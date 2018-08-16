#!/bin/bash

USAGE="
Usage:
   $0 OUTDIR MUTRATE
"

if [ $# -lt 2 ]
then
    echo "$USAGE"
    exit 0
fi

OUTDIR="$1"
MUTRATE="$2"

export OUTDIR
sbatch -o $OUTDIR/slurm_divergences_${TAG}.out -e $OUTDIR/slurm_divergences_${TAG}.out ./run_divergences.sbatch 
