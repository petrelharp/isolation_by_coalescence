#!/bin/bash

USAGE="
Usage:
   $0 OUTDIR MUTRATE WIDTH HEIGHT
"

if [ $# -lt 4 ]
then
    echo "$USAGE"
    exit 0
fi

OUTDIR="$1"
MUTRATE="$2"
WIDTH="$3"
HEIGHT="$4"

export OUTDIR
export MUTRATE
export WIDTH
export HEIGHT
sbatch -o $OUTDIR/slurm_divergences_${TAG}.out -e $OUTDIR/slurm_divergences_${TAG}.out ./run_divergences.sbatch 
