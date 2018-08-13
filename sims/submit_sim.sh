#!/bin/bash

USAGE="
Usage:
   $0 TYPE
where TYPE should be 'flat', 'biased' or 'barrier'
"

if [ $# -lt 1 ]
then
    echo "$USAGE"
    exit 0
fi

TYPE="$1"

NAME="run"
TAG=$(printf "%06d" $RANDOM); 
OUTDIR="${NAME}_${TAG}"
mkdir -p $OUTDIR
echo "Directory: $OUTDIR"

echo "defineConstant(\"TYPE\", \"${TYPE}\");" >> $OUTDIR/parameters.slim

export OUTDIR
sbatch -o $OUTDIR/slurm_${TAG}.out -e $OUTDIR/slurm_${TAG}.out ./run_sim.sbatch 
# ./run_sim.sbatch 
