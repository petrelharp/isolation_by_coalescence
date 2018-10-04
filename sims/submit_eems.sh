#!/bin/bash

USAGE="
Usage:
   $0 DATABASE
where DATABASE is the file base of the data.
"

if [ $# -lt 1 ]
then
    echo "$USAGE"
    exit 0
fi

DATABASE="$1"
DATADIR="$(dirname $DATABASE)"

NAME="$(basename $DATABASE)_eems"
TAG=$(printf "%06d" $RANDOM); 
OUTDIR="$DATADIR/${NAME}_${TAG}"
mkdir -p $OUTDIR
echo "Directory: $OUTDIR"

OUTBASE="${OUTDIR}/eems_data"
SAMPLEFILE="${DATABASE}.samples.tsv"
DISTFILE="${DATABASE}.divs.tsv"

# eemsify
Rscript eemsify.R $DISTFILE $SAMPLEFILE $OUTBASE

NINDIVS=$(wc -l ${OUTBASE}.coord | cut -f 1 -d ' ')

INIFILE=$OUTDIR/eems.ini
touch $INIFILE

echo "datapath = $OUTBASE" >> $INIFILE
echo "mcmcpath = $OUTDIR" >> $INIFILE
echo "nIndiv = $NINDIVS" >> $INIFILE
echo "nSites = 30000" >> $INIFILE
echo "nDemes = 24" >> $INIFILE
echo "diploid = false" >> $INIFILE
echo "numMCMCIter = 500000" >> $INIFILE
echo "numBurnIter = 100000" >> $INIFILE
echo "numThinIter = 9999" >> $INIFILE

export OUTDIR
sbatch -o $OUTDIR/slurm_${TAG}.out -e $OUTDIR/slurm_${TAG}.out ./run_eems.sbatch 
