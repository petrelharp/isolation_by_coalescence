#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=48G
#SBATCH --ntasks-per-core=1

cd $SLURM_SUBMIT_DIR

: ${OUTDIR?Must define OUTDIR}
if [ -z "$OUTDIR" ]
then
    echo "Must define OUTDIR (is empty)."
    exit 1
fi

INIFILE="$OUTDIR/eems.ini"

echo "Running:"
echo "   eems/runeems_snps/src/runeems_snps --params $INIFILE &> $OUTDIR/eems_run.log"

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    eems/runeems_snps/src/runeems_snps --params $INIFILE &> $OUTDIR/eems_run.log

echo "Done!"

