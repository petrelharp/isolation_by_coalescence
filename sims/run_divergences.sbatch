#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=48G
#SBATCH --ntasks-per-core=1

module use /projects/apps/shared/modulefiles
module load python3 tskit SLiM

cd $SLURM_SUBMIT_DIR

: ${OUTDIR?Must define OUTDIR}
if [ -z "$OUTDIR" ]
then
    echo "Must define OUTDIR (is empty)."
    exit 1
fi

: ${MUTRATE?Must define MUTRATE}
if [ -z "$MUTRATE" ]
then
    echo "Must define MUTRATE (is empty)."
    exit 1
fi

: ${WIDTH?Must define WIDTH}
if [ -z "$WIDTH" ]
then
    echo "Must define WIDTH (is empty)."
    exit 1
fi

: ${HEIGHT?Must define HEIGHT}
if [ -z "$HEIGHT" ]
then
    echo "Must define HEIGHT (is empty)."
    exit 1
fi

echo "Running:"
echo "  ./get_divergences.py -o $OUTDIR -u $MUTRATE -n $WIDTH -m $HEIGHT"

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
     ./get_divergences.py -o $OUTDIR -u $MUTRATE -n $WIDTH -m $HEIGHT &> $OUTDIR/run_divergences.log

echo "Done!"

