#!/bin/bash
#PBS -q small

# Usage: qsub -F "/path/to/genome.fasta /path/to/reads.fasta" nanopolish_methylation_run.sh

set -x
GENOME=$PBS_O_WORKDIR/$1
FASTA=$PBS_O_WORKDIR/$2
FMT=$(echo $FASTA | sed 's/.*\.//')
N=10000
TMP_DIR="$PBS_O_HOME/tmp/$(dirname $FASTA)"
SCRIPTS_DIR=$PBS_O_HOME/nanopore-scripts

mkdir -p $TMP_DIR
cd $TMP_DIR
if [ ! -f $(basename $FASTA).1.$FMT ]; then
  python $SCRIPTS_DIR/split_fasta.py $FASTA $N $FMT
fi
ARRAY_ID=$(qsub -F "$GENOME $FASTA $TMP_DIR" -t 1-$(ls -1 $TMP_DIR/$(basename $FASTA).*.$FMT | wc -l) $SCRIPTS_DIR/nanopolish_methylation.sh)
qsub -W "depend=afteranyarray:$ARRAY_ID" -F "$FASTA $TMP_DIR" $SCRIPTS_DIR/nanopolish_methylation_clean.sh

