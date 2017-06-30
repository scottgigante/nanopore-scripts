#!/bin/bash
#PBS -q small

# Usage: qsub -F "/path/to/genome.fasta /path/to/reads.fasta" nanopolish_methylation_run.sh

set -x
GENOME=$1
FASTA=$2
FMT=$(echo $FASTA | sed 's/.*\.//')
N=10000
TMP_DIR="$HOME/tmp/$(dirname $FASTA)"
SCRIPTS_DIR=$HOME/nanopore-scripts

mkdir -p $TMP_DIR
cd $TMP_DIR
python $SCRIPTS_DIR/split_fasta.py $FASTA $N $FMT
ARRAY_ID=$(qsub -F "$GENOME $FASTA $TMP_DIR" -t 1-$(ls -1 $TMP_DIR/$(basename $FASTA).*.$FMT | wc -l) $SCRIPTS_DIR/nanopolish_methylation.sh)
qsub -W "depend=afteranyarray:$ARRAY_ID" -F "$FASTA $TMP_DIR" $SCRIPTS_DIR/nanopolish_methylation_clean.sh

