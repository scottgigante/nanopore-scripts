#!/bin/bash
#PBS -q small

# Usage: qsub -F "wehi_home/grpu_mritchie_0/AGRF_Data/2017_05_09/reads/fasta/albacore_no_phage.fasta" gmap_run.sh

set -x
FASTA=$PBS_O_WORKDIR/$1
FMT=$(echo $FASTA | sed 's/.*\.//')
N=2000
FASTA_DIR=$(dirname $FASTA)
TMP_DIR=$PBS_O_HOME/tmp/gmap/$FASTA_DIR
SCRIPTS_DIR=$PBS_O_HOME/nanopore-scripts
mkdir -p $TMP_DIR
cp $FASTA $TMP_DIR/gmap.$FMT
(cmp --silent $FASTA $TMP_DIR/gmap.$FMT && echo "Copied successfully to HPC") || (echo "ERROR: Copying fasta to HPC failed" && exit 1)

cd $TMP_DIR
python $SCRIPTS_DIR/split_fasta.py gmap.$FMT $N $FMT
NUM_JOBS=$(ls -1 gmap.${FMT}.*.$FMT | wc -l)
qsub -t 1-$NUM_JOBS -F "$TMP_DIR $FMT" $SCRIPTS_DIR/gmap.sh

