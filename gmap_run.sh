#!/bin/bash
#PBS -q small

# Usage: qsub -F "wehi_home/grpu_mritchie_0/AGRF_Data/2017_05_09/reads/fasta/albacore_no_phage.fasta" gmap_run.sh

set -x
cd $PBS_O_WORKDIR
FASTA=$(realpath $1)
if [ $(echo $FASTA | grep -c -e "a$") -gt 0 ]; then
  FMT="fasta"
elif [ $(echo $FASTA | grep -c -e "q$") -gt 0 ]; then
  FMT="fastq"
else
  echo "ERROR: $FASTA format not recognised"
  exit 1
fi
N=2000
FASTA_DIR=$(dirname $FASTA)
TMP_DIR=$PBS_O_HOME/tmp/gmap/$FASTA_DIR
SCRIPTS_DIR=$PBS_O_HOME/nanopore-scripts
mkdir -p $TMP_DIR
cp $FASTA $TMP_DIR/gmap.$FMT
(cmp --silent $FASTA $TMP_DIR/gmap.$FMT && echo "Copied successfully to HPC") || (echo "ERROR: Copying fasta to HPC failed" && exit 1)

cd $TMP_DIR
python $SCRIPTS_DIR/split_fasta.py gmap.$FMT $N
NUM_JOBS=$(ls -1 gmap.${FMT}.*.$FMT | wc -l)
qsub -t 1-$NUM_JOBS -F "$TMP_DIR $FMT" $SCRIPTS_DIR/gmap.sh

