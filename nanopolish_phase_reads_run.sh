#!/bin/bash
#PBS -q small

# Usage: qsub -F "/path/to/genome.fasta /path/to/reads.fasta /path/to/variants.vcf" nanopolish_phase_reads_run.sh

set -x
cd $PBS_O_WORKDIR
GENOME=$(realpath $1)
FASTA=$(realpath $2)
VCF=$(realpath $3)
FMT=$(echo $FASTA | sed 's/.*\.//')
N=10000
TMP_DIR="$PBS_O_HOME/tmp/$(dirname $FASTA)"
SCRIPTS_DIR=$PBS_O_HOME/nanopore-scripts

mkdir -p $TMP_DIR
cd $TMP_DIR
if [ ! -f $(basename $FASTA).1.$FMT ]; then
  python $SCRIPTS_DIR/split_fasta.py $FASTA $N $FMT
fi
ARRAY_ID=$(qsub -F "$GENOME $FASTA $VCF $TMP_DIR" -t 1-$(ls -1 $TMP_DIR/$(basename $FASTA).*.$FMT | wc -l) $SCRIPTS_DIR/nanopolish_phase_reads.sh)
qsub -W "depend=afteranyarray:$ARRAY_ID" -F "$FASTA $TMP_DIR" $SCRIPTS_DIR/nanopolish_phase_reads_clean.sh

