#!/bin/bash
#PBS -q small

# Run both methylation and phase_reads without colliding on fasta/bam split.
# Usage: qsub -F "/path/to/masked_genome.fasta /path/to/unmasked_genome.fasta /path/to/reads.fasta /path/to/variants.vcf [/path/to/alignment.bam]" nanopolish_allelic_methylation.sh
#
# NB: nanopolish relies on being able to access the fast5 files that created your fastq.
# /path/to/reads.fasta (or fastq) must be generated using either `poretools` or `nanopolish extract`
# The paths contained in the fasta header must be accessible (without symlinks!) from ~/tmp
# This is done using `rsync -a /path/to/fast5 ~/tmp/path/to/fast5`
# It is advised to call `nanopolish extract` with relative paths, not absolute paths, for this reason

module load nanopolish
set -x
cd $PBS_O_WORKDIR
MASKED_GENOME=$(realpath $1)
UNMASKED_GENOME=$(realpath $2)
FASTA=$(realpath $3)
VCF=$(realpath $4)
BAM=$(realpath $5)
FAST5_DIR=$(realpath $6)

if [ $(echo $FASTA | grep -c -e "a$") -gt 0 ]; then
  FMT="fasta"
elif [ $(echo $FASTA | grep -c -e "q$") -gt 0 ]; then
  FMT="fastq"
else
  echo "ERROR: $FASTA format not recognised"
  exit 1
fi
N=10000
TMP_DIR="$PBS_O_HOME/tmp/$(dirname $FASTA)"
SCRIPTS_DIR=$PBS_O_HOME/nanopore-scripts

mkdir -p $TMP_DIR
cd $TMP_DIR
if [ ! -f $(basename $FASTA).1.$FMT ]; then
  if [ "$#" -gt 4 ]; then
    python $SCRIPTS_DIR/split_bam_and_fasta.py -b $BAM -f $FASTA --prefix $(basename $FASTA) --bam-suffix fastq.sorted.bam -n $N &
  else
    python $SCRIPTS_DIR/split_fasta.py $FASTA $N &
  fi
fi

if [ ! -f ${FASTA}.fa.gz ]; then
  nanopolish index -d $FAST5_DIR $FASTA
fi
wait

$SCRIPTS_DIR/nanopolish_methylation_run.sh $MASKED_GENOME $FASTA
$SCRIPTS_DIR/nanopolish_phase_reads_run.sh $MASKED_GENOME $UNMASKED_GENOME $FASTA $VCF

