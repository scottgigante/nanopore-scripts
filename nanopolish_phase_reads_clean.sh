#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=32gb

module load samtools
set -x
FASTA=$1
TMP_DIR=$2
if [ $(echo $FASTA | grep -c -e "a$") -gt 0 ]; then
  FMT="fasta"
elif [ $(echo $FASTA | grep -c -e "q$") -gt 0 ]; then
  FMT="fastq"
else
  echo "ERROR: $FASTA format not recognised"
  exit 1
fi

samtools merge ${FASTA}.phased.sorted.bam $TMP_DIR/$(basename $FASTA).*.${FMT}.phased.sorted.bam || echo "ERROR: samtools merge failed"
samtools index ${FASTA}.phased.sorted.bam

echo "COMPLETE"

