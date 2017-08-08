#!/bin/bash
#PBS -q medium

module load samtools
set -x
FASTA=$1
FASTA_DIR=$(dirname $FASTA)
TMP_DIR=tmp/gmap/$FASTA_DIR
if [ $(echo $FASTA | grep -c -e "a$") -gt 0 ]; then
  FMT="fasta"
elif [ $(echo $FASTA | grep -c -e "q$") -gt 0 ]; then
  FMT="fastq"
else
  echo "ERROR: $FASTA format not recognised"
  exit 1
fi

find $TMP_DIR -name "gmap.${FMT}.*.bam" | {
  read firstbam
  samtools view -h "$firstbam"
  while read bam; do
    samtools view "$bam"
  done
} | samtools view -ubS - | samtools sort -@ $(nproc) -T $TMP_DIR/gmap.samtools.tmp -o $TMP_DIR/gmap.sorted.bam || echo "ERROR: samtools merge"

BAM_DIR=$FASTA_DIR/../alignment
mkdir -p $BAM_DIR
cp $TMP_DIR/gmap.sorted.bam $BAM_DIR/gmap.sorted.bam
cmp --silent $TMP_DIR/gmap.sorted.bam $BAM_DIR/gmap.sorted.bam || (echo "ERROR: Copying bam to wehisan failed" && exit 1)
samtools index $BAM_DIR/gmap.sorted.bam
echo "COMPLETE"
