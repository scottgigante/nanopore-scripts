#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=16gb
module load bwa
module load nanopolish
module load samtools

set -x
GENOME=$1
FASTA=$2
VCF=$3
TMP_DIR=$4
if [ $(echo $FASTA | grep -c -e "a$") -gt 0 ]; then
  FMT="fasta"
elif [ $(echo $FASTA | grep -c -e "q$") -gt 0 ]; then
  FMT="fastq"
else
  echo "ERROR: $FASTA format not recognised"
  exit 1
fi

INPUT_FASTA=$TMP_DIR/$(basename $FASTA).${PBS_ARRAYID}.$FMT

cd ~/tmp
if [ ! -f ${INPUT_FASTA}.sorted.bam ]; then
  bwa mem -x ont2d $GENOME $INPUT_FASTA | samtools sort -T ${INPUT_FASTA}.tmp -o ${INPUT_FASTA}.sorted.bam
fi
if [ ! ${INPUT_FASTA}.sorted.bam.bai ]; then
  samtools index ${INPUT_FASTA}.sorted.bam || echo "ERROR: samtools index failed"
fi
nanopolish phase-reads -r $INPUT_FASTA -b ${INPUT_FASTA}.sorted.bam -g $GENOME $VCF | samtools sort -T ${INPUT_FASTA}.tmp -o ${INPUT_FASTA}.phased.sorted.bam || echo "ERROR: nanopolish failed"
samtools view ${INPUT_FASTA}.phased.sorted.bam > /dev/null || echo "ERROR: samtools view failed, bam corrupt"
echo "COMPLETE"

