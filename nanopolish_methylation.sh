#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=16gb
module load bwa
module load nanopolish
module load samtools

set -x
GENOME=$1
FASTA=$2
TMP_DIR=$3
FMT=$(echo $FASTA | sed 's/.*\.//')

INPUT_FASTA=$TMP_DIR/$(basename $FASTA).${PBS_ARRAYID}.$FMT

cd ~/tmp
if [ ! -f ${INPUT_FASTA}.sorted.bam ]; then
  bwa mem -x ont2d $GENOME $INPUT_FASTA | samtools sort -T ${INPUT_FASTA}.tmp -o ${INPUT_FASTA}.sorted.bam
fi
if [ ! -f ${INPUT_FASTA}.sorted.bam.bai ]; then
  samtools index ${INPUT_FASTA}.sorted.bam || echo "ERROR: samtools index failed"
fi
nanopolish call-methylation -r $INPUT_FASTA -b ${INPUT_FASTA}.sorted.bam -g $GENOME > ${INPUT_FASTA}.methylation.tsv || echo "ERROR: nanopolish failed"
echo "COMPLETE"

