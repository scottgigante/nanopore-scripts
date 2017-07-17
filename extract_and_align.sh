#!/bin/bash
#PBS -q huge
#
# Usage: ./extract_and_align.sh /path/to/reads RUN_ID /path/to/genome [phage_reads_subdir]
# Only use phage_reads_subdir if a calibration strand has been added

module load bwa
module load samtools
module load nanopolish
module load parallel
module load bamtools
set -x
cd $PBS_O_WORKDIR
READS_DIR=$1
READ_GROUP=$2
GENOME=$3
PARENT_DIR=$(dirname $READS_DIR)
OUT_NAME=$(basename $READS_DIR)
SCRIPTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
NPROC=$(nproc)

if [ "$#" -gt 3 ]; then
  PHAGE_DIR=$PARENT_DIR/$4
  mkdir -p $PHAGE_DIR

  # extract fasta
  FASTA=$PHAGE_DIR/${OUT_NAME}.fasta
  find $READS_DIR -name "*.fast5" | \
    parallel -X nanopolish extract {} > $FASTA

  # align phage
  PHAGE_BAM=$PHAGE_DIR/${OUT_NAME}.sorted.bam
  bwa mem -x ont2d -t $((NPROC-1)) $SCRIPTS_DIR/calibration_strand.fa $FASTA | samtools sort -@ 4 -T $PARENT_DIR/samtools.tmp -o $PHAGE_BAM
  samtools index $PHAGE_BAM
  
  # remove phage
  bamtools split -in $PHAGE_BAM -mapped
  python $SCRIPTS_DIR/remove_phage_reads.py $PHAGE_DIR/${OUT_NAME}.sorted.MAPPED.bam $FASTA $PHAGE_DIR
fi

# extract fastq
mkdir -p $PARENT_DIR/fasta
FASTQ=$PARENT_DIR/fasta/${OUT_NAME}.fastq
find $READS_DIR -name "*.fast5" | \
  parallel -X nanopolish extract -q {} > $FASTQ

# align reads
mkdir -p $PARENT_DIR/alignment
BAM=$PARENT_DIR/alignment/${OUT_NAME}.sorted.bam
bwa mem -x ont2d -t $((NPROC-1)) -R '@RG\tID:$READ_GROUP' $GENOME $FASTQ | samtools sort -@ 4 -T $PARENT_DIR/samtools.tmp -o $BAM
samtools index $BAM
