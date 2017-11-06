#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=32gb

module load samtools
set -x
MASKED_GENOME=$1
UNMASKED_GENOME=$2
FASTA=$3
VCF=$4
TMP_DIR=$5
if [ $(echo $FASTA | grep -c -e "a$") -gt 0 ]; then
  FMT="fasta"
elif [ $(echo $FASTA | grep -c -e "q$") -gt 0 ]; then
  FMT="fastq"
else
  echo "ERROR: $FASTA format not recognised"
  exit 1
fi
SCRIPTS_DIR=$PBS_O_HOME/nanopore-scripts

RERUN=false
ERROR=true
ARRAY=""
for i in $(eval echo "{1..$(ls -1 $TMP_DIR/${basename $FASTA}.*.${FMT} | wc -l)}"); do
  F="$TMP_DIR/${basename $FASTA}.${i}.${FMT}.phased.sorted.bam"
  if [ $(find $TMP_DIR -name "${basename $FASTA}.${i}.${FMT}.phased.sorted.bam" -size -4096c) $(samtools view $F || echo "ERROR") ]; then
    if $RERUN; then
      ARRAY="${ARRAY},"
    else
      RERUN=true
    fi
    ARRAY="${ARRAY}${i}"
    rm $F
  else
    ERROR=false
  fi
done
if $ERROR; then
  echo "ERROR: No valid phasing data found."
  exit 1
elif $RERUN; then
  ARRAY_ID=$(qsub -F "$MASKED_GENOME $UNMASKED_GENOME $FASTA $VCF $TMP_DIR" -t $ARRAY $SCRIPTS_DIR/nanopolish_phase_reads.sh)
  qsub -W "depend=afteranyarray:$ARRAY_ID" -F "$MASKED_GENOME $UNMASKED_GENOME $FASTA $VCF $TMP_DIR" $SCRIPTS_DIR/nanopolish_phase_reads_clean.sh
  exit 0
fi

samtools merge ${FASTA}.sorted.bam $TMP_DIR/$(basename $FASTA).*.${FMT}.sorted.bam || echo "ERROR: samtools merge failed"
samtools index ${FASTA}.sorted.bam

samtools merge ${FASTA}.phased.sorted.bam $TMP_DIR/$(basename $FASTA).*.${FMT}.phased.sorted.bam || echo "ERROR: samtools merge failed"
samtools index ${FASTA}.phased.sorted.bam

echo "COMPLETE"

