#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=64gb

module load nanopolish
module load samtools
set -x
GENOME=$1
FASTA=$2
TMP_DIR=$3
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
for i in 1..$(ls -1 $TMP_DIR/${basename $FASTA}.*.${FMT} | wc -l}; do
  F="$TMP_DIR/${basename $FASTA}.${i}.${FMT}.methylation.tsv"
  if [ ! -s $F || $(tail -c 1 $F) ]; then
    if [ $RERUN ]; then
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
if [ $ERROR ]; then
  echo "ERROR: No valid methylation data found."
  exit 1
elif [ $RERUN ]; then
  ARRAY_ID=$(qsub -F "$GENOME $FASTA $TMP_DIR" -t $ARRAY $SCRIPTS_DIR/nanopolish_methylation.sh)
  qsub -W "depend=afteranyarray:$ARRAY_ID" -F "$GENOME $FASTA $TMP_DIR" $SCRIPTS_DIR/nanopolish_methylation_clean.sh
  exit 0
fi

# pull out the header
HEADER=$TMP_DIR/nanopolish.methylation.header.tsv
head -n1 $TMP_DIR/$(basename $FASTA).1.${FMT}.methylation.tsv > $HEADER

# concatenate files
$SCRIPTS_DIR/merge_tsv.sh --skip 1 $TMP_DIR/$(basename $FASTA).*.${FMT}.methylation.tsv | tail -n +2 | sort -k1,1 -k2n,2n -k4,4 | cat $TMP_DIR/nanopolish.methylation.header.tsv - >${FASTA}.methylation.tsv &

#merge the split bam files into one single bam
samtools merge ${FASTA}.sorted.bam $TMP_DIR/$(basename $FASTA).*.${FMT}.sorted.bam || echo "ERROR: samtools merge failed"
samtools index ${FASTA}.sorted.bam


# summarise methylation at genomic locations
# python $NANOPOLISH_SCRIPTS/calculate_methylation_frequency.py -c 0 -i ${FASTA}.methylation.tsv > ${FASTA}.methylation.summary.tsv || echo "ERROR: nanopolish calculate_methylation_frequency.py failed. Check for memory failure" &
# python $NANOPOLISH_SCRIPTS/calculate_methylation_frequency.py -c 2 -i ${FASTA}.methylation.filtered.tsv > ${FASTA}.methylation.filtered.summary.tsv || echo "ERROR: nanopolish calculate_methylation_frequency.py failed. Check for memory failure" &
wait
echo "COMPLETE"

