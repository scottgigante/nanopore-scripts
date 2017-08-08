#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=64gb

module load nanopolish
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
SCRIPTS_DIR=$PBS_O_HOME/nanopore-scripts

# pull out the header
HEADER=$TMP_DIR/nanopolish.methylation.header.tsv
head -n1 $TMP_DIR/$(basename $FASTA).1.${FMT}.methylation.tsv > $HEADER
BINARY_HEADER=$TMP_DIR/nanopolish.methylation.header.with_binary.tsv
sed 's/$/\tmethylated/' $HEADER > $BINARY_HEADER

# concatenate files
$SCRIPTS_DIR/merge_tsv.sh $TMP_DIR/$(basename $FASTA).*.${FMT}.methylation.tsv > ${FASTA}.methylation.tsv

tail -n +2 ${FASTA}.methylation.tsv | awk '{ if (sqrt($5^2) > 2) { print; } }' | cat $HEADER - > ${FASTA}.methylation.filtered.tsv

# add + and - for bismark-style viewing in seqmonk
tail -n +2 ${FASTA}.methylation.tsv | awk '{ if ($5 > 0) { print $0 "\t+"; } else { print $0 "\t-"; } }' | cat $BINARY_HEADER - > ${FASTA}.methylation.with_binary.tsv &
tail -n +2 ${FASTA}.methylation.filtered.tsv | awk '{ if ($5 > 0) { print $0 "\t+"; } else { print $0 "\t-"; } }' | cat $BINARY_HEADER - > ${FASTA}.methylation.filtered.with_binary.tsv &

# summarise methylation at genomic locations
python $NANOPOLISH_SCRIPTS/calculate_methylation_frequency.py -c 0 -i ${FASTA}.methylation.tsv > ${FASTA}.methylation.summary.tsv || echo "ERROR: nanopolish calculate_methylation_frequency.py failed. Check for memory failure" &
python $NANOPOLISH_SCRIPTS/calculate_methylation_frequency.py -c 2 -i ${FASTA}.methylation.filtered.tsv > ${FASTA}.methylation.filtered.summary.tsv || echo "ERROR: nanopolish calculate_methylation_frequency.py failed. Check for memory failure" &
wait
echo "COMPLETE"

