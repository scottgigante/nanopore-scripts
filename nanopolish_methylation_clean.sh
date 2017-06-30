#!/bin/bash
#PBS -q small

set -x
FASTA=$1
TMP_DIR=$2
FMT=$(echo $FASTA | sed 's/.*\.//')

# pull out the header
HEADER=$TMP_DIR/nanopolish.methylation.header.tsv
head -n1 $TMP_DIR/$(basename $FASTA).1.$FMT > $HEADER
BINARY_HEADER=$TMP_DIR/nanopolish.methylation.header.with_binary.tsv
sed 's/$/\tmethylated/' $HEADER > $BINARY_HEADER

# concatenate files
cat $HEADER > ${FASTA}.methylation.tsv
for i in $TMP_DIR/$(basename $FASTA).*.${FMT}.methylation.tsv; do
  tail -n +2 $i >> ${FASTA}.methylation.tsv
done

tail -n +2 ${FASTA}.methylation.tsv | awk '{ if (sqrt($5^2) > 2) { print; } }' | cat $HEADER - > ${FASTA}.methylation.filtered.tsv

# add + and - for bismark-style viewing in seqmonk
tail -n +2 ${FASTA}.methylation.tsv | awk '{ if ($5 > 0) { print $0 "\t+"; } else { print $0 "\t-"; } }' | cat $BINARY_HEADER - > ${FASTA}.methylation.with_binary.tsv &
tail -n +2 ${FASTA}.methylation.filtered.tsv | awk '{ if ($5 > 0) { print $0 "\t+"; } else { print $0 "\t-"; } }' | cat $BINARY_HEADER - > ${FASTA}.methylation.filtered.with_binary.tsv &

# summarise methylation at genomic locations
python $NANOPOLISH_SCRIPTS/calculate_methylation_frequency.py -c 0 -i ${FASTA}.methylation.tsv > ${FASTA}.methylation.summary.tsv &
python $NANOPOLISH_SCRIPTS/calculate_methylation_frequency.py -c 2 -i ${FASTA}.methylation.filtered.tsv > ${FASTA}.methylation.filtered.summary.tsv &
wait
echo "COMPLETE"

