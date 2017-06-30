#PBS -q small

module load samtools
set -x
FASTA=$1
FASTA_DIR=$(dirname $FASTA)
TMP_DIR=tmp/gmap/$FASTA_DIR
FMT=$(echo $FASTA | sed 's/.*\.//')
SCRIPTS_DIR=$HOME/nanopore-scripts

for bam in $TMP_DIR/gmap.${FMT}.*.sam; do
  GMAP_FASTA=$(echo $bam | sed 's/.sam$//')
  if [ $(samtools view $bam 2>&1 1> /dev/null | wc -l) -gt 0 ]; then
    qsub -l nodes=1:ppn=2,mem=160gb -F "$GMAP_FASTA" $SCRIPTS_DIR/gmap.sh
  elif [ $(samtools view $bam | wc -l) -eq 0 ]; then
    qsub -l nodes=1:ppn=2,mem=160gb -F "$GMAP_FASTA" $SCRIPTS_DIR/gmap.sh
  elif [ ! -e ${GMAP_FASTA}.sorted.bam ]; then
    samtools sort -T ${GMAP_FASTA}.samtools.tmp -o ${GMAP_FASTA}.sorted.bam ${GMAP_FASTA}.sam || echo "ERROR: samtools sort failed"
  fi
done
