#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=30gb
module load gmap-gsnap
module load samtools
./monitor_cpu.sh "gmap.avx" &

set -x
FASTA_DIR=$1
FMT=$2
FASTA=$FASTA_DIR/gmap.${FMT}.${PBS_ARRAYID}.$FMT
gmap -d hg38 -D tmp/gmap/GMAP_database -A -f samse -B 4 -t 1 $FASTA > ${FASTA}.sam || echo "ERROR: gmap failed"
kill %1

samtools sort -T ${FASTA}.samtools -o ${FASTA}.sorted.bam ${FASTA}.sam || echo "ERROR: samtools failed to sort sam file"
echo "COMPLETE"
