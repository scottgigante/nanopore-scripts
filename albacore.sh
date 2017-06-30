#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=16gb

# Usage: ./albacore_run.sh [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CONFIG] /path/to/fast5

module load anaconda3
module load parallel

set -x
PARENT=$1
RAW=$2
READS=$3
CONFIG=$4
SUBDIR=$(ls -1 $PARENT/$RAW | sed "${PBS_ARRAYID}q;d")
RAW_DIR="$PARENT/$RAW/$SUBDIR"
READS_DIR="$PARENT/$READS/$SUBDIR"
TMP_RAW_DIR="tmp/$RAW_DIR"
TMP_READS_DIR="tmp/$READS_DIR"
mkdir -p $TMP_RAW_DIR
mkdir -p $TMP_READS_DIR
mkdir -p $READS_DIR

N_RAW=$(find $RAW_DIR -name "*.fast5" | wc -l)
echo "[$(date +'%y-%m-%d %H:%M:%S')] Copying files from $RAW_DIR to $TMP_RAW_DIR..."
find $RAW_DIR -name "*.fast5" | parallel -X cp {} $TMP_RAW_DIR
N_MOVED=$(find $TMP_RAW_DIR -name "*.fast5" | wc -l)
if [ $N_MOVED -lt $N_RAW ]; then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR: Copy to HPC failed."
  echo "$RAW_DIR: $N_RAW"
  echo "$TMP_RAW_DIR: $N_MOVED"
fi

echo "[$(date +'%y-%m-%d %H:%M:%S')] Basecalling fast5 into $TMP_READS_DIR..."
read_fast5_basecaller.py -o fast5,fastq -i $TMP_RAW_DIR -c $CONFIG -t 1 -s $TMP_READS_DIR

N_CALLED=$(find $TMP_READS_DIR/workspace -mindepth 2 -name "*.fast5" | wc -l)
if [ $N_CALLED -lt $N_MOVED ]; then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR:Basecalling failed."
  echo "$TMP_RAW_DIR: $N_MOVED"
  echo "$TMP_READS_DIR: $N_CALLED"
fi
rm -rf $TMP_RAW_DIR

echo "[$(date +'%y-%m-%d %H:%M:%S')] Copying files from $TMP_READS_DIR to $READS_DIR..."
find $TMP_READS_DIR/workspace -mindepth 2 -name "*.fast5" | parallel -X cp {} $READS_DIR
find $TMP_READS_DIR/workspace -name "*.fastq" -exec cat {} \; > ${READS_DIR}.fastq
N_MOVED=$(find $READS_DIR -name "*.fast5" | wc -l)
if [ $N_MOVED -lt $N_CALLED ]; then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR:Copy to wehisan failed."
  echo "$TMP_READS_DIR: $N_CALLED"
  echo "$READS_DIR: $N_MOVED"
else
  rm -rf $TMP_READS_DIR
fi
echo "[$(date +'%y-%m-%d %H:%M:%S')] COMPLETE"

