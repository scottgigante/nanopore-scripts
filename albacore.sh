#!/bin/bash
#PBS -l nodes=1:ppn=2,mem=16gb

# Usage: ./albacore_run.sh [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CONFIG] [-5] [-q] /path/to/fast5
# Subprogram usage: qsub -F "/absolute/path/to/fast5 raw albacore r94_450bps_linear.cfg fast5,fastq" albacore.sh 

# Load anaconda3 module, where albacore is installed
module load anaconda3
# Load parallel for copying
module load parallel

# Start printing the commands run to stderr for debugging
set -x
# Read command line args
PARENT=$1
RAW=$2
READS=$3
CONFIG=$4
MODE=$4
# Choose the nth subdirectory of /path/to/fast5/raw - n is specified in $PBS_ARRAYID
SUBDIR=$(ls -1 $PARENT/$RAW | sed "${PBS_ARRAYID}q;d")
# Set paths on wehisan
RAW_DIR="$PARENT/$RAW/$SUBDIR"
READS_DIR="$PARENT/$READS/$SUBDIR"
# Set paths on HPCScratch
TMP_RAW_DIR="$PBS_O_HOME/tmp/$RAW_DIR"
TMP_READS_DIR="$PBS_O_HOME/tmp/$READS_DIR"
# Make directories if they do not exist
mkdir -p $TMP_RAW_DIR
mkdir -p $TMP_READS_DIR
mkdir -p $READS_DIR

# Count the number of raw reads
N_RAW=$(find $RAW_DIR -name "*.fast5" | wc -l)
echo "[$(date +'%y-%m-%d %H:%M:%S')] Copying files from $RAW_DIR to $TMP_RAW_DIR..."
# Copy raw reads to tmp
find $RAW_DIR -name "*.fast5" | parallel -X cp {} $TMP_RAW_DIR
# Count the number of moved reads
N_MOVED=$(find $TMP_RAW_DIR -name "*.fast5" | wc -l)
# Check that all reads were moved
if [ $N_MOVED -lt $N_RAW ]; then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR: Copy to HPC failed."
  echo "$RAW_DIR: $N_RAW"
  echo "$TMP_RAW_DIR: $N_MOVED"
fi

echo "[$(date +'%y-%m-%d %H:%M:%S')] Basecalling fast5 into $TMP_READS_DIR..."
# This is it! Basecall the reads
read_fast5_basecaller.py -o $MODE -i $TMP_RAW_DIR -c $CONFIG -t $(($(nproc)-1)) -s $TMP_READS_DIR

# Count the number of basecalled reads - even reads without basecalls should have a file associated with them
N_CALLED=$(find $TMP_READS_DIR/workspace -mindepth 2 -name "*.fast5" | wc -l)
# Check that all reads were basecalled
if [ $N_CALLED -lt $N_MOVED ]; then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR:Basecalling failed."
  echo "$TMP_RAW_DIR: $N_MOVED"
  echo "$TMP_READS_DIR: $N_CALLED"
fi
# Remove temporary raw files
rm -rf $TMP_RAW_DIR


echo "[$(date +'%y-%m-%d %H:%M:%S')] Copying files from $TMP_READS_DIR to $READS_DIR..."
# Copy basecalled reads to wehisan
find $TMP_READS_DIR/workspace -mindepth 2 -name "*.fast5" | parallel -X cp {} $READS_DIR
# Copy fastq to wehisan - in future, we may want to concatenate these in cleanup
find $TMP_READS_DIR/workspace -name "*.fastq" -exec cat {} \; > ${READS_DIR}.fastq
# Count the number of fast5 files moved
N_MOVED=$(find $READS_DIR -name "*.fast5" | wc -l)
# Check if all basecalled files were moved
if [ $N_MOVED -lt $N_CALLED ]; then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR:Copy to wehisan failed."
  echo "$TMP_READS_DIR: $N_CALLED"
  echo "$READS_DIR: $N_MOVED"
else
  # All good, we can remove the temporary basecalled files. If copying failed, user may want to retrieve them.
  rm -rf $TMP_READS_DIR
fi
echo "[$(date +'%y-%m-%d %H:%M:%S')] COMPLETE"

