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
INPUT=$2
OUTPUT=$3
CONFIG=$4
MODE=$5
# Choose the nth subdirectory of /path/to/fast5/raw - n is specified in $PBS_ARRAYID
SUBDIR=$(ls -1 $PARENT/$INPUT | sed "${PBS_ARRAYID}q;d")
# Set paths on wehisan
INPUT_DIR="$PARENT/$INPUT/$SUBDIR"
OUTPUT_DIR="$PARENT/$OUTPUT/$SUBDIR"
# Set paths on HPCScratch
TMP_INPUT_DIR="$PBS_O_HOME/tmp/$INPUT_DIR"
TMP_OUTPUT_DIR="$PBS_O_HOME/tmp/$OUTPUT_DIR"
# Make directories if they do not exist
mkdir -p $TMP_INPUT_DIR
mkdir -p $TMP_OUTPUT_DIR
mkdir -p $OUTPUT_DIR

# Count the number of raw reads
N_INPUT=$(find $INPUT_DIR -name "*.fast5" | wc -l)
echo "[$(date +'%y-%m-%d %H:%M:%S')] Copying files from $INPUT_DIR to $TMP_INPUT_DIR..."
# Copy raw reads to tmp
find $INPUT_DIR -name "*.fast5" | parallel -X cp {} $TMP_INPUT_DIR
# Count the number of moved reads
N_MOVED=$(find $TMP_INPUT_DIR -name "*.fast5" | wc -l)
# Check that all reads were moved
if [ $N_MOVED -lt $N_INPUT ]; then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR: Copy to HPC failed."
  echo "$INPUT_DIR: $N_INPUT"
  echo "$TMP_INPUT_DIR: $N_MOVED"
fi

echo "[$(date +'%y-%m-%d %H:%M:%S')] Basecalling fast5 into $TMP_OUTPUT_DIR..."
# This is it! Basecall the reads
read_fast5_basecaller.py -o $MODE -i $TMP_INPUT_DIR -c $CONFIG -t $(($(nproc)-1)) -s $TMP_OUTPUT_DIR

# Count the number of basecalled reads - even reads without basecalls should have a file associated with them
N_CALLED=$(find $TMP_OUTPUT_DIR/workspace -mindepth 2 -name "*.fast5" | wc -l)
if [[ $MODE == *"fastq"* ]]; then
  N_FASTQ_CALLED=$(grep -c -e "^@" $TMP_OUTPUT_DIR/workspace/*.fastq)
else
  N_FASTQ_CALLED=0
fi
# Check that all reads were basecalled
if ([[ $MODE == *"fast5"* ]] && [ $N_CALLED -lt $N_MOVED ]) || ( [[ $MODE == *"fastq"* ]] && [ $N_FASTQ_CALLED -lt $N_MOVED ]); then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR:Basecalling failed."
  echo "$TMP_INPUT_DIR: $N_MOVED"
  echo "$TMP_OUTPUT_DIR/*.fast5: $N_CALLED"
  echo "$TMP_OUTPUT_DIR/*.fastq: $N_FASTQ_CALLED"
fi
# Remove temporary raw files
rm -rf $TMP_INPUT_DIR

echo "[$(date +'%y-%m-%d %H:%M:%S')] Copying files from $TMP_OUTPUT_DIR to $OUTPUT_DIR..."
# Copy basecalled reads to wehisan
find $TMP_OUTPUT_DIR/workspace -mindepth 2 -name "*.fast5" | parallel -X cp {} $OUTPUT_DIR
# Count the number of fast5 files moved
N_MOVED=$(find $OUTPUT_DIR -name "*.fast5" | wc -l)
# Copy fastq to wehisan - in future, we may want to concatenate these in cleanup
if [[ $MODE == *"fastq"* ]]; then
  find $TMP_OUTPUT_DIR/workspace -name "*.fastq" -exec cat {} \; > ${OUTPUT_DIR}.fastq
  N_FASTQ_MOVED=$(grep -c -e "^@" ${OUTPUT_DIR}.fastq)
else
  N_FASTQ_MOVED=0
fi
# Check if all basecalled files were moved
if [ $N_MOVED -lt $N_CALLED ] || [ $N_FASTQ_MOVED -lt $N_FASTQ_CALLED ]; then
  echo "[$(date +'%y-%m-%d %H:%M:%S')] ERROR:Copy to wehisan failed."
  echo "$TMP_OUTPUT_DIR/*.fast5: $N_CALLED"
  echo "$OUTPUT_DIR/*.fast5: $N_MOVED"
  echo "$TMP_OUTPUT_DIR/*.fastq: $N_FASTQ_CALLED"
  echo "${OUTPUT_DIR}.fastq: $N_FASTQ_MOVED"
else
  # All good, we can remove the temporary basecalled files. If copying failed, user may want to retrieve them.
  rm -rf $TMP_OUTPUT_DIR
fi
echo "[$(date +'%y-%m-%d %H:%M:%S')] COMPLETE"

