#!/bin/bash

# Usage: ./albacore_run.sh [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CONFIG] [-5] [-q] /path/to/fast5
# Defaults: -i raw -o albacore -c r95_450bps_linear.cfg
# Directory structure: 
# Raw reads:  /path/to/fast5/raw
# Basecalled: /path/to/fast5/albacore

# Find the absolute path to this and related scripts (eg ~/nanopore-scripts)
SCRIPTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
# Load anaconda3 module which contains albacore
module load anaconda3
# Check that albacore is installed - if not, quit and throw and error message
if [ $(pip freeze | grep -c "ont-albacore") -eq 0 ]; then
  echo "ERROR:Albacore not installed. Please install from wheel with"
  echo "pip3 install [--user] ont_albacore-x.x.x-cp35-cp35m-manylinux1_x86_64.whl"
  exit 1
fi

# Set default arguments
INPUT="raw"
OUTPUT="albacore"
CONFIG="r95_450bps_linear.cfg"
MODE="fast5"
FAST5=false
FASTQ=false
HELP="albacore_run.sh: Albacore wrapper for Milton
Scott Gigante, 2017-08-01
Usage: ./albacore_run.sh [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CONFIG] [-5] [-q] /path/to/fast5"

# Parse command line arguments. If an argument is missing, it will be left as default
while getopts 'hi:o:c:5q' arg
do
  case ${arg} in
    h) echo "$HELP"; exit 0;;
    i) INPUT=${OPTARG};;
    o) OUTPUT=${OPTARG};;
    c) CONFIG=${OPTARG};;
    5) FAST5=true;;
    q) FASTQ=true;;
    *) echo "$HELP" >&2; exit 1 # illegal option
    esac
done
# Set the output mode
if $FAST5 && $FASTQ; then
  MODE="fast5,fastq"
elif $FASTQ; then
  MODE="fastq"
elif $FAST5; then
  MODE="fast5"
fi

# Shift the beginning argument index so we only have what hasn't been parsed by getopt (eg, /path/to/fast5)
shift $(($OPTIND - 1))
# Check there is exactly one argument left
if [ "$#" -eq 0 ]; then
  echo "ERROR:Directory argument is mandatory."
  echo "$HELP" >&2
  exit 1
elif [ "$#" -gt 1 ]; then
  echo "ERROR:Unrecognised argument: ${@:2}"
  echo "$HELP" >&2
  exit 1
fi
# Set the fast5 parent directory to the absolute path of /path/to/fast5
PARENT=$(realpath $1)

# Check how many subfolders are in /path/to/fast5/raw
NUM_RUNS=$(ls -1 $PARENT/$INPUT | wc -l)
# Make the output directory, if it doesn't already exist
mkdir -p $PARENT/$OUTPUT
# Print some run information
echo "albacore_run.sh: Albacore wrapper for Milton
Scott Gigante, 2017-06-23

Input: $PARENT/$INPUT
Output: $PARENT/$OUTPUT
Config: $CONFIG
Output: $MODE
Jobs: $NUM_RUNS"
# Run the array job
qsub -t 1-$NUM_RUNS -F "$PARENT $INPUT $OUTPUT $CONFIG $MODE" $SCRIPTS_DIR/albacore.sh

