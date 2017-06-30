#!/bin/bash

# Usage: ./albacore_run.sh [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CONFIG] /path/to/fast5

module load anaconda3
if [ $(pip freeze | grep -c "ont-albacore") -eq 0 ]; then
  echo "ERROR:Albacore not installed. Please install from wheel with"
  echo "pip3 install [--user] ont_albacore-x.x.x-cp35-cp35m-manylinux1_x86_64.whl"
  exit 1
fi

RAW="raw"
READS="albacore"
CONFIG="r94_450bps_linear.cfg"
HELP="albacore_run.sh: Albacore wrapper for Milton
Scott Gigante, 2017-06-23
Usage: ./albacore_run.sh [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CONFIG] /path/to/fast5"

while getopts 'hi:o:c:' arg
do
  case ${arg} in
    h) echo "$HELP"; exit 0;;
    i) RAW=${OPTARG};;
    o) READS=${OPTARG};;
    c) CONFIG=${OPTARG};;
    *) echo "$HELP" >&2; exit 1 # illegal option
    esac
done
shift $(($OPTIND - 1))

if [ "$#" -eq 0 ]; then
  echo "ERROR:Directory argument is mandatory."
  echo "$HELP" >&2
  exit 1
fi
PARENT=$(pwd)/$1

SCRIPTS_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
NUM_RUNS=$(ls -1 $PARENT/$RAW | wc -l)
mkdir -p $PARENT/$READS
echo "albacore_run.sh: Albacore wrapper for Milton
Scott Gigante, 2017-06-23

Input: $PARENT/$RAW
Output: $PARENT/$READS
Config: $CONFIG
Jobs: $NUM_RUNS"
qsub -t 1-$NUM_RUNS -F "$PARENT $RAW $READS $CONFIG" $SCRIPTS_DIR/albacore.sh

