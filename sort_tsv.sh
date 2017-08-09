#!/bin/bash

# Sort a tsv file without destroying headers
# Assumes headers start with a hash or length of header is provided by --skip
# Usage: ./sort_tsv.sh infile.tsv [--skip n] <args for sort>
# Help: `man sort`
INFILE=$1
if [ $2 == "--skip" ]; then
  i=$3
  shift 3
else
  # check header length
  i=0
  while [ true ]; do
    if [ $(head -n $i $INFILE | grep -c -e "^#") -lt $i ]; then
      i=$(($i-1))
      break
    fi
    i=$(($i+1))
  done
  shift 1
fi
ARGS=$@

tail -n +$(($i+1)) $INFILE | sort $ARGS | cat <(head -n $i $INFILE) -
