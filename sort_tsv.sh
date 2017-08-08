#!/bin/bash

# Sort a tsv file without destroying headers
# Assumes headers start with a hash
# Usage: ./sort_tsv.sh infile.tsv <args for sort>
# Help: `man sort`
INFILE=$1
shift 1
ARGS=$@

# check header length
i=0
while [ true ]; do
  if [ $(head -n $i $INFILE | grep -c -e "^#") -lt $i ]; then
    i=$(($i-1))
    break
  fi
  i=$(($i+1))
done

set -x
tail -n +$(($i+1)) $INFILE | sort $ARGS | cat <(head -n $i $INFILE) -
