#!/bin/bash
# Merge tsv files
# Assumes headers start with a hash or header is length n
# Retains only header from first file
# Usage: ./merge_tsv.sh [--skip n] in1.tsv in2.tsv [in3.tsv ...]

# check header length
if [ $1 == "--skip" ]; then
  i=$2
  shift 2
else
  i=0
  while [ true ]; do
    if [ $(head -n $i $1 | grep -c -e "^#") -lt $i ]; then
      i=$(($i-1))
      break
    fi
    i=$(($i+1))
  done
fi

cat <(head -n $i $1) <(for f in $@; do tail -n +$(($i+1)) $f; done)
