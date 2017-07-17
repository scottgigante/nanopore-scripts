#!/bin/bash

# Usage: ./merge_tsv.sh in1.tsv in2.tsv [in3.tsv ...]

cat <(head -n1 $1) <(for i in $@; do tail -n +2 $i; done)
