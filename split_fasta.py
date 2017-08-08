#!/usr/bin/env python

# Usage: ./split_fasta.py /path/to/fasta.fa 2000
# ./split_fasta.py /path/to/reads.fastq 100 fastq
from __future__ import print_function
from Bio import SeqIO
import sys
import os

fasta = sys.argv[1]
n = int(sys.argv[2])
if len(sys.argv) > 3:
    fmt = sys.argv[3]
else:
    fmt = "fasta"
out_fasta = os.path.join(os.getcwd(), os.path.basename(fasta))

with open(fasta, 'r') as in_handle:
    out_handle = open("{}.1.{}".format(out_fasta, fmt), 'w') # start from 1 so PBS likes us
    i = 0
    for record in SeqIO.parse(in_handle, fmt):
        if i % n == 0 and i > 0:
            out_handle.close()
            out_handle = open("{}.{}.{}".format(out_fasta, (i // n) + 1, fmt), 'w')
        SeqIO.write(record, out_handle, fmt)
        i += 1
out_handle.close()
print("Successfully split {} into {} files.".format(fasta, i // n))
