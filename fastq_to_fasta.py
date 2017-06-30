#!/usr/bin/env python

# Usage: ./fastq_to_fasta.py path/to/reads.fastq
from Bio import SeqIO
import sys

in_filename = sys.argv[1]
if not in_filename[-1] == 'q':
    print "Input must be a fastq file"
else:
    out_filename = in_filename[:-1] + 'a'
    with open(in_filename, 'r') as in_handle, open(out_filename, 'w') as out_handle:
        for record in SeqIO.parse(in_handle, 'fastq'):
            SeqIO.write(record, out_handle, 'fasta')

