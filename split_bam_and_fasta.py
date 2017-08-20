#!/usr/bin/env python
# Splits a BAM file and the FASTA/FASTQ used to generate the alignment into matched chunks
# Every read in chunk i of the the BAM file will have its original FASTA/FASTQ record stored in chunk i of the FASTA/Q
# Usage: ./split_bam_and_fasta.py in.bam in.fasta/q chunk_size 

from __future__ import print_function
import sys
import pysam
from Bio import SeqIO

bam_fn = sys.argv[1]
fasta_fn = sys.argv[2]
chunk_size = int(sys.argv[3])
fmt = "fastq" if fasta_fn[-1] == "q" else "fasta"
out_fn = bam_fn.replace(".bam", "")

chunk_idx = dict()
with open(fasta_fn, 'r') as fasta:
	i = 0
	out_fasta = open("{}.{}.{}".format(out_fn, i // chunk_size, fmt), 'w')
	for read in SeqIO.parse(fasta, fmt):
		chunk_idx[read.id] = i // chunk_size
		SeqIO.write(read, out_fasta, fmt)
		i += 1
		if i % chunk_size == 0:
			out_fasta.close()
			out_fasta = open("{}.{}.{}".format(out_fn, i // chunk_size, fmt), 'w')
	out_fasta.close()

with pysam.AlignmentFile(bam_fn) as bam:
	out_bams = [pysam.AlignmentFile("{}.{}.bam".format(out_fn, j), 'wb', template=bam) for j in range(i // chunk_size + 1)]
	for read in bam:
		out_bams[chunk_idx[read.query_name]].write(read)
	for out_bam in out_bams:
		out_bam.close()
