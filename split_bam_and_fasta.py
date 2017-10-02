#!/usr/bin/env python
# Splits a BAM file and the FASTA/FASTQ used to generate the alignment into matched chunks
# Every read in chunk i of the the BAM file will have its original FASTA/FASTQ record stored in chunk i of the FASTA/Q
# Usage: ./split_bam_and_fasta.py in.bam in.fasta/q chunk_size [suffix]
# Output files: $(pwd)/in.#.fastq, $(pwd)/in.#[suffix].bam

from __future__ import print_function
import sys
import os
import pysam
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bam", required=True)
parser.add_argument("-f", "--fasta", required=True)
parser.add_argument("-n", "--chunk-size", type=int, default=1000)
parser.add_argument("--prefix", default=None)
parser.add_argument("--fasta-suffix", default=None)
parser.add_argument("--bam-suffix", default=None)
args = parser.parse_args()

bam_fn = args.bam
fasta_fn = args.fasta
chunk_size = args.chunk_size
fmt = "fastq" if fasta_fn[-1] == "q" else "fasta"
prefix = args.prefix if args.prefix is not None else os.path.basename(bam_fn).replace(".bam", "")
prefix = os.path.join(os.getcwd(), prefix)
bam_suffix = args.bam_suffix if args.bam_suffix is not None else "bam"
fasta_suffix = args.fasta_suffix if args.fasta_suffix is not None else fmt

chunk_idx = dict()
n_chunks = 1
with open(fasta_fn, 'r') as fasta:
	i = 0
	out_fasta = open("{}.{}.{}".format(prefix, i // chunk_size, fasta_suffix), 'w')
	for read in SeqIO.parse(fasta, fmt):
		chunk_idx[read.id] = i // chunk_size
		SeqIO.write(read, out_fasta, fmt)
		i += 1
		if i % chunk_size == 0:
			out_fasta.close()
			out_fasta = open("{}.{}.{}".format(prefix, i // chunk_size, fasta_suffix), 'w')
			n_chunks += 1
	out_fasta.close()

with pysam.AlignmentFile(bam_fn) as bam:
	out_bams = [pysam.AlignmentFile("{}.{}.{}".format(prefix, j, bam_suffix), 'wb', template=bam) for j in range(n_chunks)]
	for read in bam:
		out_bams[chunk_idx[read.query_name]].write(read)
	for out_bam in out_bams:
		out_bam.close()
