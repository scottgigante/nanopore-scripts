#!/usr/bin/env python
# Splits a BAM file and the FASTA/FASTQ used to generate the alignment into matched chunks
# Every read in chunk i of the the BAM file will have its original FASTA/FASTQ record stored in chunk i of the FASTA/Q
# Usage: ./split_bam_and_fasta.py in.bam in.fasta/q chunk_size 

from __future__ import print_function
import sys
import os
import csv
import pysam
from Bio import SeqIO

def build_fqi(fasta_fn, idx_fn, fmt):
	start_token = ">" if fmt == "fasta" else "@"
	idx = dict()
	with open(fasta_fn, 'r') as fasta, open(idx_fn, 'w') as handle:
		pos = 0
		line = fasta.readline()
		while len(line) > 0:
			if line[0] == start_token:
				# new record
				try:
					print("{}\t{}".format(read_name, read_pos), file=handle)
					idx[read_name] = read_pos
				except NameError:
					# first line, read_name not yet set
					pass
				read_pos = pos
				read_name = line.strip(start_token).split(" ")[0]
			pos = fasta.tell()
			line = fasta.readline()
	return idx

def load_fqi(idx_fn):
	idx = dict()
	with open(idx_fn, 'r') as handle:
		reader = csv.reader(handle, delimiter="\t")
		for row in reader:
			idx[row[0]] = int(row[1])
	return idx

bam_fn = sys.argv[1]
fasta_fn = sys.argv[2]
chunk_size = int(sys.argv[3])
fmt = "fastq" if fasta_fn[-1] == "q" else "fasta"
out_fn = bam_fn.replace(".bam", "")

idx_fn = "{}.fqi".format(fasta_fn) # not actually an fai; but close.
if not os.path.isfile(idx_fn):
	idx = build_fqi(fasta_fn, idx_fn, fmt)
else:
	idx = load_fqi(idx_fn)

with pysam.AlignmentFile(bam_fn) as bam, open(fasta_fn, 'r') as fasta:
	i = 0
	out_bam = pysam.AlignmentFile("{}.{}.bam".format(out_fn, i // chunk_size), 'wb', template=bam)
	out_fasta = open("{}.{}.{}".format(out_fn, i // chunk_size, fmt), 'w')
	for alignment in bam:
		try:
			pos = idx[alignment.query_name]
			fasta.seek(pos)
			record = next(SeqIO.parse(fasta, fmt))
			out_bam.write(alignment)
			SeqIO.write(record, out_fasta, fmt)
			i += 1	
			if i % chunk_size == 0:
				out_bam.close()
				out_fasta.close()
				out_bam = pysam.AlignmentFile("{}.{}.bam".format(out_fn, i // chunk_size), 'wb', template=bam)
				out_fasta = open("{}.{}.{}".format(out_fn, i // chunk_size, fmt), 'w')
		except KeyError:
			# unmapped read
			pass
	out_bam.close()
	out_fasta.close()
