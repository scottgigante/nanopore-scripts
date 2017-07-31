# Usage: python mask_genome_variants.py reference.fa reference.masked.fa variants.vcf
from __future__ import print_function
from Bio import SeqIO
from Bio import Seq
import sys
import csv

def load_vcf(handle):
        header = next(handle)
        while header.startswith("##"):
                header = next(handle)
        header = header.strip("#").strip("\n").split("\t")
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
        return reader

def get_ambiguity_base(bases):
	bases = set(bases)
	if bases == {'A', 'C', 'G', 'T'}:
		return 'N'
	elif bases == {'A', 'C', 'G'}:
		return 'V'
	elif bases == {'A', 'C', 'T'}:
		return 'H'
	elif bases == {'A', 'G', 'T'}:
		return 'D'
	elif bases == {'C', 'G', 'T'}:
		return 'B'
	elif bases == {'A', 'C'}:
		return 'M'
	elif bases == {'G', 'T'}:
		return 'K'
	elif bases == {'A', 'G'}:
		return 'R'
	elif bases == {'C', 'T'}:
		return 'Y'
	elif bases == {'A', 'T'}:
		return 'W'
	elif bases == {'C', 'G'}:
		return 'S'
	elif bases == {'A'}:
		return 'A'
	elif bases == {'C'}:
		return 'C'
	elif bases == {'G'}:
		return 'G'
	elif bases == {'T'}:
		return 'T'
	else:
		raise Exception("Unrecognises bases: ['" + "','".join(list(bases)) + "']")

in_fn = sys.argv[1]
out_fn = sys.argv[2]
vcf_fn = sys.argv[3]

genome = dict()
with open(in_fn, 'r') as handle:
	for record in SeqIO.parse(handle, "fasta"):
		record.seq = record.seq.tomutable()
		genome[record.id] = record

with open(vcf_fn, 'r') as handle:
	vcf = load_vcf(handle)
	for row in vcf:
		bases = row['ALT'].split(',') + [row['REF']]
		
		snp = True
		for base in bases:
			if len(base) > 1:
				snp = False
		if not snp:
			# don't touch indels
			continue
		
		pos = int(row['POS']) - 1
		record = genome[row['CHROM']]
		record.seq[pos] = get_ambiguity_base(bases)
		genome[row['CHROM']] = record

contigs = list()
max_autosome=0
for contig in genome.keys():
	contigs.append(contig)
	try:
		if int(contig) > max_autosome:
			max_autosome = int(contig)
	except ValueError:
		pass
def contig_to_int(contig):
	try:
		return(int(contig))
	except ValueError:
		if contig == 'X':
			return max_autosome + 1
		elif contig == 'Y':
			return max_autosome + 2
		elif contig == 'MT':
			return max_autosome + 3
		else:
			return max_autosome + 4
def cmp_contigs(contig):
	return(contig_to_int(contig), contig)
contigs = sorted(contigs, key=cmp_contigs)

with open(out_fn, 'w') as handle:
	for contig in contigs:
		SeqIO.write(genome[contig], handle, "fasta")
