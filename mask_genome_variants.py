# Usage: python mask_genome_variants.py reference.fa reference.masked.fa variants.vcf
from __future__ import print_function
from Bio import SeqIO
from Bio import Seq
import sys
import csv
import argparse

def load_vcf(handle):
        header = next(handle)
        while header.startswith("##"):
                header = next(handle)
        header = header.strip("#").strip("\n").split("\t")
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
        return reader

__ambiguity_codes = {
	'N' : {'A', 'C', 'G', 'T'},
	'V' : {'A', 'C', 'G'},
	'H' : {'A', 'C', 'T'},
	'D' : {'A', 'G', 'T'},
	'B' : {'C', 'G', 'T'},
	'M' : {'A', 'C'},
	'K' : {'G', 'T'},
	'R' : {'A', 'G'},
	'Y' : {'C', 'T'},
	'W' : {'A', 'T'},
	'S' : {'C', 'G'},
	'A' : {'A'},
	'C' : {'C'},
	'G' : {'G'},
	'T' : {'T'}
}

def get_ambiguity_base(bases):
	bases = set(bases)
	for ambiguity_base, possible_bases in __ambiguity_codes.items():
		if bases == possible_bases:
			return ambiguity_base
	raise Exception("Unrecognises bases: ['" + "','".join(list(bases)) + "']")

def disambiguate_base(base):
	base = str(base)
	return list(__ambiguity_codes[base])

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True, help='filename of input genome fasta')
	parser.add_argument('-o', '--output', required=True, help='filename of output masked genome fasta')
	parser.add_argument('-v', '--vcf', required=True, help='filename of vcf listing snps')
	parser.add_argument('--alternate-only', default=False, action='store_true', help='replace reference base with alternate allele, rather than retaining both (Default: false)')
	parser.add_argument('--update', default=False, action='store_true', help='Update an already masked genome with another VCF file (Default: false)')
	return parser.parse_args()

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

args = parse_args()

genome = dict()
with open(args.input, 'r') as handle:
	for record in SeqIO.parse(handle, "fasta"):
		record.seq = record.seq.tomutable()
		genome[record.id] = record

with open(args.vcf, 'r') as handle:
	vcf = load_vcf(handle)
	for row in vcf:
		bases = row['ALT'].split(',')
		if not args.alternate_only:
			bases = bases + [row['REF']]
		snp = True
		for base in bases:
			if len(base) > 1:
				snp = False
		if not snp:
			# don't touch indels
			continue
		
		pos = int(row['POS']) - 1
		record = genome[row['CHROM']]
		if args.update:
			bases = bases + disambiguate_base(record.seq[pos])
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

contigs = sorted(contigs, key=cmp_contigs)

with open(args.output, 'w') as handle:
	for contig in contigs:
		SeqIO.write(genome[contig], handle, "fasta")
