# Usage: python check_alignment_identity.py in1.bam in2.bam [...] out.tsv

import re
import pysam
import sys
import multiprocessing
import itertools
import traceback

def idmm(read):
    read_length = read.reference_length
    try:
        md = [t[1] for t in read.get_tags() if t[0] == 'MD'][0]
    except IndexError:
        # no MD tag - skip
        return 0, 0, 0, 0
    cigar = read.cigartuples

    deletions = sum([len(i) for i in re.findall("\^([A|C|G|T]*)", md)])
    mismatches = len(re.findall("A|C|G|T", re.sub("\^[A|G|C|T]*", "", md)))
    insertions = sum([c[1] for c in cigar if c[0] == 1])
    # insertions = sum([int(i) for i in re.findall("I([0-9]*)",read.cigar)])

    return read_length, insertions, deletions, mismatches

def analyse_bam(filename):
    try:
        aligned_len = 0.0
        insertion = 0
        deletion = 0
        mismatch = 0
        with pysam.AlignmentFile(filename) as bamfile:
            for read in bamfile:
                result = idmm(read)
                aligned_len += result[0]
                insertion += result[1]
                deletion += result[2]
                mismatch += result[3]
        insertion = insertion/aligned_len
        deletion = deletion/aligned_len
        mismatch = mismatch/aligned_len
        return filename, aligned_len, insertion, deletion, mismatch, insertion + deletion + mismatch
    except Exception:
        traceback.print_exc()
        print()
        raise

outfilename = sys.argv[-1]
filenames = sys.argv[1:-1]
if re.search("\.(bam|sam)$", outfilename):
    print "Error: output file cannot ends with .bam or .sam"
    print "Usage: python check_alignment_identity.py in1.bam in2.bam [...] out.tsv"
    exit()

p = multiprocessing.Pool(len(filenames))
with open(outfilename, 'w') as outfile:
    outfile.write("filename\tbases\tinsertion\tdeletion\tmismatch\ttotal\n")
    result = p.map(analyse_bam, filenames)
    for r in result:
        outfile.write("\t".join([str(i) for i in r]) + "\n")
p.close()
p.join()

