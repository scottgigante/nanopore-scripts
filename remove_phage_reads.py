#!/usr/bin/python

# Usage: python remove_phage_reads.py in.REF_PHAGE.bam all_reads.fasta path/to/phage/

"""
A script for separating phage from other reads
"""

import pysam
import sys
import os
import shutil
import multiprocessing
import functools
import itertools
from Bio import SeqIO

def loadRef(fasta):
    """
    Create a dictionary linking read names to fast5 file paths

    :param fasta: Filename of the fasta file created by poretools

    :returns: Dictionary linking read names to fast5 file paths
    """
    refs = dict()
    readsPath = '/'.join(fasta.split('/')[:-1])
    handle = open(fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        #print "%s with %d bases\n" % (record.id, len(record.seq))
        fast5Path = record.description.split(' ')[-1]
        refs[record.id] = os.path.abspath(fast5Path)
    handle.close()
    return refs

def move(src, dst):
    try:
        shutil.move(src, dst)
    except Exception as e:
        if os.path.isfile(src):
            print str(e)
        pass

_, bamfile, fastafile, out_path = sys.argv
# the bamfile should be only those reads which align to the phage. Align to a chimeric chromosome and use bamtools split -in in.bam -reference
# the fastafile should be generated using poretools

refs = loadRef(fastafile)

try:
    os.mkdir(out_path)
except OSError:
    pass

with pysam.AlignmentFile(bamfile) as bam:
    pool = multiprocessing.Pool()
    pool.map(functools.partial(move, dst=out_path), [refs[read.query_name] for read in bam if 'XA' not in itertools.izip(*read.get_tags()).next()])
    pool.close()
    pool.join()


