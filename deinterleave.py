#!/usr/bin/env python3


__author__ = "Christine Foxx"
__copyright__ = "Last updated September 30, 2020"
__credits__ = ["Christine Foxx"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Christine Foxx"
__email__ = "cfoxx@iastate.edu"

'''
Checked _R1.fastq and _R2.fastq sequence length == 0.5(file.fastq)
'''

import os
import Bio
import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='''deinterleave.py takes an input directory with paired-end *.fastq file(s) and outputs 2 files containing de-interleaved reads in the format 
        (*_R1.fastq containing forward reads and *_R2.fastq containing reverse reads).''',
        epilog='''The BioPython cookbook: http://biopython.org/DIST/docs/tutorial/Tutorial.html$sec367'''
    )
    
    parser.add_argument('--input-path', dest='path', required=True,
                        help='Directory containing paired-end *.fastq file(s) labeled with conventional Illumina headers')
                        
    for file in os.listdir(path):
        if file.endswith(".fastq"):
            input = path + "/" + file
            fwd_out = input.replace(".fastq", "_R1.fastq")
            rev_out = input.replace(".fastq", "_R2.fastq")

            # Writing out forward and reverse reads as separate files:
            fwd_handle = open(fwd_out, "w")
            rev_handle = open(rev_out, "w")
            fwd_short = os.path.basename(fwd_out)
            rev_short = os.path.basename(rev_out)
            sys.stdout.write("De-interleaving sequences from %s..." % file)
            sys.stdout.flush()
            fwd_count = 0
            rev_count = 0

            with open(input) as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    if "1:N:0" or "/1" in title:
                        fwd_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        fwd_count += 1
                        sys.stdout.write("Saved %i records to %s" % (fwd_count, fwd_short))
                        sys.stdout.flush()
                    elif "2:N:0" or "/2" in title:
                        rev_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        rev_count += 1
                        sys.stdout.write("Saved %i records to %s" % (rev_count, rev_short))
                        sys.stdout.flush()
                    else:
                        sys.stdout.write("Couldn't de-interleave sequence %s! Are you sure you have paired-end reads in this file?" % title)
                        sys.stdout.flush()
                        break
            rev_handle.close()
            fwd_handle.close()
    exit(0)

if __name__ == "__main__":
    main()
