#!/usr/bin/env python3


"""
Last updated September 30, 2020 by Christine Foxx (cfoxx@iastate.edu)
Checked _R1.fastq and _R2.fastq sequence length == 0.5(file.fastq)
"""

import os
import Bio
import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='''deinterleave.py takes an input directory with paired-end *.fastq file(s) and outputs 2 files \
        containing de-interleaved reads in the format (*_R1.fastq containing forward reads and *_R2.fastq containing \
        reverse reads).''',
        epilog='''The BioPython cookbook: http://biopython.org/DIST/docs/tutorial/Tutorial.html$sec367'''
    )
    
    parser.add_argument('--input-path', dest='path', required=True,
                        help="Directory containing paired-end *.fastq file(s) labeled with conventional Illumina \
                        headers ('1:N:0' or '2:N:0') or conventional Sanger headers ('/1' or '/2') for each sequence")
    parser.set_defaults(append=False)
    args = parser.parse_args()

    for file in os.listdir(args.path):
        if file.endswith(".fastq"):
            path = args.path
            input_file = path + "/" + file
            fwd_out = input_file.replace(".fastq", "_R1.fastq")
            rev_out = input_file.replace(".fastq", "_R2.fastq")

            # Writing out forward and reverse reads as separate files:
            fwd_handle = open(fwd_out, "w")
            rev_handle = open(rev_out, "w")
            fwd_short = os.path.basename(fwd_out)
            rev_short = os.path.basename(rev_out)
            sys.stdout.write("De-interleaving sequences from %s...\n" % file)
            sys.stdout.flush()
            fwd_count = 0
            rev_count = 0

            with open(input) as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    if "1:N:0" or "/1" in title:
                        fwd_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        fwd_count += 1
                        sys.stdout.write("Saved %i records to %s\r" % (fwd_count, file))
                        sys.stdout.flush()
                    elif "2:N:0" or "/2" in title:
                        rev_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        rev_count += 1
                    else:
                        raise ValueError("Couldn't de-interleave sequence %s! Are you sure you have paired-end reads \
                        in %s?" % (title, file))
            rev_handle.close()
            fwd_handle.close()
            sys.stdout.write("Saved %i records to %s\r" % (fwd_count, fwd_short))
            sys.stdout.write("Saved %i records to %s\r" % (rev_count, rev_short))
    exit(0)


if __name__ == "__main__":
    main()
