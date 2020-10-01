#!/usr/bin/env python3


'''
File last updated October 1, 2020 by Christine Foxx (cfoxx@iastate.edu)
'''

import os
import Bio
import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='''remove_reads.py takes an input directory containing de-interleaved paired-end reads in .fastq format that \
        have been classified with Kraken 2 and contain kraken2:taxid|XXXXX information in each sequence header, processes a list of \
        taxonomic IDs to exclude, and outputs filtered de-interleaved paired end reads in .fastq format.''',
        epilog='''The BioPython cookbook: http://biopython.org/DIST/docs/tutorial/Tutorial.html$sec367'''
    )

    parser.add_argument('--input-path', dest='path', required=True,
                        help='Directory containing paired-end *.fastq file(s) labeled that contain kraken2:taxid|XXXXX information \
                        in each sequence header')
    parser.add_argument('--exclude-taxid', dest='taxid', required=True, default='9606',
                        help='List of Kraken 2 or NCBI taxonomic ID-compatible codes (separated by spaces) to exclude from each \
                        *.fastq file(s) [default: '9606', human]')

    for file in os.listdir(path):
        if file.endswith(".fastq"):
            in_file = path + "/" + file
            out_file = input.replace(".fastq", "_filtered.fastq")

            # Writing out filtered files:
            out_handle = open(out_file, "w")
            out_short = os.path.basename(out_handle)
            sys.stdout.write("Filtering sequences from %s..." % file)
            sys.stdout.flush()
            count = 0

            with open(in_file) as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    if taxid in title:
                        count += 1
                        sys.stdout.write("Found %i records matching to %s" % (count, taxid)
                    elif taxid" not in title:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                    else:
                        raise ValueError("Couldn't find any sequences matching to taxid %i in %s!" % (taxid, in_handle)
            out_handle.close()
        exit(0)
    
if __name__ == "__main__":
    main()
