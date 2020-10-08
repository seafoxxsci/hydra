#!/usr/bin/env python3


"""
File last updated October 8, 2020 by Christine Foxx (cfoxx@iastate.edu)
"""

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
        epilog='''The BioPython cookbook: http://biopython.org/DIST/docs/tutorial/Tutorial.html'''
    )
    parser.add_argument('--input-path', dest='path', required=True,
                        help='Directory containing paired-end *.fastq file(s) labeled that contain kraken2:taxid|XXXXX information \
                        in each sequence header')
    parser.add_argument('--exclude-taxid', dest='taxid', required=True, nargs='+',
                        help='List of Kraken 2 or NCBI taxonomic ID-compatible codes (separated by spaces) to exclude from each \
                        *.fastq file(s)')
    args = parser.parse_args()

    key = '_classified_'
    dir = args.path

    for file in os.listdir(dir):
        if file.endswith(".fastq") and (key in file):
            in_file = dir + "/" + file
            out_file = in_file.replace(".fastq", "_filtered.fastq")
            taxa = ["kraken:taxid|" + x for x in args.taxid]  # Make this a list, loop through list elements

            # Writing out filtered files:
            out_handle = open(out_file, "w")
            sys.stdout.write("Filtering sequences from %s..." % file)
            sys.stdout.flush()
            count = 0

            with open(in_file) as in_handle:
                 for title, seq, qual in FastqGeneralIterator(in_handle):
                     for taxon in taxa:
                        if taxon in title:
                            count += 1
                            sys.stdout.write("Found %i records matching to %s\r" % (count, taxon))
                            sys.stdout.flush()
                        elif taxon not in title:
                            out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        else:
                            raise ValueError("Are you sure your files have been classified by Kraken 2?")
            out_handle.close()
    exit(0)


if __name__ == "__main__":
    main()
