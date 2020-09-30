#!/usr/bin/env python3


'''
Written: September 30, 2020 by Christine Foxx (cfoxx@iastate.edu)
Checked output cleaned_R1.fastq and cleaned_R2.fastq for any matches to human reads (taxid:9606)
'''

import os
import Bio
import sys
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='''remove_reads.py takes an input directory with de-interleaved paired-end *.fastq file(s) and outputs 2 files containing de-interleaved reads in the format 
        (*_R1.fastq containing forward reads and *_R2.fastq containing reverse reads).''',
        epilog='''The BioPython cookbook: http://biopython.org/DIST/docs/tutorial/Tutorial.html$sec367'''
    )
    exit(0)

if __name__ == "__main__":
    main()
