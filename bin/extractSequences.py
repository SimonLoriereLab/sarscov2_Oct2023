#!/usr/bin/env python3
# coding: utf-8

import sys, getopt
import csv
import logging
import gzip
import os
import lzma
import io
import random
import subprocess
import re
from Bio import SeqIO

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:n:o:")
    except getopt.GetoptError:
        print('extractSequences.py -n <name file> -o <outfile> -i <fasta file>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('extractSequences.py -n <name file> -o <outfile> -i <fasta file>')
            sys.exit(0)
        elif opt in ("-i"):
            fastafile = arg
        elif opt in ("-o"):
            outfile = arg
        elif opt in ("-n"):
            namefile = arg

            
    logging.info(f'Parsing name file {namefile}')
    namelist=dict()
    with open(namefile) as file:
        tsv_file = csv.reader(file, delimiter=",")
        for line in tsv_file:
            namelist[line[0]]=1

    
    logging.info(f'Parsing sequence file {fastafile}')
    # We then parse the sequence file to print each sequence in the proper output file
    
    with gzip.open(outfile,'wt') as outf:
        with gzip.open(fastafile,'rt') as file:
            for record in SeqIO.parse(file, "fasta"):
                id=record.description
                if record.description in namelist:
                    outf.write(f'>{record.description}\n')
                    outf.write(f'{record.seq}\n')
                    
if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
