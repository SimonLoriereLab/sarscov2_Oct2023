#!/usr/bin/env python3
# coding: utf-8

import sys, getopt
import csv
import logging
import os
from Bio import SeqIO

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:s:o:m:c:p:")
    except getopt.GetoptError:
        print('renameSeq.py -i <emergen convert file> -s <emergen seq file> -o <output seq file>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('renameSeq.py -i <emergen convert file> -s <emergen seq file> -o <output seq file>')
            sys.exit(0)
        elif opt in ("-i"):
            inputfile = arg
        elif opt in ("-s"):
            seqfile = arg
        elif opt in ("-o"):
            outseq = arg
            
    logging.info(f'Emergen file is {inputfile}')
    logging.info(f'Emergen sequence file is {seqfile}')
    logging.info(f'Output sequence file is {outseq}')

    # We parse the conversion file
    rename=dict()
    with open(inputfile) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            rename[line[0]]=line[1]

    with open(outseq, 'w') as outfile:
        # We parse the fasta file again to rename the sequences
        with open(seqfile) as file:
            for record in SeqIO.parse(file, "fasta"):
                print(f'>{rename.get(record.description)}',file=outfile)
                print(f'{record.seq}',file=outfile)

if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
