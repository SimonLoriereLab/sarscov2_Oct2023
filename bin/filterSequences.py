#! /usr/bin/env python3
# coding: utf-8
import sys, getopt
import gzip
import csv
import logging
import os
import subprocess
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor

#
# This script takes an input fasta file, and outputs a fasta file without sequences having more than 1% N
# Moreover it renames the sequences : spaces are replaced with --
#

def cleanName(name):
    n=name.replace(' ','--')
    n=n.split("|")[0]
    return n

def getNames(f):
    names = {}
    with open(f) as file:
        for line in file:
            names[line.rstrip()] = 1
    return names

def compress(file):
    print(f'START COMPRESS {file}')
    subprocess.call(f'gzip {file}', shell=True)
    print(f'END COMPRESS {file}')

def main(argv):
    split=1000
    cpus=3
    try:
        opts, args = getopt.getopt(argv,"hs:o:S:B:b:C:c:t:")
    except getopt.GetoptError:
        print('filterSequences.py -s <sequences fasta> -o <output sequences> -S <Split by number of sequences> -B <BA.2.86 names file> -b <BA.2.86 output> -t <cpus>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('filterSequences.py -s <sequences fasta> -o <output sequences> -S <Split by the number of sequences> -B <BA.2.86 names file> -b <BA.2.86 output> -t <cpus>') 
            sys.exit(1)
        elif opt in ("-s"):
            fastafile = arg
        elif opt in ("-o"):
            fastaoutfile = arg
        elif opt in ("-S"):
            split=int(arg)
        elif opt in ("-B"):
            ba286 = arg
        elif opt in ("-b"):
            ba286output = arg
        elif opt in ("-C"):
            context = arg
        elif opt in ("-c"):
            contextoutput = arg
        elif opt in ("-t"):
            cpus = int(arg)
            

    logging.info(f'Input Fasta file : {fastafile}')
    logging.info(f'Output Fasta file : {fastaoutfile}')
    logging.info(f'Split by : {split}')
    logging.info(f'BA.2.86 name file : {ba286}')
    logging.info(f'BA.2.86 out file : {ba286output}')
    logging.info(f'Context name file : {context}')
    logging.info(f'Context out file : {contextoutput}')

    if cpus < 3 :
        cpus=3
    
    executor =  ThreadPoolExecutor(max_workers=cpus)
    ba286names = getNames(ba286)
    contextnames = getNames(context)
    nbatch=0
    nseq=0
    outfile=None

    executor =  ThreadPoolExecutor(max_workers=cpus)
    
    # We parse the fasta file again to remove sequences
    if fastafile == '-':
        file = sys.stdin
    else:
        file = open(fastafile)
    outfileba286 = gzip.open(f'{ba286output}', 'wt')
    outfilecontext = gzip.open(f'{contextoutput}', 'wt')
    filename=f'{fastaoutfile}_{nbatch}.fasta'
    futures_list=[]
    for record in SeqIO.parse(file, "fasta"):
        if len(record.seq)>29000 and float(record.seq.count('N')+record.seq.count('n'))/float(len(record.seq))<0.01:
            if nseq%split==0 :
                if outfile!=None:
                    outfile.close()
                    futures = executor.submit(compress, filename)
                    futures_list.append(futures)
                filename=f'{fastaoutfile}_{nbatch}.fasta'
                outfile = open(filename, 'w')
                nbatch+=1
            cname = cleanName(record.description)
            if cname in ba286names:
                print(f'>{cname}\n{record.seq}',file=outfileba286)
            elif cname in contextnames :
                print(f'>{cname}\n{record.seq}',file=outfilecontext)
            else:
                print(f'>{cname}\n{record.seq}',file=outfile)
            nseq+=1

    if outfile!=None:
        outfile.close()
        futures = executor.submit(compress, filename)
        futures_list.append(futures)
    file.close()
    outfileba286.close()
    outfilecontext.close()
    for future in futures_list:
        try:
            future.result()
        except Exception:
            results.append(None)

if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    main(sys.argv[1:])
