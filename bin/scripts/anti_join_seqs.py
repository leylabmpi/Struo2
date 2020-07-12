#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import uuid
import argparse
import logging

desc = 'Filtering out the intersection of two fasta files'
epi = """DESCRIPTION:
Taking the setdiff of the first fasta. 

Output written to STDOUT.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('fasta1', metavar='fasta1', type=str,
                    help='The first fasta file')
parser.add_argument('fasta2', metavar='fasta2', type=str,
                    help='The second fasta file')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def make_index(fasta):
    regex = re.compile(r' .+')
    if fasta.endswith('.gz'):
        _openR = lambda x: gzip.open(x, 'rb')
    else:
        _openR = lambda x: open(x, 'r')        
    
    idx = {}
    with _openR(fasta) as inF:
        for line in inF:
            if fasta.endswith('.gz'):
                line = line.decode('utf8')            
            if line.startswith('>'):
                line = line.lstrip('>').rstrip()
                idx[regex.sub('', line)] = 0  
    return set(idx.keys())

def filter_fasta(fasta, idx):
    if fasta.endswith('.gz'):
        _openR = lambda x: gzip.open(x, 'rb')
    else:
        _openR = lambda x: open(x, 'r')        
    

    found = {}
    hit = False
    with _openR(fasta) as inF:
        for line in inF:
            line = line.rstrip()
            if fasta.endswith('.gz'):
                line = line.decode('utf8')
            if line.startswith('>'):
                line = line.lstrip('>')
                # filter is already seen
                try:
                    found[line]
                    continue
                except KeyError:
                    pass
                # is seq in index?
                try:                    
                    found[line] = idx[line]
                    hit = True
                except KeyError:
                    hit = False
                    continue
                print(seq_name = '>' + idx[line])
            else:
                if hit:
                    print(line)
                        
    logging.info('File written: {}'.format(output))
    logging.info('Number of seqs written: {}'.format(len(found.keys())))
        
def main(args):    
    # creating the seq header index
    seq_idx = make_index(args.fasta1) & make_index(args.fasta2)
    
    # filtering the fasta files
    filter_fasta(args.fasta1, seq_idx)
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
