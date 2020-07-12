#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import uuid
import argparse
import logging

desc = 'Adding user-provided sequences'
epi = """DESCRIPTION:
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('fasta', metavar='fasta', type=str, nargs='+',
                    help='Fasta files (nuc, then prot)')
parser.add_argument('--in-fna', type=str, default='genes.fna',
                    help='Nucleotide output')
parser.add_argument('--in-faa', type=str, default='genes.faa',
                    help='Amino acid output')
parser.add_argument('--in-txt', type=str, default='genes.txt',
                    help='Names index output')
parser.add_argument('--out-fna', type=str, default='wUser_genes.fna',
                    help='Nucleotide output')
parser.add_argument('--out-faa', type=str, default='wUser_genes.faa',
                    help='Amino acid output')
parser.add_argument('--out-txt', type=str, default='wUser_genes.txt',
                    help='Names index output')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def read_fasta(infile):
    seqs = {}
    seq_name = ''
    with open(infile) as inF:
        for line in inF:
            line = line.rstrip()
            if line.startswith('>'):
                seq_name = line.lstrip('>')
            else:
                try:
                    seqs[seq_name] += line
                except KeyError:
                    seqs[seq_name] = line
    return seqs            

def make_index(input_fna, input_faa):
    # loading fasta files
    if input_fna.lower() != 'skip':
        fna = read_fasta(input_fna)
    else:
        fna = {}
    if input_faa.lower() != 'skip':
        faa = read_fasta(input_faa)
    else:
        faa = {}
    # union of names
    seq_names = set(fna.keys()) & set(faa.keys())
    seq_idx = {x:str(uuid.uuid4()) for x in seq_idx}
    # return
    return fna, faa, seq_idx

def seq_cat(seqs, in_fasta, out_fasta):
    with open(in_fasta) as inF, open(out_fasta, 'w') as outF:
        for line in inF:
            outF.write(line)
        for seqid,seq in seqs.items():
            outF.write('>' + seqid + '\n' + seq + '\n')
    logging.info('File written: {}'.format(out_fasta))

def names_cat(seq_idx, in_txt, out_txt):
    with open(in_txt) as inF, open(out_txt, 'w') as outF:
        for line in inF:
            outF.write(line)
        for uuid,seqid in seq_idx.items():
            outF.write(uuid + '\t' + seqid + '\n')
    logging.info('File written: {}'.format(out_txt))
            
def main(args):
    # getting overlap of user-provided nuc & prot gene names
    fna,faa,seq_idx = make_index(args.fasta[0], args.fasta[1])

    # combining sequences
    seq_cat(fna, args.in_fna, args.out_fna)
    seq_cat(faa, args.in_faa, args.out_faa)
    names_cat(seq_idx, args.in_txt, args.out_txt)
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
