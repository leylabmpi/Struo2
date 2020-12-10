#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import gzip
import bz2
import argparse
import logging
import pickle
from pickle import UnpicklingError

desc = 'Transferring UniRef annotations by cluster cutoff'
epi = """DESCRIPTION:
Using an index of how UniRef50 clusters map to UniRef90 clusters.
Renaming annotations with the other UniRef cluster level
(eg., UniRef90_Q8WZ42-5 => UniRef50_Q8WZ42-5).
Due to the multi-mapping when going from UniRef50 to UniRef90,
it makes much more sense to go from UniRef90 to UniRef50.
For multi-mappings, one mapping will be selected.

The input file format is determined by the input file extension.
Output is written to STDOUT.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('index_file', metavar='index_file', type=str,
                    help='UniRef 50 <=> 90 index')
parser.add_argument('--in-nuc', type=str, default='', 
                    help='Input nucleotide fasta file path (default: %(default)s)')
parser.add_argument('--in-prot', type=str, default='', 
                    help='Input amino acid fasta file path (default: %(default)s)')
parser.add_argument('--in-tsv', type=str, default='', 
                    help='Input gene metadata table path (default: %(default)s)')
parser.add_argument('--out-nuc', type=str, default='', 
                    help='Output nucleotide fasta file path (default: %(default)s)')
parser.add_argument('--out-prot', type=str, default='', 
                    help='Output amino acid fasta file path (default: %(default)s)')
parser.add_argument('--out-tsv', type=str, default='', 
                    help='Output gene metadata table path (default: %(default)s)')
parser.add_argument('-d', '--direction', type=str, default='90=>50',
                    choices = ['90=>50', '50=>90'],
                    help='Changing annotations from X to Y (default: %(default)s)')
parser.add_argument('-p', '--pickle-idx', type=str, default='',
                    help='Write the index as a pickled dict for faster loading (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def _open(infile, mode='rb'):
    """
    Openning of input, regardless of compression
    """
    if infile.endswith('.bz2'):
        return bz2.open(infile, mode)
    elif infile.endswith('.gz'):
        return gzip.open(infile, mode)
    else:
        return open(infile)

def _decode(line, infile):
    """
    Decoding input, depending on the file extension
    """
    if os.path.isfile(infile) and (infile.endswith('.gz') or infile.endswith('.bz2')):
        line = line.decode('utf-8')
    return line

def read_index(infile, direction):
    logging.info('Loading file: {}'.format(infile))
    # if pickle
    try:
        idx = pickle.load(open(infile, 'rb'))
        return idx
    except UnpicklingError:
        pass
    
    idx = {}
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            line = _decode(line, infile)
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            if len(line) < 2:
                msg = 'Line{}: <2 values'
                raise ValueError(msg.format(i))
            # assuming UniRef50<tabl>UniRef90
            if direction == '50=>90':
                idx[line[0]] = line[1]
            else:
                idx[line[1]] = line[0]
    logging.info('  No. of index keys: {}'.format(len(idx.keys())))
    return idx

def which_ext(infile):
    infile = infile.rstrip('.gz')
    for x in ['.fasta', '.fna', '.faa', '.fa']:
        if infile.endswith(x):
            return 'fasta'
    for x in ['.txt', '.tsv']:
        if infile.endswith(x):
            return 'txt'
    msg = 'Cannot determine file type from extension: {}'
    raise IOError(msg.format(infile))

def make_dir(outfile):
    D = os.path.split(outfile)[0]
    if not os.path.isdir(D):
        os.makedirs(D)

def rename_fasta(infile, outfile, idx):
    logging.info('Processing file: {}'.format(infile))
    make_dir(outfile)
    to_write = False
    stats = {'in_index' : 0, 'not_in_index' : 0}
    with _open(infile) as inF, open(outfile, 'w') as outF:
        for line in inF:
            line = _decode(line, infile)
            # header
            if line.startswith('>'):
                line = line.lstrip('>').split('|')
                try:
                    line[0] = idx[line[0]]
                    to_write = True
                    stats['in_index'] += 1
                except KeyError:
                    msg = 'Cannot find "{}" in index'
                    logging.warning(msg.format(line[0]))
                    to_write = False
                    stats['not_in_index'] += 1                    
                if to_write is True:
                    line = '>' + '|'.join(line)
            # body
            if to_write is True:
                outF.write(line)
    # status
    msg = '  No. of genes found in UniRef50<=>90 index: {}'
    logging.info(msg.format(stats['in_index']))
    msg = '  No. of genes NOT found in UniRef50<=>90 index: {}'
    logging.info(msg.format(stats['not_in_index']))
    if stats['in_index'] == 0:
        raise ValueError('No genes were present in the UniRef50<=>90 index!')
                    
def rename_txt(infile, outfile, idx):
    logging.info('Processing file: {}'.format(infile))
    make_dir(outfile)    
    header = {}
    to_write = False
    stats = {'in_index' : 0, 'not_in_index' : 0}
    with _open(infile) as inF, open(outfile, 'w') as outF:
        for i,line in enumerate(inF):
            line = _decode(line, infile)
            line = line.rstrip().split('\t')
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
            else:
                try:
                    line[header['annotation']] = idx[line[header['annotation']]]
                    to_write = True
                    stats['in_index'] += 1       
                except KeyError:
                    to_write = False
                    stats['not_in_index'] += 1
            if i == 0 or to_write is True:
                outF.write('\t'.join(line) + '\n')
    # status
    msg = '  No. of genes found in UniRef50<=>90 index: {}'
    logging.info(msg.format(stats['in_index']))
    msg = '  No. of genes NOT found in UniRef50<=>90 index: {}'
    logging.info(msg.format(stats['not_in_index']))
                
def pickle_idx(idx, outfile):
    logging.info('Pickling index to {}'.format(outfile))
    with open(outfile, 'wb') as outF:
        pickle.dump(idx, outF)
    logging.info('  File pickled. Exiting')
    sys.exit()
            
def main(args):
    # loading UniRef cluster index
    idx = read_index(args.index_file, args.direction)
    # pickle
    if args.pickle_idx != '':
        pickle_idx(idx, args.pickle_idx)
    # renaming input
    for inpath,outpath in zip([args.in_nuc, args.in_prot, args.in_tsv],
                              [args.out_nuc, args.out_prot, args.out_tsv]):
        if which_ext(inpath) == 'fasta':
            rename_fasta(inpath, outpath, idx)
        elif which_ext(inpath) == 'txt':
            rename_txt(inpath, outpath, idx)
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
