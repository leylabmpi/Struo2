#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import gzip
import argparse
import logging

desc = 'Simple check of gene input'
epi = """DESCRIPTION:
Check of overlapping UUIDs for fasta files & metadata table.
Also, check that required data is provided in the metadata table.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('prot_fasta', metavar='prot_fasta', type=str,
                    help='Protein sequence fasta')
parser.add_argument('metadata', metavar='metadata', type=str,
                    help='gene metadata')
parser.add_argument('-n', '--nuc-fasta', type=str, default=None,
                    help='Nucleotide sequence fasta (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def get_open(infile):
    if infile.endswith('.gz'):
        _open = lambda x: gzip.open(x, 'rb')
    else:
        _open = lambda x: open(x, 'r')
    return _open

def read_fasta(infile):
    logging.info('Reading file: {}'.format(infile))
    _open = get_open(infile)
    seqs = []
    with _open(infile) as inF:
        for line in inF:
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            if line.startswith('>'):
                seqs.append(line.lstrip('>').rstrip())
    return set(seqs)

def read_meta(infile):
    logging.info('Reading file: {}'.format(infile))
    _open = get_open(infile)
    meta = {}
    header = {}
    req_cols = ['seq_uuid', 'seq_orig_name', 'domain', 'phylum',
                'class', 'order', 'family', 'genus', 'species',
                'taxid', 'genome_name', 'genome_length_bp']
    non_empty_cols = ['seq_uuid', 'seq_orig_name', 'genus', 'species']
    entry_cnt = 0
    with _open(infile) as inF:
        for i,line in enumerate(inF):
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            line = line.rstrip().split('\t')
            if line == '':
                continue
            # header
            if i == 0:
                header = {x:ii for ii,x in enumerate(line)}
                missing = []
                for x in req_cols:                    
                    try:
                        _ = header[x]
                    except KeyError:
                        missing.append(x)
                if len(missing) > 0:
                    msg = 'Missing columns in metadata table: {}'
                    raise ValueError(msg.format(','.join(missing)))
                else:
                    logging.info('The metadata table has all required columns')
                continue
            # body
            entry_cnt += 1
            for col in non_empty_cols:
                if line[header[col]] == '':
                    msg = 'Line {}: Column "{}" cannot be empty'
                    raise ValueError(msg.format(i+1, col))
            seq_uuid = line[header['seq_uuid']] 
            seq_orig_name = line[header['seq_orig_name']]
            try:
                _ = meta[seq_uuid]
                raise ValueError('Value duplicated: {}'.format(seq_uuid))
            except KeyError:
                pass
            meta[seq_uuid] = seq_orig_name 
            
    return meta
            
def main(args):
    # aa fasta
    aa_seqs = read_fasta(args.prot_fasta)
    # nuc fasta
    if args.nuc_fasta is not None:        
        nuc_seqs = read_fasta(args.nuc_fasta)
        if not len(aa_seqs) == len(nuc_seqs):
            msg = 'WARNIGN: prot. & nuc. seqs differ in length!'
            logging.warning(msg)
            # compare fasta files
            logging.info('Comparing fasta files...')
            ## just in prot.
            just_aa = list(aa_seqs - nuc_seqs)
            if len(just_aa) > 0:
                just_aa = '\n  '.join(just_aa)
                print('Genes just in the prot. fasta:\n  {}'.format(just_aa))
            ## just in prot.
            just_nuc = list(nuc_seqs - aa_seqs)
            if len(just_nuc) > 0:
                just_nuc = '\n  '.join(just_nuc)
                print('Genes just in the nuc. fasta:{}  \n'.format(just_nuc))
            raise ValueError('Exiting due to mismatches')
    # metadata
    meta = read_meta(args.metadata)
    ## comparing to fasta
    if len(meta.keys()) != len(aa_seqs):
        msg = 'WARNIGN: No. of metadata entries does not match no. of prot. seqs'
        logging.warning(msg)
        ### just in protein
        just_aa = list(aa_seqs - set(meta.keys()))
        if len(just_aa) > 0:
            just_aa = '\n  '.join(just_aa)
            print('Genes just in the prot. fasta:\n  {}'.format(just_aa))
        just_txt = list(set(meta.keys()) - aa_seqs)
        if len(just_txt) > 0:
            msg = 'Genes just in the metadata table:\n  {}'
            print(msg.format('\n  '.join(just_txt)))
            just_txt = '\n  '.join([meta[x] for x in just_txt])
            print('Genes just in the metadata table (gene_orig_name):\n  {}'.format(just_txt))
        raise ValueError('Exiting due to mismatches')
        
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
