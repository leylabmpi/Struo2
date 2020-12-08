#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import gzip
import uuid
import argparse
import logging

desc = 'Filtering two fasta files down to just intersection'
epi = """DESCRIPTION:
Filtering 2 fasta files down to the intersection of their sequencines.
The sequence headers must perfectly match.
If any duplicate headers, only the first will be selected.

Output columns:
* seq UUID
* seq original name
* domain
* phylum
* class
* order
* family
* genus
* species
* taxid
* genome ID
* genome length

Output written to STDOUT.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('fasta1', metavar='fasta1', type=str,
                    help='The first fasta file')
parser.add_argument('fasta2', metavar='fasta2', type=str,
                    help='The second fasta file')
parser.add_argument('fasta1_output', metavar='fasta1_output', type=str,
                    help='Name of the output fasta1 file')
parser.add_argument('fasta2_output', metavar='fasta2_output', type=str,
                    help='Name of the output fasta2 file')
parser.add_argument('--taxonomy', type=str, default='',
                    help='genome taxonomy')
parser.add_argument('--taxID', type=str, default='',
                    help='genome taxonomy')
parser.add_argument('--accession', type=str, default='',
                    help='genome accession')
parser.add_argument('--genome-file', type=str, default='',
                    help='genome fasta file (to get genome length)')
parser.add_argument('--gzip', action='store_true', default=False,
                    help='gzip output')
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

def format_taxonomy(tax, taxID):
    """Formatting taxonomy string
    """
    logging.info('Taxonomy string provided: {}'.format(tax))
    logging.info('TaxID provided: {}'.format(taxID))
    
    try:
        taxID = int(float(taxID.strip()))
    except ValueError:
        msg = 'WARNING: taxID "{}" is not an integer!'
        logging.warning(msg.format(taxID))        
        #raise ValueError(msg)
    tax = [re.sub('[^A-Za-z0-9-_;]+', '_', x) for x in tax.split(';')]
    regex_d = re.compile(r'[Dd]__.+')
    regex_p = re.compile(r'[Pp]__.+')
    regex_c = re.compile(r'[Cc]__.+')
    regex_o = re.compile(r'[Oo]__.+')
    regex_f = re.compile(r'[Ff]__.+')    
    regex_g = re.compile(r'[Gg]__.+')
    regex_s = re.compile(r'[Ss]__.+')

    domain = 'd__unclassified'
    phylum = 'p__unclassified'
    cls = 'c__unclassified'
    order = 'o__unclassified'
    family = 'f__unclassified'
    genus = 'g__unclassified'
    species = 's__unclassified'
    for lev in tax:
        if regex_d.match(lev):
            domain = lev
        elif regex_p.match(lev):
            phylum = lev
        if regex_c.match(lev):
            cls = lev
        elif regex_o.match(lev):
            order = lev
        if regex_f.match(lev):
            family = lev
        if regex_g.match(lev):
            genus = lev
        elif regex_s.match(lev):
            species = lev
            
    tax = [domain, phylum, cls, order, family, genus, species, str(taxID)]
    logging.info('Converted taxonomy string to {}'.format(';'.join(tax)))
    return tax

def filter_fasta(fasta, idx, output, gzip_out=False):
    """
    Filtering fasta to just those in idx
    """
    if fasta.endswith('.gz'):
        _openR = lambda x: gzip.open(x, 'rb')
    else:
        _openR = lambda x: open(x, 'r')        
    
    if gzip_out is True:
        _openW = lambda x: gzip.open(x, 'wb')
    else:
        _openW = lambda x: open(x, 'w')

    found = {}
    hit = False
    regex = re.compile(r' .+')
    with _openR(fasta) as inF, _openW(output) as outF:
        for line in inF:
            if fasta.endswith('.gz'):
                line = line.decode('utf8')
            if line.startswith('>'):
                line = regex.sub('', line.lstrip('>').rstrip())
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
                seq_name = '>' + idx[line] + '\n'
                try:
                    outF.write(seq_name)
                except TypeError:
                    outF.write(seq_name.encode('utf-8'))
            else:
                if hit:
                    try:
                        outF.write(line)
                    except TypeError:
                        outF.write(line.encode('utf-8'))
                        
    logging.info('File written: {}'.format(output))
    logging.info('Number of seqs written: {}'.format(len(found.keys())))
    return found

def idx_overlap(idx1, idx2, verbose=True):
    """
    Getting overlapping keys
    """
    idx = {}
    for x in set(idx1.keys()) & set(idx2.keys()):
        idx[x] = idx1[x]
    if verbose:
        logging.info('No. of seqIDs in idx1: {}'.format(len(idx1.keys())))
        logging.info('No. of seqIDs in idx2: {}'.format(len(idx2.keys())))
        logging.info('No. of overlapping seqIDs: {}'.format(len(idx.keys())))
    return idx

def write_name_idx(idx, tax, genome_id, genome_len):
    """
    Writing gene metadata
    """
    header = ['seq_uuid', 'seq_orig_name', 'domain', 'phylum',
              'class', 'order', 'family', 'genus', 'species',
              'taxid', 'genome_name', 'genome_length_bp']
    print('\t'.join(header))
    for k,v in idx.items():
        print('\t'.join([v, k] + tax + [genome_id, str(genome_len)]))        

def get_genome_length(infile):
    """
    Getting the length of the genome
    """
    if infile.endswith('.gz'):
        _open = lambda x: gzip.open(x, 'rb')
    else:
        _open = lambda x: open(x)
    seq_len = 0
    with _open(infile) as inF:
        for line in inF:
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            if not line.startswith('>'):
                seq_len += len(line.rstrip())
    return seq_len
                
def main(args):
    """
    Main interface
    """
    tax = format_taxonomy(args.taxonomy, args.taxID)
    genome_len = get_genome_length(args.genome_file)
    genomeID = args.accession
    
    # creating the seq header index
    seq_idx = make_index(args.fasta1) & make_index(args.fasta2)
    seq_idx = {x:str(uuid.uuid4()).replace('-', '') for x in seq_idx}
    
    # filtering the fasta files
    idx = filter_fasta(args.fasta1, seq_idx,
                       args.fasta1_output, gzip_out=args.gzip)
    idx = filter_fasta(args.fasta2, seq_idx,
                       args.fasta2_output, gzip_out=args.gzip)
    # creating name index
    write_name_idx(idx, tax, genomeID, genome_len)
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
