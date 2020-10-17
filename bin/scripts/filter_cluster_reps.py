#!/usr/bin/env python
from __future__ import print_function
import sys,os
import argparse
import logging

desc = 'Filtering the cluster reps to just those that lack current annotations'
epi = """DESCRIPTION:
Filtering cluster reps to just those that are lacking any 
annotation data.
Output (filtered fasta) written to STDOUT
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('cluster_membership', metavar='cluster_membership', type=str,
                    help='mmseqs cluster membership file (format: cluster_rep<tab>cluster_member)')
parser.add_argument('query_hits', metavar='query_hits', type=str,
                    help='blast-formatted table of hits (cluster_rep <=> target_db_seqs)')
parser.add_argument('cluster_reps_fasta', metavar='cluster_reps_aa', type=str,
                    help='mmseqs cluster representatives fasta file')
parser.add_argument('--hit-columns', type=str,
                    default='qseqid,sseqid,evalue,pident,alnlen,slen',
                    help='Hit table output columns (default: %(default)s)')      
parser.add_argument('--min-pident', type=float, default=50,
                    help='Min % identity of hit (default: %(default)s)')
parser.add_argument('--min-cov', type=float, default=80,
                    help='Min % alignment coverage of subject sequence length (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def _open(infile):
    if infile.endswith('.gz'):
        return gzip.open(infile, 'rb')
    else:
        return open(infile)

def read_membership(infile):
    """
    Reading in cluster membership table: cluster_rep<tab>cluster_member
    return: 
      dict: {cluster_member : cluster_id} 
    """
    logging.info('Reading file: {}'.format(infile))
    mem = {}
    with _open(infile) as inF:
        for line in inF:
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            line = line.rstrip().split('\t')
            if len(line) < 2:
                continue
            mem[line[1]] = line[0]
    return mem

def read_hits(infile, mem, colnames, min_pident=0, min_cov=0):
    """ Loading query hits
    Return: set(cluster rep)  # clusters with hits
    """
    logging.info('Loading hits table...')    
    clusts = []
    idx = {x:i for i,x in enumerate(colnames.split(','))}
    with _open(infile) as inF:
        for i,line in enumerate(inF):            
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            line = line.rstrip().split('\t')
            if line[0] == '':
                continue
            if len(line) < 2:
                msg = 'Line {}: <2 values in hits table'
                raise ValueError(msg.format(i+1))
            else:
                # filtering
                ## percent identity
                pident = 0
                try:
                    pident = float(line[idx['pident']])
                except KeyError:
                    pass
                if pident < min_pident:
                    continue
                ## coverage of target seq
                cov = 0
                try:
                    cov = float(line[idx['slen']]) / float(line[idx['alnlen']]) * 100
                except KeyError:
                    pass
                if cov < min_cov:
                    continue
                ## getting cluster id
                qseqid = line[idx['qseqid']]
                try:
                    clusts.append(mem[qseqid])
                except KeyError:
                    msg = 'Cannot find "{}" in cluster membership'
                    raise KeyError(msg.format(qseqid))
    clusts = set(clusts)
    logging.info('No of clusters w/ annotations: {}'.format(len(clusts)))
#    print('2ff13be4663246e69aaf365463447d82' in clusts); sys.exit()
    return clusts

def filter_fasta(infile, clust_w_annot):
    cnts = {'all' : 0, 'filtered' : 0, 'kept' : 0}
    to_keep = False
    with _open(infile) as inF:
        for line in inF:            
            if infile.endswith('.gz'):
                line = line.decode('utf-8')
            if line.startswith('>'):
                cnts['all'] += 1
                # already has annotation?
                line = line.lstrip('>').rstrip().split(' ')[0]
                if line in clust_w_annot:
                    to_keep = False
                    cnts['filtered'] += 1
                else:
                    to_keep = True
                    cnts['kept'] += 1
                    print('>' + line, end='')
            elif to_keep == True:
                print(line, end='')
    # status
    logging.info('No. of total seqs: {}'.format(cnts['all']))
    logging.info('No. of filt. seqs: {}'.format(cnts['filtered']))
    logging.info('No. of kept seqs: {}'.format(cnts['kept']))
    
def main(args):
    # determining which clusters already have annotations
    clust_w_annot = read_hits(args.query_hits,
                              read_membership(args.cluster_membership),
                              colnames = args.hit_columns,
                              min_pident = args.min_pident,
                              min_cov = args.min_cov)
    # filtering fasta to just those lacking annotations
    filter_fasta(args.cluster_reps_fasta, clust_w_annot)
    
    
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

